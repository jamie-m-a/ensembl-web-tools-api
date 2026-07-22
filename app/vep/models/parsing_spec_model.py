"""Static, strongly-typed model of the CSQ parsing spec.

The spec describes how to turn one plugin's CSQ columns into structured output,
replacing a hand-written `_parse_*` function. It is *data* — the `parsing`
section of the merged JSON under `vep/specs/`, later served by the annotation
API — so it is validated hard on arrival: this model is the contract with that
data.

Deliberately strict (`extra="forbid"`): a spec with an unknown key is a spec we
do not understand, and failing loudly at load time is much cheaper than silently
producing empty annotations at parse time.
"""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

# "raw" is not a coercion: it captures the element's own source text, so a value
# survives verbatim even where the named fields misread it. ProtVar uses this as
# a hedge — its column layout is only known best-effort.
ValueType = Literal["string", "float", "int", "raw"]

# Transforms understood by the interpreter. Kept deliberately small; this set was
# derived by enumerating the existing `_parse_*` functions rather than invented.
Transform = Literal[
    "scalar", "list", "first", "zip", "regex", "pattern_map", "chunk", "positional",
    "key_value",
]


class FieldSpec(BaseModel):
    """One output field of a composite value (e.g. one column of a `zip`, a named
    group of a `regex`, or a slot of a `positional`/`chunk`).

    A `raw`-typed field consumes no positional slot — it reports the source text
    of the element it sits in.
    """

    model_config = ConfigDict(extra="forbid")

    field: str
    type: ValueType = "string"
    # String tidying, applied in order after coercion. VEP escapes spaces as
    # underscores in free text (GO term names, phenotype labels), so undoing
    # that is a general need rather than a GO quirk.
    replace: dict[str, str] | None = None
    strip: bool = False


class DropWhen(BaseModel):
    """When to discard a produced element. Exactly one mode.

    all_null  every field of the element came out null (MaveDB: a position where
              neither score nor urn is real).
    null      the named field came out null (OpenTargets: a row with no disease
              is not an association, whatever else it carries).
    """

    model_config = ConfigDict(extra="forbid")

    all_null: bool = False
    null: str | None = None

    @model_validator(mode="after")
    def _exactly_one_mode(self) -> "DropWhen":
        if bool(self.all_null) == bool(self.null):
            raise ValueError("drop_when needs exactly one of `all_null` or `null`")
        return self


class PostOp(BaseModel):
    """An operation over the whole produced list, applied in order.

    dedup  drop elements identical to an earlier one (the OpenTargets plugin
           currently emits duplicate rows).
    sort   order by `by`. `nulls` places elements whose key is null, and is
           independent of `desc` — "strongest first, unscored last" needs
           desc + nulls: last.
    """

    model_config = ConfigDict(extra="forbid")

    op: Literal["dedup", "sort"]
    by: str | None = None
    desc: bool = False
    nulls: Literal["first", "last"] = "last"

    @model_validator(mode="after")
    def _check_op_shape(self) -> "PostOp":
        if self.op == "sort" and not self.by:
            raise ValueError("sort requires `by`")
        if self.op == "dedup" and self.by:
            raise ValueError("dedup takes no `by`")
        return self


class WhenSpec(BaseModel):
    """A condition on another CSQ column, gating whether a target is built.

    `includes` tests membership of the '&'-split list, not a substring of the raw
    value — ClinVar surfaces its breakdown only when the classification list
    contains exactly "Conflicting_classifications_of_pathogenicity", and a
    substring test would also fire on a value that merely embedded that text.
    """

    model_config = ConfigDict(extra="forbid")

    field: str
    includes: str


class TargetSpec(BaseModel):
    """How to build one output field from one or more CSQ columns.

    `from` names the source column(s); `field` names the output. Transforms:
      scalar       one column -> one value
      list         one column -> '&'-split list, empties and 'NA' dropped
      first        one column -> first real item of a '&'-split list
      zip          N aligned '&'-lists -> list of objects (positions preserved,
                   so 'NA' placeholders still occupy a slot and keep the columns
                   aligned with each other)
      regex        one column -> object(s) from named groups; `each` applies the
                   pattern per '&'-item, otherwise to the whole value. Items
                   that do not match are skipped.
      pattern_map  columns matching `from_pattern` -> dict keyed by the
                   wildcard. The columns are discovered from the CSQ header at
                   runtime, so the field set need not be known up front (this is
                   how gnomAD's per-ancestry AF columns work).
      chunk        one column -> list of objects, taking `size` '&'-items per
                   object (ProtVar's interaction interfaces are partner & score
                   repeating).
      positional   one column -> one object, `as` assigned to '&'-items strictly
                   by index. Items beyond `as` are ignored; missing ones are
                   null. Use `wrap: "list"` where the output is a
                   single-element list.
      key_value    one column -> dict, splitting on `pair_delimiter` then
                   `kv_delimiter`. Order-independent by construction — for a
                   value whose pair order is not meaningful (or, as observed in
                   UTRAnnotator's 5UTR_annotation, not stable), this is the
                   correct read; a plain scalar copies whatever order the
                   plugin happened to emit. A piece without `kv_delimiter` is
                   dropped rather than raising, since malformed/legacy pieces
                   should not break parsing of an otherwise-good value.
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    field: str
    # `from` and `as` are Python keywords, hence the aliases.
    source: str | list[str] | None = Field(default=None, alias="from")
    transform: Transform
    type: ValueType = "string"
    # `zip` / `regex`: the output fields. For zip they match `from` positionally;
    # for regex each `field` names the regex group to read.
    as_fields: list[FieldSpec] | None = Field(default=None, alias="as")
    # `zip` only: whether to iterate to the longest or shortest input column.
    # The existing parsers disagree — MaveDB pads to the longest, OpenTargets
    # truncates to the shortest — so it has to be explicit.
    align: Literal["max", "min"] = "max"
    # `zip` / `chunk`: discard produced elements, then reshape the list.
    drop_when: DropWhen | None = None
    post: list[PostOp] | None = None
    # `regex` only.
    pattern: str | None = None
    each: bool = False
    # `pattern_map` only: a column-name pattern with one `{placeholder}`, e.g.
    # "gnomAD_exomes_AF_{pop}", plus any matching columns to leave out (the
    # overall-AF column can itself match the pattern).
    from_pattern: str | None = None
    exclude: list[str] | None = None
    # `chunk` only: how many '&'-items make up one object.
    size: int | None = None
    # `positional` only: emit the single object inside a list.
    wrap: Literal["list"] | None = None
    # `key_value` only.
    pair_delimiter: str | None = None
    kv_delimiter: str | None = None
    # Build this target only when the condition holds; otherwise it comes out
    # empty (ClinVar's breakdown is only read for conflicting classifications).
    when: WhenSpec | None = None
    # For a target whose value is a list of objects (zip/regex/chunk/...): the
    # keys each element carries. Purely declarative — it does not change parsing;
    # it lets the display spec's `list` blocks reference an element's fields
    # (e.g. a MaveDB assay's `urn`/`score`) and have those refs validated at load
    # time, the list-item analogue of the top-level `field` refs.
    item_fields: list[str] | None = None

    @model_validator(mode="after")
    def _check_transform_shape(self) -> "TargetSpec":
        if self.transform == "zip":
            if not isinstance(self.source, list):
                raise ValueError("zip requires `from` to be a list of columns")
            if not self.as_fields:
                raise ValueError("zip requires `as`")
            if len(self.as_fields) != len(self.source):
                raise ValueError("zip requires one `as` entry per `from` column")
        elif self.transform == "regex":
            if not isinstance(self.source, str):
                raise ValueError("regex requires `from` to be a single column")
            if not self.pattern:
                raise ValueError("regex requires `pattern`")
            if not self.as_fields:
                raise ValueError("regex requires `as` naming the groups to read")
        elif self.transform == "pattern_map":
            if not self.from_pattern:
                raise ValueError("pattern_map requires `from_pattern`")
            if "{" not in self.from_pattern or "}" not in self.from_pattern:
                raise ValueError("pattern_map `from_pattern` needs a {placeholder}")
            if self.source is not None:
                raise ValueError("pattern_map uses `from_pattern`, not `from`")
        elif self.transform == "chunk":
            if not isinstance(self.source, str):
                raise ValueError("chunk requires `from` to be a single column")
            if not self.as_fields:
                raise ValueError("chunk requires `as`")
            if not self.size or self.size < 1:
                raise ValueError("chunk requires a positive `size`")
        elif self.transform == "positional":
            if not isinstance(self.source, str):
                raise ValueError("positional requires `from` to be a single column")
            if not self.as_fields:
                raise ValueError("positional requires `as`")
        elif self.transform == "key_value":
            if not isinstance(self.source, str):
                raise ValueError("key_value requires `from` to be a single column")
            if not self.pair_delimiter or not self.kv_delimiter:
                raise ValueError("key_value requires `pair_delimiter` and `kv_delimiter`")
        else:
            if not isinstance(self.source, str):
                raise ValueError(f"{self.transform} requires `from` to be a single column")
            if self.as_fields:
                raise ValueError(
                    f"`as` is only valid for zip/regex/chunk/positional, not {self.transform}"
                )
        return self


class PluginSpec(BaseModel):
    """How to parse one plugin's contribution to a CSQ entry.

    Two independent "nothing here" rules, mirroring the hand-written parsers:
      csq_fields        which columns this plugin owns. If none are in the CSQ
                        header, the plugin did not run — skip it entirely.
      require_any_input the columns are present, but this record has no value in
                        any of them -> no annotation. Note this tests raw
                        presence, so a literal 'NA' counts as present (matching
                        the current parsers).
      require_any_output built the output, but the fields that carry the payload
                        came out empty -> no annotation.
    """

    model_config = ConfigDict(extra="forbid")

    plugin: str
    scope: Literal["allele", "transcript"]
    # Where the result attaches on the response model, e.g. "mavedb".
    output: str
    csq_fields: list[str]
    require_any_input: list[str] | None = None
    require_any_output: list[str] | None = None
    targets: list[TargetSpec]


class ParsingSpec(BaseModel):
    """A whole parsing-spec document: every plugin, for one genome."""

    model_config = ConfigDict(extra="forbid")

    # Content digest of this document; pins a job to the ruleset that produced
    # its options (see the sidecar written at submission). Optional so a
    # ParsingSpec can nest inside the merged document, which owns the single
    # digest and stamps this to match (merged_spec_model.py); still computed by
    # spec_loader when a ParsingSpec is the whole loaded document.
    spec_version: str = ""
    genome: dict | None = None
    plugins: list[PluginSpec]

    def plugin(self, name: str) -> PluginSpec | None:
        """The spec for one plugin by name, or None."""
        return next((p for p in self.plugins if p.plugin == name), None)
