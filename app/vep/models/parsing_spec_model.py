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
    "key_value", "records", "stack",
]


class DropEntries(BaseModel):
    """Drop the entries of a packed field that match a pattern.

    A field whose value is several entries joined by a separator (IntAct packs
    its interaction participants as `uniprotkb:P42336_and_ensembl:ENSP...`).
    Where `replace` rewrites fixed text, this removes whole entries by what they
    look like, and is the only field-level operation that needs a regex — an
    identifier is never the same twice, so a literal cannot name it.

    If every entry goes, so does the field: an empty packed string is no value,
    and the column drops as any absent field does.
    """

    model_config = ConfigDict(extra="forbid")

    # The separator the entries are joined by, as it appears in the CSQ.
    sep: str
    # A Python regular expression, matched against each entry (`re.search`).
    matching: str
    # Why these entries are dropped, in prose. Carried by the model rather than
    # left to a comment because the spec is JSON and cannot hold one, and a rule
    # that silently removes data is exactly the kind that needs to say why — and
    # to say whether it is meant to be permanent.
    why: str | None = None


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
    # Extra tokens meaning "no value" for this field, on top of '' and 'NA'.
    # ClinVar writes '.' where a condition has no ontology ids, and a VCF-derived
    # source generally uses '.' for an absent value — but '.' is a real value
    # elsewhere, so this is declared per field rather than made global.
    null_values: list[str] | None = None
    # See DropEntries. Applied last, after the tidying above.
    drop_entries: DropEntries | None = None


class Match(BaseModel):
    """One equality on a field of a produced element.

    The right-hand side is either a literal (`equals`) or the value of a CSQ
    column (`equals_column`). It is the same predicate either way, so it is
    compared the same way either way -- which it had not been: this was two
    models with two comparators, and they disagreed. `only_if` compared the raw
    values with `!=` while a join's `where` stringified first, so `equals: "1"`
    against an int-typed field worked in one and was silently inert in the
    other. See `_same` in spec_interpreter for the comparator they now share.
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    field: str
    equals: str | None = None
    # `column` was this field's whole name back when it was `ColumnMatch`, and a
    # spec pinned to a job before the merge still spells it that way. Accepted
    # as an alias so those keep loading -- see the sidecar compatibility tests.
    equals_column: str | None = Field(default=None, alias="column")
    # Applied to the *column's* value before comparing, taking the comparable
    # part from a `key` group -- the same device `RowScope.item_pattern` and a
    # join's `right_key_pattern` use, and needed for the same reason.
    #
    # VEP writes a versioned gene id in `Gene` (`ENSG00000121879.8`) while a
    # plugin's own payload carries the bare accession (`ENSG00000121879`). Those
    # are the same gene and must compare equal; without this they never do, and
    # the failure is silent in the worst direction -- `_same` returns False, so
    # an `unless_matches` rule would drop every element rather than none.
    column_pattern: str | None = None

    @model_validator(mode="after")
    def _one_right_hand_side(self) -> "Match":
        if (self.equals is None) == (self.equals_column is None):
            raise ValueError(
                "a match compares against either a literal (`equals`) or a CSQ "
                f"column (`equals_column`), not both or neither; field "
                f"{self.field!r}"
            )
        if self.column_pattern is not None:
            if self.equals_column is None:
                raise ValueError(
                    "`column_pattern` extracts part of a CSQ column's value, so "
                    f"it needs `equals_column`; field {self.field!r}"
                )
            if "(?P<key>" not in self.column_pattern:
                raise ValueError(
                    "`column_pattern` needs a `key` group naming the part to "
                    f"compare; field {self.field!r}"
                )
        return self


class DropWhen(BaseModel):
    """When to discard a produced element. Exactly one mode.

    all_null       every field of the element came out null (MaveDB: a position
                   where neither score nor urn is real).
    null           the named field came out null (OpenTargets: a row with no
                   disease is not an association, whatever else it carries).
    unless_matches the named field does not equal the given CSQ column. VEP
                   repeats a whole plugin value on every alt allele's CSQ row, so
                   a value that is really per-allele has to be narrowed to the
                   row's own allele: a phenotype association is kept only where
                   its risk allele is this allele. An element whose field is null
                   never matches, so this also drops associations carrying no
                   risk allele at all.

    `only_if` narrows any mode to elements matching a condition on one of their
    own fields, for a requirement that holds for some kinds of element but not
    others — the allele rule applies to a "Variation" association, while a "Gene"
    one legitimately carries no allele and must survive.
    """

    model_config = ConfigDict(extra="forbid")

    all_null: bool = False
    null: str | None = None
    unless_matches: Match | None = None
    only_if: Match | None = None

    @model_validator(mode="after")
    def _exactly_one_mode(self) -> "DropWhen":
        modes = [bool(self.all_null), self.null is not None,
                 self.unless_matches is not None]
        if sum(modes) != 1:
            raise ValueError(
                "drop_when needs exactly one of `all_null`, `null` or "
                "`unless_matches`"
            )
        return self


class PostOp(BaseModel):
    """An operation over the whole produced list, applied in order.

    dedup   drop elements identical to an earlier one (the OpenTargets plugin
            currently emits duplicate rows).
    sort    order by `by`. `nulls` places elements whose key is null, and is
            independent of `desc` — "strongest first, unscored last" needs
            desc + nulls: last.
    exclude drop elements whose `by` field is one of `values` — placeholder terms
            a source emits in place of real data ("ClinVar:_phenotype_not
            _specified" is not a phenotype). Compared case-insensitively against
            the parsed value, so the spec need not mirror a source's casing.

    lookup  add a field derived from another via a shipped reference table —
            GO's aspect, which the plugin's output does not carry (it emits an
            id and a name only). `by` is the field to look up, `into` the field
            to write, `table` the table's name in `vep/data/`. An id the table
            does not know writes null rather than failing: GO releases and the
            annotation files move independently, so an unknown term is a normal
            state, not a broken spec.

    split_field  the inverse of `concat`: cut a field at the first `sep` and
            write one side of it to `into`, keeping the original intact. A
            source that packs two identifiers into one column needs both —
            MaveDB's accession is `urn:mavedb:00000045-a-1#2010`, the score set
            and the variant within it, and the link wants the score set in its
            path and the whole accession in its query. `by` is the field to cut,
            `keep` which side ("before" by default). A value with no `sep` in it
            is all "before" and no "after", so the missing half writes null and
            whatever depends on it drops.

    concat  join two or more fields into a new one — OpenTargets publishes a
            p-value as a mantissa and an exponent in separate columns, and
            `3.32` + `e` + `-28` is the number a reader wants. `fields` names
            the parts in order, `sep` goes between them, `into` is the field to
            write. Any part being null makes the result null rather than a
            malformed string: a row with no p-value must show nothing, not
            "e" on its own.

    curie_link  turn a CURIE list into one link. A source that names a thing in
            several ontologies at once (ClinVar's CLNDISDB —
            `MeSH:D030342,MedGen:C0950123`) has no single id; `prefer` picks
            which authority to trust, in order, falling back to whatever is
            there. `templates` maps the source prefix to its URL, `{id}` being
            the bare accession. `by` is the CURIE-list field, `into` the URL
            field to write, `label_into` (optional) the chosen CURIE itself.
            A list with nothing matching any template writes null, so a
            condition ClinVar gives no usable id for renders as plain text
            rather than a dead link.

    mapped_link  build a link whose URL shape depends on which source an element
            came from. `by` names the field that decides (for phenotypes, the
            `source`), `templates` maps each of its values to a URL, and `into`
            is the URL field to write. Unlike `curie_link` the key is a plain
            field value rather than a prefix inside a CURIE, and the template is
            filled from the element's *own* fields — `{id}` and `{external_id}`
            are both available, which is the point: the Phenotypes plugin gives
            every association both, and which one addresses the record differs
            by source. G2P and the GWAS catalogue are keyed by `id` (an Ensembl
            gene id, an rs id), while OMIM and Orphanet are keyed by
            `external_id`, their own accession, with `id` holding the gene the
            association hangs off.

            A source absent from `templates`, or one whose template names a
            field this element left empty, writes null — so an association from
            a source with no known URL renders as plain text rather than a dead
            link, exactly as `curie_link` does.

    collapse  merge elements differing *only* in the named fields, gathering
            those fields into a nested list. ClinVar files one submission
            against several conditions at once, so a variant's table can carry
            five rows that are one classification by one submitter under five
            disease names — identical in every column but the condition, and
            read as five findings when they are one. `fields` names what may
            differ, `into` the field the gathered sets land in. First-seen
            order, and a group of one collapses too, so the display reads one
            shape rather than two.

    derive_if_empty  build a list from a packed field, but only where the list
            is still empty. ClinVar's conditions table takes its names from
            CLNDN, and a handful of RCV records name a condition CLNDN does not
            carry at all — the record itself has it, as
            `MedGen:C0338106:Colon_adenocarcinoma`, so the cell need not be
            blank. `by` is the packed field, `sep` what separates its entries,
            `pattern` a regex whose named groups become each element's fields,
            and `into` the list to fill. Entries the pattern does not match are
            skipped.

    default  fill a field where the source left it empty. ClinVar accepts a
            submission with no classification at all ("no classification
            provided"), and grouping by classification would otherwise have
            nothing to file it under and drop it — so it becomes an explicit
            "none" that the reader can see, hung off its own condition and RCV
            like any other. `by` is the field, `value` what to put there.

    only_if_differs  keep a nested element's field only where the row does not
            already show that value, writing it to `into` (null otherwise, so
            the display drops it). ClinVar's submitters file against their own
            wording — "Renal tubular dysgenesis of genetic origin" where the
            aggregate says "Renal tubular dysgenesis" — and that is worth
            showing exactly when it is not what the reader can already see.
            `in` is the dotted path to the nested elements, `field` what to
            compare, `against` the dotted path to the values the row shows.
            Elements are copied, never mutated: a join attaches the same object
            to every row it matched.

    `exclude` is a post-op rather than a `drop_when` mode because `drop_when`
    takes exactly one mode and the phenotype targets already spend theirs on the
    per-allele rule.
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    op: Literal[
        "dedup", "sort", "exclude", "lookup", "concat", "curie_link", "collapse",
        "only_if_differs", "default", "derive_if_empty", "mapped_link",
        "split_field",
    ]
    by: str | None = None
    desc: bool = False
    nulls: Literal["first", "last"] = "last"
    values: list[str] | None = None
    # `lookup` and `concat`.
    into: str | None = None
    # `lookup` only.
    table: str | None = None
    # `concat`, and `collapse` (there the fields that may differ).
    fields: list[str] | None = None
    # `default` only: what to put where the source left the field empty.
    value: str | None = None
    # `split_field` only: which side of the separator to write.
    keep: Literal["before", "after"] = "before"
    # `derive_if_empty` only: the regex whose named groups become the fields.
    pattern: str | None = None
    # `only_if_differs` only.
    in_: str | None = Field(default=None, alias="in")
    field: str | None = None
    against: str | None = None
    sep: str = ""
    # `curie_link` only.
    prefer: list[str] | None = None
    templates: dict[str, str] | None = None
    label_into: str | None = None
    # The separator between CURIEs in the source list.
    curie_sep: str = ","

    @model_validator(mode="after")
    def _check_op_shape(self) -> "PostOp":
        if self.op == "sort" and not self.by:
            raise ValueError("sort requires `by`")
        if self.op == "dedup" and self.by:
            raise ValueError("dedup takes no `by`")
        if self.op == "exclude" and not (self.by and self.values):
            raise ValueError("exclude requires `by` and `values`")
        if self.op != "exclude" and self.values is not None:
            raise ValueError("`values` belongs to exclude")
        if self.op == "lookup" and not (self.by and self.into and self.table):
            raise ValueError("lookup requires `by`, `into` and `table`")
        if self.op == "curie_link" and not (self.by and self.into and self.templates):
            raise ValueError("curie_link requires `by`, `into` and `templates`")
        if self.op == "mapped_link" and not (self.by and self.into and self.templates):
            raise ValueError("mapped_link requires `by`, `into` and `templates`")
        if self.op not in ("curie_link", "mapped_link") and self.templates:
            raise ValueError("`templates` belongs to curie_link or mapped_link")
        if self.op != "curie_link" and self.prefer:
            raise ValueError("`prefer` belongs to curie_link")
        if self.op != "lookup" and self.table:
            raise ValueError("`table` belongs to lookup")
        if self.op == "concat" and not (self.fields and self.into):
            raise ValueError("concat requires `fields` and `into`")
        if self.op == "concat" and len(self.fields) < 2:
            raise ValueError("concat needs at least two `fields`")
        if self.op == "derive_if_empty" and not (
            self.by and self.into and self.pattern
        ):
            raise ValueError(
                "derive_if_empty requires `by`, `into` and `pattern`"
            )
        if self.op != "derive_if_empty" and self.pattern:
            raise ValueError("`pattern` belongs to derive_if_empty")
        if self.op == "split_field" and not (self.by and self.into and self.sep):
            raise ValueError("split_field requires `by`, `into` and `sep`")
        if self.op == "default" and not (self.by and self.value):
            raise ValueError("default requires `by` and `value`")
        if self.op != "default" and self.value is not None:
            raise ValueError("`value` belongs to default")
        if self.op == "only_if_differs" and not (
            self.in_ and self.field and self.against and self.into
        ):
            raise ValueError(
                "only_if_differs requires `in`, `field`, `against` and `into`"
            )
        if self.op != "only_if_differs" and (self.in_ or self.field or self.against):
            raise ValueError("`in`/`field`/`against` belong to only_if_differs")
        if self.op == "collapse" and not (self.fields and self.into):
            raise ValueError("collapse requires `fields` and `into`")
        if self.op not in ("concat", "collapse") and self.fields is not None:
            raise ValueError("`fields` belongs to concat or collapse")
        if self.op not in (
            "lookup", "concat", "curie_link", "collapse", "only_if_differs",
            "derive_if_empty", "mapped_link", "split_field",
        ) and self.into:
            raise ValueError(
                "`into` belongs to lookup, concat, curie_link, mapped_link, "
                "split_field, collapse or only_if_differs"
            )
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
    # The separator and key-extraction `RowScope` grew and this did not, though
    # both are the same membership test. Without them a `when` could only be
    # written against an '&'-separated column of bare values -- the enriched
    # ClinVar columns are '+'-separated, and some carry decorated entries.
    sep: str = "&"
    item_pattern: str | None = None


class StackGroup(BaseModel):
    """One source group of a `stack`: a `zip` over its own columns, tagged.

    `const` is what makes the stacked list usable: the rows of different groups
    are the same shape but not the same thing, and the tag is the only record of
    which group a row came from. ClinVar states the same facts three times over,
    once per classification type, in three sets of columns that carry the type
    only in their *names* — CLNDN vs ONCDN vs SCIDN. Tagging on the way in turns
    that back into data, so one list can be filtered, joined and split by type
    rather than three lists having to be kept in step.
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    # `from` and `as` are Python keywords, hence the aliases.
    source: list[str] = Field(alias="from")
    as_fields: list[FieldSpec] = Field(alias="as")
    # Fields given the same value on every row this group produces.
    const: dict[str, str] = Field(default_factory=dict)
    sep: str = "&"
    # Read each column whole rather than splitting it, so the group contributes
    # exactly one row. For a group of scalar columns whose values may themselves
    # contain the separator: ClinVar's aggregate germline classification is one
    # assertion even when it reads `Conflicting_classifications_of_pathogenicity
    # +risk_factor`, and splitting it would invent a second classification with
    # no review status of its own.
    split: bool = True
    align: Literal["max", "min"] = "max"

    @model_validator(mode="after")
    def _one_as_per_column(self) -> "StackGroup":
        if len(self.as_fields) != len(self.source):
            raise ValueError("a stack group needs one `as` entry per `from` column")
        return self


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
                   repeating). `record_sep` marks the boundary between objects
                   where the source states one; without it the items are one
                   flat run, cut every `size`.
      positional   one column -> one object, `as` assigned to '&'-items strictly
                   by index. Items beyond `as` are ignored; missing ones are
                   null. Use `wrap: "list"` where the output is a
                   single-element list.
      records      one column -> list of objects, two levels of separator:
                   `sep` between records, `item_sep` between a record's fields,
                   which are then assigned to `as` by index (like `positional`,
                   but repeating). This is the shape of a source that packs whole
                   sub-records into one column — ClinVar's per-submitter (15
                   fields) and per-RCV (5 fields) data.
      stack        several groups of columns -> one list. Each group in `of` is
                   a `zip` over its own columns, tagged with that group's
                   `const` fields, and the groups' rows are concatenated in
                   order. For a source that publishes the same shape several
                   times over in differently-named columns (see StackGroup).
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
    # The separator between items, for the transforms that split a column
    # (`list`, `first`, `zip`, `chunk`, `positional`, and `regex` with `each`).
    # '&' is what VEP writes — it rewrites both ',' and '|' to '&' — so a source
    # carrying structure below that level must use a delimiter VEP leaves alone:
    # the enriched ClinVar VCF uses '~' between subfields and '+' between
    # repeats. Defaults to '&', so every existing target is unaffected.
    sep: str = "&"
    # Percent-decode the produced value's string leaves. Off by default: only a
    # source that escapes its own separators needs it (the enriched ClinVar VCF
    # escapes '% , ; = | & ~ +' inside values). Applied *after* every split, so
    # an encoded '%2C' can never be read as a delimiter.
    decode: bool = False
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
    # `stack` only: the source groups, concatenated in order.
    of: list[StackGroup] | None = None
    # `chunk` only: how many '&'-items make up one object.
    size: int | None = None
    # `chunk` only: a separator *between* objects, above `sep`.
    #
    # Without it a chunk target reads one flat run of items and cuts it every
    # `size`, which is right only while every object is exactly that long — a
    # source that ever writes a short one silently shifts every later object's
    # fields along by the difference. With it, each record is chunked on its own,
    # so a malformed one damages only itself.
    #
    # Both are read, deliberately: ProtVar moved from the flat run to '+' between
    # pockets, and a job submitted before that change is still being read for a
    # week afterwards (see the class docstring on VEP's separator rewriting —
    # '+' is one of the few characters it leaves alone).
    record_sep: str | None = None
    # `positional` only: emit the single object inside a list.
    wrap: Literal["list"] | None = None
    # `key_value` only.
    # `records` only: the separator *within* one record, below `sep`.
    item_sep: str = "~"
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
    # Built for the joins to draw from, and dropped once they have run.
    #
    # ClinVar's submissions and RCV records are parsed as their own lists so a
    # join can file each one under the condition it belongs to. Nothing displays
    # them at that level -- they are read through the conditions -- and
    # `_apply_joins` attaches the very same objects, so leaving them in place
    # shipped every submission twice. That was 40% of ClinVar's payload.
    join_source: bool = False

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
        elif self.transform == "records":
            if not isinstance(self.source, str):
                raise ValueError("records requires `from` to be a single column")
            if not self.as_fields:
                raise ValueError("records requires `as` naming each record's fields")
        elif self.transform == "stack":
            if not self.of:
                raise ValueError("stack requires `of` naming its source groups")
            if self.source is not None:
                raise ValueError("stack reads its columns from `of`, not `from`")
        else:
            if not isinstance(self.source, str):
                raise ValueError(f"{self.transform} requires `from` to be a single column")
            if self.as_fields:
                raise ValueError(
                    f"`as` is only valid for zip/regex/chunk/positional, not {self.transform}"
                )
        # Checked outside the chain above, which each branch leaves early.
        if self.record_sep is not None:
            if self.transform != "chunk":
                raise ValueError(
                    f"`record_sep` is only valid for chunk, not {self.transform}"
                )
            if self.record_sep == self.sep:
                raise ValueError("chunk `record_sep` must differ from `sep`")
        return self


class JoinSpec(BaseModel):
    """Merge one of a plugin's produced lists into another, matching on a key.

    Every transform reads a single column, so a source that spreads one logical
    table across several columns cannot be assembled by them alone — ClinVar
    names a condition in one column, its per-submitter classifications in
    another, and its RCV records in a third, all keyed by the condition's name.
    This runs after the targets are built and stitches them together.

    `right_key_pattern` is what makes it general rather than a ClinVar special
    case: two lists routinely key on the same thing while one writes it
    decorated (ClinVar's RCV condition is `MedGen:C4540192:<name>` where the
    condition list has the bare `<name>`). A regex with a `key` group extracts
    the comparable part; without one the value is used as it stands.

    `count_by` summarises the matches instead of attaching them: grouped by that
    field, in first-seen order, as `[{<count_by>, count}]`. That is the usual
    shape for "how many submitters said what", and keeps the counting out of the
    display layer. `nest_as` additionally keeps each group's own members under
    it, so the summary and the rows behind it stay together — a count the reader
    can open is a count plus its evidence, and pairing them here is what stops
    the display having to re-derive which rows a count was made of.
    """

    # populate_by_name so a serialised spec round-trips: dumping writes the field
    # names (`source`/`as_field`), the document uses the aliases (`from`/`as`),
    # and a pinned sidecar has to load back either way.
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    # The list to enrich, and the list to draw from (both `field` names of this
    # plugin's targets).
    into: str
    # `from` is a Python keyword.
    source: str = Field(alias="from")
    left_key: str
    right_key: str
    # Applied to the right-hand key before comparing; needs a `key` group.
    right_key_pattern: str | None = None
    # The same two, for the *left* key. A left row can be about several things
    # at once: one ClinVar RCV record covers up to five conditions, listed in
    # one '+'-joined field as `MedGen:C0266313:Renal_tubular_dysgenesis`, so the
    # separator splits them and the pattern takes the name off each.
    left_key_sep: str | None = None
    left_key_pattern: str | None = None
    # When one right-hand row belongs under several left-hand rows, the
    # separator its key list uses. ClinVar files one submission (or one RCV
    # record) against several conditions at once, '+'-joined; the row then
    # appears under each of them. Without this it would match none of them.
    right_key_sep: str | None = None
    # ClinVar's submitters write the same condition in different cases.
    case_insensitive: bool = False
    # Consider only right-hand rows holding this value: the count of what a
    # source itself says counts, rather than a count we infer. ClinVar flags
    # which submissions produced the aggregate classification, and no rule we
    # could write over the terms would agree with it — an expert-panel review
    # makes one submission the aggregate and the other 43 not.
    where: Match | None = None
    # Further equalities a match must satisfy, as {left field: right field}.
    #
    # A single key is not always enough to identify a row: ClinVar lists the
    # same condition under more than one classification type (Rosette-forming
    # glioneuronal tumor is both a germline and a somatic one), so joining on
    # the name alone files a somatic submission under the germline condition.
    # Compared like the key, so `case_insensitive` applies here too.
    also_match: dict[str, str] | None = None
    # The field added to each left row. Unused by a `count_into` join, whose
    # product is the number rather than the rows.
    as_field: str | None = Field(default=None, alias="as")
    # Attach how many rows matched rather than the rows themselves. For a count
    # the reader sees as a number ("39 of 44"), where a list would only be
    # counted again at render time.
    count_into: str | None = None
    # Summarise rather than attach: group the matches by this field and count.
    count_by: str | None = None
    # With `count_by`, the field each group carries its own members under.
    nest_as: str | None = None
    @model_validator(mode="after")
    def _nesting_needs_something_to_nest_under(self) -> "JoinSpec":
        if self.nest_as and not self.count_by:
            raise ValueError(
                "`nest_as` names the field a *group's* members hang off, so it "
                f"needs `count_by` to group by; join into {self.into!r} has none"
            )
        return self

    @model_validator(mode="after")
    def _writes_exactly_one_field(self) -> "JoinSpec":
        if bool(self.as_field) == bool(self.count_into):
            raise ValueError(
                "a join writes either the matched rows (`as`) or how many there "
                f"were (`count_into`), not both or neither; join into {self.into!r}"
            )
        return self

    # Deleted as a concept -- it was write-only, and both declarations had
    # drifted from the targets they described (what a joined-in row carries is
    # derived now, see `MergedSpec._list_element_fields`). Still *accepted*,
    # because a spec pinned to a job before that change carries it, and
    # `extra="forbid"` would otherwise reject the whole sidecar and leave that
    # job with no annotations at all. Never read.
    item_fields: list[str] | None = None

    def produced_fields(self) -> list[str]:
        """The fields this join adds to each row of `into`."""
        return [field for field in (self.as_field, self.count_into) if field]


class RowScope(BaseModel):
    """Which of a variant's CSQ rows a plugin's annotation actually belongs to.

    VEP repeats a custom's columns on *every* CSQ row of the variant, so an
    annotation about one gene is served against every gene the variant touches.
    ClinVar's record for 22:23834143 is about SMARCB1, but DERL3's transcripts
    overlap the same position, and the classification was appearing under both.

    `column` is the row's own (SYMBOL); `listed_in` is the plugin's column naming
    what it is about (ClinVar_GENEINFO, `SMARCB1:6598&WARS2-AS1:101929147`).
    `item_pattern` takes the comparable part of an entry via a `key` group, as a
    join's `right_key_pattern` does.

    Narrowing applies only when there is something to narrow by: a row with no
    value of its own, or an annotation naming nothing, is left alone. Dropping
    those would trade a wrong attribution for a missing one — an intergenic row
    has no symbol to match, and the annotation is still true of the variant.
    """

    model_config = ConfigDict(extra="forbid")

    column: str
    listed_in: str
    sep: str = "&"
    item_pattern: str | None = None


class JoinedPostOp(PostOp):
    """A post-op over a target, applied *after* the joins have run.

    A target's own `post` cannot see what a join added, because the joins stitch
    the targets together afterwards. Ordering a list by a joined-in value needs
    this later pass: whether a ClinVar condition has a submission behind the
    aggregate classification is only known once the submissions have been
    matched to it.
    """

    target: str


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
    # Fields whose value is a URL built from the variant rather than read from a
    # CSQ column, as `{field: template}`. The tokens are `chromosome`,
    # `position`, `reference` and `alternative`; a template naming anything else,
    # or a variant missing a part, yields no field.
    #
    # For a resource that keys its pages on the variant when the plugin emits no
    # URL of its own. The result is an ordinary parsed field, so a display row
    # points at it with `link_from` exactly as it does at Geno2MP's own
    # `Geno2MP_URL` column.
    variant_links: dict[str, str] | None = None
    # Where the result attaches on the response model, e.g. "mavedb".
    output: str
    csq_fields: list[str]
    # Restrict this plugin to the CSQ rows it is really about (see RowScope).
    applies_to: RowScope | None = None
    require_any_input: list[str] | None = None
    require_any_output: list[str] | None = None
    targets: list[TargetSpec]
    # Applied after every target is built (see JoinSpec).
    joins: list[JoinSpec] | None = None
    # Applied to a named target after the joins (see JoinedPostOp).
    post_joins: list[JoinedPostOp] | None = None


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
