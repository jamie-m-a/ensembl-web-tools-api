"""Static, strongly-typed model of the *display* spec: how one option's parsed
annotation is laid out in the results annotation detail.

The parsing spec says how a plugin's CSQ columns become structured data; this
says how that data is presented — the labels, order, headings, number formats
and placeholders that were, until now, twelve hand-written `case` bodies in the
frontend's `VepResultsAnnotationDetail`. Moving them here makes the backend the
single owner of the option contract end to end (which options exist, how they
are parsed, how they are shown) and lets the frontend render generically.

It is authored per genome, so unlike the per-job display *panels* it lives
inside the merged spec document as a third sibling section, under the same
content digest, and is pinned to a job for free.

Deliberately small: every field here maps 1:1 onto a rendering primitive the
frontend already has (`RowSpec` / `renderRowGroup` / `renderRowBlock`). Nothing
in this model invents new rendering behaviour, and options whose output is
interactive or derived (ClinVar, OpenTargets, ProtVar, ...) are deliberately
*not* expressible — they stay as frontend overrides.
"""

import re
from typing import Annotated, Iterator, Literal, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator

# The value formats the frontend's `formatValue` understands. `text` is the
# default (stringify as-is); the rest are the existing formatter functions.
# `humanize_join` humanises each element of a list then joins them — ClinVar's
# significance terms, shown as one comma-separated value.
RowFormat = Literal[
    "text", "num", "humanize", "phenotype", "join", "humanize_join"
]

# `{field}` placeholders in a link template — the item fields interpolated into
# the URL (e.g. ".../term/{id}").
_TEMPLATE_FIELD = re.compile(r"\{(\w+)\}")


class ComposeSpec(BaseModel):
    """A row value built from more than one field.

    Only one shape exists today: `with_score`, the frontend's `withScore` —
    "Likely benign (0.07)" from a classification plus its score. AlphaMissense
    and EVE both need it, and both drop the row entirely when the
    *classification* is absent, whatever the score says.
    """

    model_config = ConfigDict(extra="forbid")

    format: Literal["with_score"]
    classification: str
    score: str

    def field_refs(self) -> list[str]:
        return [self.classification, self.score]


class SubOption(BaseModel):
    """The form sub-option a row's value comes from.

    Lets "Show all" list a sub-option that ran but produced nothing as a dash
    (the default view drops the empty row instead). `default` mirrors the form
    default: a sub-option left at a default-on value isn't written to the
    submitted parameters, so the frontend treats "absent" as its default (see
    `subOptionRan`). The id is a form option id — the hand-synced seam with
    `form_panels`, like the top-level `option_id`; not a `plugin.field` ref, so
    the display↔parsing check does not touch it.
    """

    model_config = ConfigDict(extra="forbid")

    id: str
    default: bool = False


class WhenSpec(BaseModel):
    """A condition gating whether a block renders, tested against one field.

    `present` -> render only when the field has content; `empty` -> only when it
    is absent (null / '' / empty list). ClinVar uses it to flip between a bare
    "Clinical significance" row (no conflicting breakdown) and a headed block
    (breakdown present). The field is a `<plugin>.<field>` reference like a row's
    `from`, resolved against the parsing spec at load like the rest.
    """

    model_config = ConfigDict(extra="forbid")

    present: str | None = None
    empty: str | None = None

    @model_validator(mode="after")
    def _exactly_one(self) -> "WhenSpec":
        if bool(self.present) == bool(self.empty):
            raise ValueError("when needs exactly one of `present` or `empty`")
        return self

    @property
    def field_ref(self) -> str:
        # exactly one is set (validated above)
        return self.present or self.empty  # type: ignore[return-value]


class DisplayRow(BaseModel):
    """One label/value row.

    `from` is a `<plugin>.<field>` reference into the *parsing* spec — the
    plugin id and one of its declared target fields. Which entity that plugin is
    read from (allele or transcript consequence) is deliberately not stated
    here: it already lives on the parsing plugin's `scope`, and is derived at
    serve time (see `MergedSpec.plugin_scopes`).
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    # React list key. Optional: absent means "use the row's position", which is
    # stable for these fixed lists.
    key: str | None = None
    label: str
    # `from` is a Python keyword, hence the alias (as in TargetSpec).
    source: str | None = Field(default=None, alias="from")
    compose: ComposeSpec | None = None
    format: RowFormat | None = None
    mono: bool = False
    # What to show when the value is absent. Unset drops the row entirely; set
    # keeps it and shows this (SpliceAI's deltas always read as a set of eight).
    placeholder: str | None = None
    # Help text for a (?) button beside the label. The text is data; the button
    # is a frontend primitive.
    help: str | None = None
    # The sub-option this row's value comes from. Only affects "Show all": a
    # selected-but-empty sub-option shows a dash there; the default view still
    # drops it. Rows without one behave exactly as before.
    sub_option: SubOption | None = None

    @model_validator(mode="after")
    def _exactly_one_source(self) -> "DisplayRow":
        if bool(self.source) == bool(self.compose):
            raise ValueError("row needs exactly one of `from` or `compose`")
        return self

    def field_refs(self) -> list[str]:
        return [self.source] if self.source else self.compose.field_refs()


class LinkSpec(BaseModel):
    """How to turn a cell value into a link.

    `external` -> a plain anchor (`target=_blank`). `template` is a full URL with
    `{field}` placeholders filled from the item's fields (e.g. a GO term or
    MaveDB URN); `builder` names a frontend link builder for URLs that aren't a
    simple template (ProtVar's algorithmic URL). `app_popup` -> an in-app
    "View in" popup, which is always a named `builder` (it needs the job's genome
    and the consequence, not just the annotation field) — e.g. the protein id.
    """

    model_config = ConfigDict(extra="forbid")

    kind: Literal["external", "app_popup"]
    template: str | None = None
    builder: str | None = None

    @model_validator(mode="after")
    def _template_xor_builder(self) -> "LinkSpec":
        if bool(self.template) == bool(self.builder):
            raise ValueError("link needs exactly one of `template` or `builder`")
        if self.kind == "app_popup" and not self.builder:
            raise ValueError("an app_popup link must use a `builder`")
        return self

    def template_fields(self) -> list[str]:
        """The item field names a `template` interpolates; empty for a builder."""
        return _TEMPLATE_FIELD.findall(self.template) if self.template else []


class CellSpec(BaseModel):
    """One cell of a repeated item (see `DisplayListBlock`).

    `from` is a field *of the list element* (not `plugin.field`) — e.g. `score`
    on a MaveDB assay. Omit it for a scalar list whose elements are the value
    themselves (phenotype strings). `link` makes the cell an anchor.
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    label: str | None = None
    # `from` is a Python keyword, hence the alias.
    source: str | None = Field(default=None, alias="from")
    format: RowFormat | None = None
    mono: bool = False
    link: LinkSpec | None = None

    def item_field_refs(self) -> Iterator[str]:
        """Every item field this cell reads: its `from` plus any `{field}`
        placeholders in a link template. Builder links contribute nothing (the
        frontend builder owns its inputs)."""
        if self.source:
            yield self.source
        if self.link:
            yield from self.link.template_fields()


class TruncateSpec(BaseModel):
    """Show the first `visible_count` items with the rest behind a show-more
    toggle (the frontend's `TruncatedList`)."""

    model_config = ConfigDict(extra="forbid")

    visible_count: int = Field(gt=0)


class DisplayItemLabel(BaseModel):
    """The label of a list element rendered as a label/value row (see
    `DisplayItemSpec.label`).

    `from` reads one item field (ClinVar's per-class significance); `template`
    interpolates item fields into text ("Pocket {pocket_id}"). `format` applies
    to a `from` value (e.g. humanize).
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    source: str | None = Field(default=None, alias="from")
    template: str | None = None
    format: RowFormat | None = None

    @model_validator(mode="after")
    def _source_xor_template(self) -> "DisplayItemLabel":
        if bool(self.source) == bool(self.template):
            raise ValueError(
                "item label needs exactly one of `from` or `template`"
            )
        return self

    def item_field_refs(self) -> Iterator[str]:
        """The item fields this label reads: its `from`, or the `{field}`
        placeholders in its `template`."""
        if self.source:
            yield self.source
        if self.template:
            yield from _TEMPLATE_FIELD.findall(self.template)


class DisplayItemSpec(BaseModel):
    """How one element of a list renders.

    Without `label`, a row of one or more inline cells (a GO id + name). With
    `label`, a label/value row instead: `label` is the row's label and the
    `cells` render as its value — ClinVar's per-class counts, ProtVar's pockets.
    """

    model_config = ConfigDict(extra="forbid")

    label: DisplayItemLabel | None = None
    cells: list[CellSpec] = Field(min_length=1)

    def item_field_refs(self) -> Iterator[str]:
        """Every item field this element reads, across its label and cells."""
        if self.label:
            yield from self.label.item_field_refs()
        for cell in self.cells:
            yield from cell.item_field_refs()


class DisplayRowsBlock(BaseModel):
    """A run of fixed rows, optionally under the option's own sub-heading.

    `heading` present -> the frontend's `renderRowBlock` (an `OptionBlock` whose
    heading only appears if a row survived); absent -> `renderRowGroup` (the
    rows on their own).

    `requires` names a plugin that must have produced an annotation at all for
    the block to render. It exists for SpliceAI, whose delta rows carry a
    placeholder: without it, a variant with no SpliceAI annotation would render
    eight dashes instead of nothing (the hand-written case returned early).
    """

    model_config = ConfigDict(extra="forbid")

    kind: Literal["rows"] = "rows"
    heading: str | None = None
    requires: str | None = None
    # A data condition on top of `requires`: render this block only when the
    # named field is present / empty (ClinVar's bare vs headed shapes).
    when: WhenSpec | None = None
    rows: list[DisplayRow]


class DisplayListBlock(BaseModel):
    """A variable-length list: one item (a row of cells) per element of a
    list-valued field, optionally truncated. Covers the options whose output is
    a repeat rather than a fixed set of rows — phenotypes, GO terms, MaveDB
    assays, ...

    `from` is the `<plugin>.<listField>` the elements come from; that field must
    be a parse-plugin target declaring the element's `item_fields`, which the
    cells' `from`/link templates reference.
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    kind: Literal["list"]
    heading: str | None = None
    requires: str | None = None
    when: WhenSpec | None = None
    source: str = Field(alias="from")
    truncate: TruncateSpec | None = None
    item: DisplayItemSpec

    def list_ref(self) -> tuple[str, str]:
        """The `(plugin, listField)` this block iterates."""
        plugin, _, field = self.source.partition(".")
        return plugin, field


class DisplayGroupBlock(BaseModel):
    """A run of sub-blocks under one optional heading, gated as a whole by `when`.

    Lets a heading span more than one block conditionally: ClinVar's conflicting
    case is a "Classification" row plus a per-class breakdown list under one
    "Clinical significance" heading, shown only when the breakdown is present.
    Distinct from `DisplayOptionSpec.heading`, which spans *every* block of the
    option unconditionally.
    """

    model_config = ConfigDict(extra="forbid")

    kind: Literal["group"] = "group"
    heading: str | None = None
    when: WhenSpec | None = None
    blocks: list["DisplayBlock"]


# A block is a fixed set of rows, a repeated list, or a group of sub-blocks,
# discriminated on `kind`.
DisplayBlock = Annotated[
    Union[DisplayRowsBlock, DisplayListBlock, DisplayGroupBlock],
    Field(discriminator="kind"),
]

# `DisplayGroupBlock.blocks` refers to the union defined just above it.
DisplayGroupBlock.model_rebuild()


class DisplayOptionSpec(BaseModel):
    """How one form option renders: a sequence of blocks.

    A sequence, not a single block, because an option can legitimately emit more
    than one: `eve` is a bare EVE row *plus* a sibling popEVE heading block.
    """

    model_config = ConfigDict(extra="forbid")

    option_id: str
    # An option-level heading wrapping *all* the option's blocks in one
    # `OptionBlock`, shown whenever the option renders anything. For an option
    # whose output spans more than one block under a single heading — MaveDB's
    # "Variant" row plus its assays list — where a per-block heading can't reach
    # across blocks. Distinct from a block's own `heading` (use one or the other).
    heading: str | None = None
    blocks: list[DisplayBlock]

    def iter_blocks(self) -> Iterator[DisplayBlock]:
        """Every block in the option, groups flattened. The group block itself is
        yielded (so its own `when` gets checked) as well as its children, so the
        consistency check can treat the tree as a flat list of blocks."""

        def walk(blocks: list[DisplayBlock]) -> Iterator[DisplayBlock]:
            for block in blocks:
                yield block
                if isinstance(block, DisplayGroupBlock):
                    yield from walk(block.blocks)

        yield from walk(self.blocks)


class DisplaySpec(BaseModel):
    """The display half of the merged document: every laid-out option."""

    model_config = ConfigDict(extra="forbid")

    options: list[DisplayOptionSpec]


class DisplayPayload(BaseModel):
    """What the results response carries: the display spec plus the plugin ->
    scope map derived from `parsing`, which the frontend needs to know whether
    to read a row's plugin from the allele or the transcript consequence.

    The scopes are derived rather than authored so there is only ever one place
    that states them (the parsing plugin), and no hand-synced copy to drift.
    """

    model_config = ConfigDict(extra="forbid")

    options: list[DisplayOptionSpec]
    plugin_scopes: dict[str, str]
