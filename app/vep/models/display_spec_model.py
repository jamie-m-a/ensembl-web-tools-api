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

from typing import Iterator, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

# The value formats the frontend's `formatValue` understands. `text` is the
# default (stringify as-is); the rest are the existing formatter functions.
RowFormat = Literal["text", "num", "humanize", "phenotype", "join"]


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

    @model_validator(mode="after")
    def _exactly_one_source(self) -> "DisplayRow":
        if bool(self.source) == bool(self.compose):
            raise ValueError("row needs exactly one of `from` or `compose`")
        return self

    def field_refs(self) -> list[str]:
        return [self.source] if self.source else self.compose.field_refs()


class DisplayBlock(BaseModel):
    """A run of rows, optionally under the option's own sub-heading.

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
    rows: list[DisplayRow]


class DisplayOptionSpec(BaseModel):
    """How one form option renders: a sequence of blocks.

    A sequence, not a single block, because an option can legitimately emit more
    than one: `eve` is a bare EVE row *plus* a sibling popEVE heading block.
    """

    model_config = ConfigDict(extra="forbid")

    option_id: str
    blocks: list[DisplayBlock]

    def field_refs(self) -> Iterator[tuple[str, str]]:
        """Every `<plugin>.<field>` this option reads, split into its two parts.
        A reference that is not `plugin.field` shaped yields an empty field name,
        which the consistency check reports."""
        for block in self.blocks:
            for row in block.rows:
                for ref in row.field_refs():
                    plugin, _, field = ref.partition(".")
                    yield plugin, field


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
