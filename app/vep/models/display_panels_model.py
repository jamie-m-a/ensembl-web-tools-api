"""The option panels pinned to a job at submission, so its results render
against the options it was actually submitted with.

`form_panels.get_visible_panels` builds these as plain dicts for the input form.
The same structure is pinned per job (a sidecar beside the merged spec) and
handed back on the results response, so the results view can lay itself out from
the submitted panels rather than whatever the live form config says now.

The models are deliberately *permissive*: the panel structure is still evolving
(options carry `category`, `sub_options`, nested `{"type": "group", ...}` nodes,
and select choices which are bare `{label, value}` pairs), and pinning must be
lossless. So every field is optional, unknown keys are kept (`extra="allow"`),
and serialisation drops the keys that were not present in the source dict. The
result round-trips `get_visible_panels()` output exactly.
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, model_serializer


class DisplayOption(BaseModel):
    """One node of a panel's option tree.

    Covers all of the shapes `form_panels` emits: a plain option, an option with
    `sub_options`, a `{"type": "group", ...}` heading whose `options` are nested
    options, and a select's `{"label", "value"}` choices (whose `value` is kept
    as an extra field). Nothing is required, so no shape is rejected.
    """

    model_config = ConfigDict(extra="allow")

    id: str | None = None
    label: str | None = None
    type: str | None = None
    default: Any = None
    category: str | None = None
    sub_options: list[DisplayOption] | None = None
    options: list[DisplayOption] | None = None

    @model_serializer(mode="wrap")
    def _drop_absent(self, handler) -> dict:
        """Emit only the keys the source dict actually had. Absent optional
        fields are None, and no panel value is legitimately None, so dropping
        them makes the dump equal to the original dict."""
        return {key: value for key, value in handler(self).items() if value is not None}


class DisplayPanel(BaseModel):
    """A single option panel (e.g. "Allele frequencies") and its options."""

    model_config = ConfigDict(extra="allow")

    id: str
    label: str
    options: list[DisplayOption] = []


def to_display_panels(panels: list[dict]) -> list[DisplayPanel]:
    """Validate raw `get_visible_panels()` output into the pinned model."""
    return [DisplayPanel.model_validate(panel) for panel in panels]


def dump_display_panels(panels: list[DisplayPanel]) -> list[dict]:
    """The panels back as plain dicts, identical to what went in."""
    return [panel.model_dump() for panel in panels]
