"""Help that belongs to the input form and not to the results.

An option's help is hung on whichever node turns out to be its visible title in
the results. For an option whose rows carry help of their own, that is a row
which already has a (?) button — two of them, saying overlapping things.
`form_only` keeps the help on the form and drops it from the results.
"""

import json
from types import SimpleNamespace

from app.vep.models.display_panels_model import (
    DisplayPanel,
    dump_display_panels,
    to_display_panels,
)
from app.vep.models.option_help_model import HelpSpec, OptionHelp
from app.vep.utils.vcf_results import _drop_form_only_help

HELP = HelpSpec(
    options=[
        OptionHelp(option_id="flagged", description="On the form only.", form_only=True),
        OptionHelp(option_id="plain", description="Everywhere."),
    ]
)
SPEC = SimpleNamespace(help=HELP)


def _panels() -> list[DisplayPanel]:
    return to_display_panels(
        [
            {
                "id": "a_panel",
                "label": "A panel",
                "options": [
                    {"id": "flagged", "label": "Flagged", "help": {"description": "On the form only."}},
                    {"id": "plain", "label": "Plain", "help": {"description": "Everywhere."}},
                ],
            }
        ]
    )


def test_the_flag_is_not_part_of_the_help_payload():
    """It says where the help goes, not what it says."""
    payload = HELP.payload_for("flagged", "Flagged")
    assert payload["description"] == "On the form only."
    assert "form_only" not in json.dumps(payload)


def test_the_flagged_options_are_nameable():
    assert HELP.form_only_options() == {"flagged"}


def test_the_results_panels_lose_only_the_flagged_option_s_help():
    served = {
        o["id"]: o for p in dump_display_panels(_drop_form_only_help(_panels(), SPEC))
        for o in p["options"]
    }
    assert "help" not in served["flagged"]
    assert served["plain"]["help"] == {"description": "Everywhere."}


def test_the_option_itself_survives():
    """Only its help goes — the option still has to render."""
    served = {
        o["id"]: o for p in dump_display_panels(_drop_form_only_help(_panels(), SPEC))
        for o in p["options"]
    }
    assert served["flagged"]["label"] == "Flagged"


def test_a_spec_flagging_nothing_is_left_alone():
    spec = SimpleNamespace(help=HelpSpec(options=[OptionHelp(option_id="plain", description="x")]))
    before = dump_display_panels(_panels())
    assert dump_display_panels(_drop_form_only_help(_panels(), spec)) == before


def test_no_panels_or_no_spec_is_not_an_error():
    """Jobs pinned before either existed still render."""
    assert _drop_form_only_help(None, SPEC) is None
    before = dump_display_panels(_panels())
    assert dump_display_panels(_drop_form_only_help(_panels(), None)) == before
