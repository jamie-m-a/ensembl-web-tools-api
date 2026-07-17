"""Differential tests: spec-driven parsing vs the hand-written `_parse_*` bank.

The hand-written parsers are the oracle. For the same CSQ entry, the interpreter
driven by `parsing_specs/human_grch38.json` must produce exactly what the
corresponding `_parse_*` produces (compared as plain data, via model_dump).

This is what proves the spec vocabulary is sufficient before anything is
rewired, so the fixtures are deliberately shared with test_csq_parsers.
"""

import pytest

from app.tests.test_csq_parsers import EMPTY, INDEX_MAP, row_list
from app.vep.models.parsing_spec_model import ParsingSpec, TargetSpec
from app.vep.utils.spec_interpreter import apply_plugin_spec
from app.vep.utils.spec_loader import SPEC_DIR, load_spec_file
from app.vep.utils.vcf_results import _parse_mavedb, _parse_mutfunc

SPEC: ParsingSpec = load_spec_file(SPEC_DIR / "human_grch38.json")


def dump(model):
    """A parser's output as plain data, or None."""
    return model.model_dump() if model is not None else None


def run(plugin: str, csq_values):
    spec = SPEC.plugin(plugin)
    assert spec is not None, f"no spec for {plugin}"
    return apply_plugin_spec(csq_values, INDEX_MAP, spec)


# --- the spec document itself ------------------------------------------------


def test_bundled_spec_validates():
    """The shipped JSON round-trips through the strict model."""
    assert SPEC.spec_version
    assert {p.plugin for p in SPEC.plugins} == {"mutfunc", "mavedb"}


def test_unknown_key_is_rejected():
    """extra=forbid: a spec we don't understand fails at load, not at parse."""
    with pytest.raises(Exception):
        ParsingSpec.model_validate(
            {"spec_version": "x", "plugins": [], "surprise": True}
        )


def test_zip_requires_matching_as_entries():
    with pytest.raises(Exception):
        TargetSpec.model_validate(
            {
                "field": "assays",
                "from": ["a", "b"],
                "transform": "zip",
                "as": [{"field": "only_one"}],
            }
        )


# --- mutfunc: four scalars ---------------------------------------------------

MUTFUNC_SCORES = dict(
    mutfunc_motif="0.1", mutfunc_int="0.2", mutfunc_mod="0.3", mutfunc_exp="0.4"
)


def test_mutfunc_matches_hand_written_parser():
    csq = row_list(**MUTFUNC_SCORES)
    assert run("mutfunc", csq) == dump(_parse_mutfunc(csq, INDEX_MAP))


def test_mutfunc_empty_matches():
    assert run("mutfunc", EMPTY) == dump(_parse_mutfunc(EMPTY, INDEX_MAP)) == None


def test_mutfunc_partial_matches():
    """Only some scores present: the rest must come back None, not be dropped."""
    csq = row_list(mutfunc_motif="0.1", mutfunc_exp="0.4")
    assert run("mutfunc", csq) == dump(_parse_mutfunc(csq, INDEX_MAP))


# --- MaveDB: positional zip, the hard case -----------------------------------

MAVEDB_MULTI = dict(
    MaveDB_score="1.5&2.5&NA",
    MaveDB_urn="urn:1&urn:2&urn:3",
    MaveDB_doi="10.1/a&NA&10.1/c",
    MaveDB_nt="c.1A>G&NA",
    MaveDB_pro="p.Lys1Arg&NA",
)


def test_mavedb_multi_assay_matches_hand_written_parser():
    csq = row_list(**MAVEDB_MULTI)
    assert run("mavedb", csq) == dump(_parse_mavedb(csq, INDEX_MAP))


def test_mavedb_multi_assay_shape_is_as_expected():
    """Pin the actual expected values, so a bug in *both* paths can't pass."""
    result = run("mavedb", row_list(**MAVEDB_MULTI))
    assert result["protein_variant"] == "p.Lys1Arg"
    assert [(a["urn"], a["score"]) for a in result["assays"]] == [
        ("urn:1", 1.5),
        ("urn:2", 2.5),
        ("urn:3", None),  # NA score, but a real urn -> assay kept
    ]


def test_mavedb_empty_matches():
    assert run("mavedb", EMPTY) == dump(_parse_mavedb(EMPTY, INDEX_MAP)) == None


def test_mavedb_uneven_columns_match():
    """Fewer scores than urns: `align: max` must pad, not truncate."""
    csq = row_list(MaveDB_score="1.5", MaveDB_urn="urn:1&urn:2")
    assert run("mavedb", csq) == dump(_parse_mavedb(csq, INDEX_MAP))


def test_mavedb_protein_variant_only_is_none():
    """pro present but no score/urn -> no assays -> whole annotation is None
    (require_any_output), matching the hand-written parser."""
    csq = row_list(MaveDB_pro="p.Lys1Arg")
    assert run("mavedb", csq) == dump(_parse_mavedb(csq, INDEX_MAP)) == None


def test_mavedb_all_na_assay_dropped():
    """A position where both score and urn are NA is dropped entirely."""
    csq = row_list(MaveDB_score="1.5&NA", MaveDB_urn="urn:1&NA")
    assert run("mavedb", csq) == dump(_parse_mavedb(csq, INDEX_MAP))
