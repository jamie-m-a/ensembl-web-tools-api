"""Structural-variant + breakend handling in the results parser.

VEP writes a type *word* (e.g. `deletion`) into the CSQ `Allele` column, which the
length heuristic mis-reads as an insertion and the UI renders as a sequence. These
cover classifying SVs from the VCF record instead (see the VCF 4.2 spec, symbolic
+ breakend alternate alleles §1.4.5, §5.3-5.4).
"""

from types import SimpleNamespace

from app.vep.utils.vcf_results import (
    _structural_info,
    _get_variant_type,
    _get_alt_allele_details,
)


def _record(alt_serialized: str, info: dict, ref: str = "N", chrom: str = "1", pos: int = 1000):
    """A minimal vcfpy-record stand-in: one alt (serializing to `alt_serialized`)
    plus an INFO dict, CHROM, POS and REF."""
    alt = SimpleNamespace(serialize=lambda: alt_serialized, value=None)
    return SimpleNamespace(ALT=[alt], INFO=info, REF=ref, CHROM=chrom, POS=pos)


# --- type + allele + detail, per SVTYPE -------------------------------------


def test_deletion_uses_svtype_not_length():
    # The regression: alt="deletion" is 8 chars so the length heuristic called it
    # an insertion. SVTYPE=DEL fixes it.
    sv = _structural_info(_record("<DEL>", {"SVTYPE": "DEL", "SVLEN": 765}))
    assert sv == {"type_word": "deletion", "allele": "<DEL>", "detail": "765 bp"}
    assert _get_variant_type("N", "deletion") == "insertion"  # the old, wrong path


def test_svtype_word_map():
    cases = {
        "INS": "insertion",
        "DUP": "duplication",
        "INV": "inversion",
        "CNV": "copy_number_variation",
    }
    for svtype, word in cases.items():
        sv = _structural_info(_record(f"<{svtype}>", {"SVTYPE": svtype, "END": 2000}))
        assert sv["type_word"] == word
        assert sv["allele"] == f"<{svtype}>"


def test_mobile_element_subtype_collapses_to_base_type():
    sv = _structural_info(
        _record("<DEL:ME:ALU>", {"SVTYPE": "DEL", "SVLEN": -168, "END": 1168})
    )
    assert sv["type_word"] == "deletion"
    assert sv["allele"] == "<DEL:ME:ALU>"  # full symbolic kept for display
    assert sv["detail"] == "168 bp"  # abs(SVLEN)


def test_symbolic_without_svtype_falls_back_to_id_token():
    sv = _structural_info(_record("<DUP>", {"SVLEN": 500}))
    assert sv["type_word"] == "duplication"


# --- span / detail rules -----------------------------------------------------


def test_detail_prefers_abs_svlen():
    sv = _structural_info(_record("<DEL>", {"SVTYPE": "DEL", "SVLEN": -434, "END": 9}))
    assert sv["detail"] == "434 bp"


def test_detail_falls_back_to_end_minus_pos_when_svlen_zero_or_absent():
    # INV: SVLEN=0 -> use END - POS
    inv = _structural_info(_record("<INV>", {"SVTYPE": "INV", "SVLEN": 0, "END": 5000}, pos=1000))
    assert inv["detail"] == "4000 bp"
    # CNV: no SVLEN -> END - POS
    cnv = _structural_info(_record("<CNV>", {"SVTYPE": "CNV", "END": 5000}, pos=1000))
    assert cnv["detail"] == "4000 bp"


def test_insertion_detail_is_svlen_when_no_end():
    ins = _structural_info(_record("<INS>", {"SVTYPE": "INS", "SVLEN": 7093}))
    assert ins["detail"] == "7093 bp"


def test_detail_none_when_neither_svlen_nor_end():
    assert _structural_info(_record("<INS>", {"SVTYPE": "INS"}))["detail"] is None


# --- breakends: shown as <BND> with the two loci, no bases -------------------


def test_breakend_shown_as_bnd_with_both_loci():
    sv = _structural_info(
        _record(
            "G[17:198982[", {"SVTYPE": "BND", "MATEID": ["bnd_X"]},
            ref="G", chrom="2", pos=321681,
        )
    )
    assert sv == {
        "type_word": "breakend",
        "allele": "<BND>",
        "detail": "2:321681 ↔ 17:198982",
    }


def test_breakend_mate_parsed_from_other_bracket_orientation():
    # `[p[t` orientation (base trailing) — mate still parsed from between brackets.
    sv = _structural_info(_record("[2:321681[C", {}, ref="C", chrom="17", pos=198982))
    assert sv["type_word"] == "breakend"
    assert sv["allele"] == "<BND>"
    assert sv["detail"] == "17:198982 ↔ 2:321681"


# --- simple variants are untouched ------------------------------------------


def test_simple_variants_are_not_structural():
    assert _structural_info(_record("T", {}, ref="A")) is None
    assert _structural_info(_record("ATTTT", {}, ref="A")) is None
    assert _get_variant_type("A", "T") == "SNV"
    assert _get_variant_type("A", "ATT") == "insertion"
    assert _get_variant_type("ATT", "A") == "deletion"


# --- _get_alt_allele_details applies the sv override ------------------------


def test_alt_allele_details_uses_sv_for_type_sequence_and_detail():
    sv = {"type_word": "deletion", "allele": "<DEL>", "detail": "765 bp"}
    allele = _get_alt_allele_details("N", "deletion", [], {}, None, sv)
    assert allele.allele_type == "deletion"
    assert allele.allele_sequence == "<DEL>"
    assert allele.structural_variant_detail == "765 bp"


def test_alt_allele_details_simple_variant_has_no_detail():
    allele = _get_alt_allele_details("A", "T", [], {}, None, None)
    assert allele.allele_type == "SNV"
    assert allele.allele_sequence == "T"
    assert allele.structural_variant_detail is None
