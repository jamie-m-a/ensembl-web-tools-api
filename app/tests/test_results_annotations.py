"""The generic spec-driven `annotations` emitted on alleles and transcript
consequences (the additive go-flat wire format).

Checks the *wiring*: that _get_alt_allele_details drives the pinned parsing
spec's plugins through apply_plugin_spec and attaches the results at the right
scope. Since the go-flat cutover these are the only annotation data on the
response; apply_plugin_spec's own correctness is covered by
test_spec_interpreter.
"""

from app.vep.utils.spec_interpreter import apply_plugin_spec
from app.vep.utils.spec_loader import load_merged_spec
from app.vep.utils.vcf_results import _gate_af_columns, _get_alt_allele_details

SPEC = load_merged_spec("human_grch38").parsing

# A CSQ header layout with structural columns plus one allele-scope frequency
# (gnomad_exomes, incl. a per-population column to exercise the pattern_map
# flat-frequency shape), one allele-scope custom (clinvar), and one
# transcript-scope simple plugin (revel).
COLUMNS = [
    "Allele", "Feature_type", "Consequence", "Feature", "Gene", "BIOTYPE",
    "CANONICAL", "STRAND", "SYMBOL",
    "REVEL", "gnomAD_exomes_AF", "gnomAD_exomes_AF_nfe", "ClinVar_CLNSIG",
]
INDEX_MAP = {column: i for i, column in enumerate(COLUMNS)}
ROW = "|".join([
    "T", "Transcript", "missense_variant", "ENST001", "ENSG001", "protein_coding",
    "YES", "1", "BRCA2",
    "0.9", "0.01", "0.02", "Pathogenic",
])
CSQ_VALUES = ROW.split("|")


def _expected(plugin_name):
    return apply_plugin_spec(CSQ_VALUES, INDEX_MAP, SPEC.plugin(plugin_name))


def test_allele_scope_annotations_attached():
    allele = _get_alt_allele_details("A", "T", [ROW], INDEX_MAP, SPEC)
    by_plugin = {a.plugin: a for a in allele.annotations}

    assert by_plugin.keys() >= {"gnomad_exomes", "clinvar"}
    assert all(a.scope == "allele" for a in allele.annotations)
    assert by_plugin["gnomad_exomes"].data == _expected("gnomad_exomes")
    assert by_plugin["clinvar"].data == _expected("clinvar")
    # the flat frequency shape carries the per-population column
    assert by_plugin["gnomad_exomes"].data["populations"] == {"nfe": 0.02}


def test_transcript_scope_annotations_attached():
    allele = _get_alt_allele_details("A", "T", [ROW], INDEX_MAP, SPEC)
    consequence = allele.predicted_molecular_consequences[0]
    by_plugin = {a.plugin: a for a in consequence.annotations}

    assert "revel" in by_plugin
    assert all(a.scope == "transcript" for a in consequence.annotations)
    assert by_plugin["revel"].data == _expected("revel")


def test_annotations_are_the_only_annotation_data():
    allele = _get_alt_allele_details("A", "T", [ROW], INDEX_MAP, SPEC)
    # What used to be the typed `clinvar` / `frequencies` fields now arrives
    # only as generic annotations, at allele scope.
    by_plugin = {a.plugin: a.data for a in allele.annotations}
    assert by_plugin["clinvar"]["significance"] == ["Pathogenic"]
    assert by_plugin["gnomad_exomes"]["overall"] == 0.01
    assert allele.predicted_molecular_consequences[0].annotations  # non-empty


def test_no_spec_means_no_generic_annotations():
    allele = _get_alt_allele_details("A", "T", [ROW], INDEX_MAP, None)
    assert allele.annotations == []
    assert allele.predicted_molecular_consequences[0].annotations == []
    # the envelope is unaffected by the absence of a spec
    assert allele.allele_sequence == "T"
    assert allele.predicted_molecular_consequences[0].gene_symbol == "BRCA2"


# --- AF-population emission gate ---------------------------------------------
# A full-cache VCF carries every ancestry; a job that selected only some AF
# populations must still show only those. The parser's pattern_map reads every
# column present, so the gate (in _with_display_panels) trims the served
# annotation to the pinned expected columns.

_GATE_COLUMNS = [
    "Allele", "Feature_type", "Consequence", "Feature", "Gene", "BIOTYPE",
    "CANONICAL", "STRAND", "SYMBOL",
    "gnomAD_exomes_AF",
    "gnomAD_exomes_AF_nfe", "gnomAD_exomes_AF_eas", "gnomAD_exomes_AF_afr",
]
_GATE_INDEX = {column: i for i, column in enumerate(_GATE_COLUMNS)}
_GATE_ROW = "|".join([
    "T", "Transcript", "missense_variant", "ENST001", "ENSG001", "protein_coding",
    "YES", "1", "BRCA2",
    "0.01", "0.02", "0.03", "0.04",
])


def test_af_columns_gated_to_selected_populations_only():
    allele = _get_alt_allele_details("A", "T", [_GATE_ROW], _GATE_INDEX, SPEC)
    gnomad = {a.plugin: a for a in allele.annotations}["gnomad_exomes"]
    # every ancestry in the (full-cache) VCF is parsed before gating
    assert set(gnomad.data["populations"]) == {"nfe", "eas", "afr"}
    assert gnomad.data["overall"] == 0.01

    # the submission selected only the nfe population column — not the overall
    _gate_af_columns([allele], SPEC, {"gnomAD_exomes_AF_nfe"})

    assert gnomad.data["populations"] == {"nfe": 0.02}
    # the overall's column wasn't selected, so it is gated too (no "All" row)
    assert gnomad.data["overall"] is None


def test_af_columns_keeps_the_overall_when_its_column_is_selected():
    allele = _get_alt_allele_details("A", "T", [_GATE_ROW], _GATE_INDEX, SPEC)
    gnomad = {a.plugin: a for a in allele.annotations}["gnomad_exomes"]
    _gate_af_columns(
        [allele], SPEC, {"gnomAD_exomes_AF", "gnomAD_exomes_AF_eas"}
    )
    assert gnomad.data["overall"] == 0.01
    assert gnomad.data["populations"] == {"eas": 0.03}


def test_af_columns_gate_is_a_no_op_without_a_spec():
    allele = _get_alt_allele_details("A", "T", [_GATE_ROW], _GATE_INDEX, SPEC)
    gnomad = {a.plugin: a for a in allele.annotations}["gnomad_exomes"]
    _gate_af_columns([allele], None, {"gnomAD_exomes_AF_nfe"})
    assert set(gnomad.data["populations"]) == {"nfe", "eas", "afr"}
    assert gnomad.data["overall"] == 0.01
