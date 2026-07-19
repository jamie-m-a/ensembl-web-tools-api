"""The generic spec-driven `annotations` emitted on alleles and transcript
consequences (the additive go-flat wire format).

Checks the *wiring*: that _get_alt_allele_details drives the pinned parsing
spec's plugins through apply_plugin_spec and attaches the results at the right
scope, additively (the typed fields stay). apply_plugin_spec's own correctness
(== the hand-written _parse_* bank) is proven separately by test_spec_real_vcfs.
"""

from app.vep.utils.spec_interpreter import apply_plugin_spec
from app.vep.utils.spec_loader import load_merged_spec
from app.vep.utils.vcf_results import _get_alt_allele_details

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


def test_emit_is_additive_typed_fields_remain():
    allele = _get_alt_allele_details("A", "T", [ROW], INDEX_MAP, SPEC)
    # The typed fields are still populated beside the generic annotations.
    assert allele.clinvar is not None
    assert allele.frequencies is not None
    assert allele.predicted_molecular_consequences[0].annotations  # non-empty


def test_no_spec_means_no_generic_annotations():
    allele = _get_alt_allele_details("A", "T", [ROW], INDEX_MAP, None)
    assert allele.annotations == []
    assert allele.predicted_molecular_consequences[0].annotations == []
    # typed fields unaffected by the absence of a spec
    assert allele.clinvar is not None
