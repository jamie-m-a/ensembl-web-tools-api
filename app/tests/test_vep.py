from io import StringIO
import gzip
import os
import tempfile
from pydantic import FilePath
import pytest

from app.vep.models import vcf_results_model as model
from app.vep.utils.vcf_results import (
    get_results_from_path,
    get_results_from_stream,
    _get_prediction_index_map,
)
from app.vep.utils.tsv_export import stream_vep_tsv, gzip_text_stream

# A representative CSQ column list, used to build a CSQ header fixture for the
# index-map test. `_get_prediction_index_map` indexes whatever columns a header
# declares, so any distinct list exercises it.
TARGET_COLUMNS = [
    "Allele", "AF", "Consequence", "Feature", "Feature_type", "BIOTYPE",
    "CANONICAL", "SYMBOL", "Gene", "STRAND", "IMPACT",
    "MANE_SELECT", "MANE_PLUS_CLINICAL",
    "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "DOMAINS",
    "ProtVar_stability", "ProtVar_int", "ProtVar_pocket",
    "IntAct_feature_ac", "IntAct_feature_type", "IntAct_interaction_ac",
    "mutfunc_motif", "mutfunc_int", "mutfunc_mod", "mutfunc_exp",
    "MaveDB_score", "MaveDB_urn", "MaveDB_doi", "MaveDB_nt", "MaveDB_pro",
]
from app.vep.utils.vcf_results import (
    _set_allele_type,
    _get_alt_allele_details,
    _get_csq_value,
)

CSQ_DESCRIPTION = "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|REF_ALLELE|UPLOADED_ALLELE|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|SIFT|PolyPhen|AF|CLIN_SIG|SOMATIC|PHENO|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS"

CSQ_1 = "T|upstream_gene_variant|MODIFIER|FAM138F|ENSG00000282591|Transcript|ENST00000631376.1|lncRNA||||||||||rs868831437|C|C/T|4978|-1||HGNC|HGNC:33581|YES|||0.4860||||||||"

CSQ_2 = (
    "A|intergenic_variant|MODIFIER|||||||||||||||rs1555675005|T|T/A|||||||||||||||||"
)

CSQ_NO_FREQ = "T|upstream_gene_variant|MODIFIER|FAM138F|ENSG00000282591|Transcript|ENST00000631376.1|lncRNA||||||||||rs868831437|C|C/T|4978|-1||HGNC|HGNC:33581|YES|||||||||||"

CSQ_NO_CON = "T||MODIFIER|FAM138F|ENSG00000282591|Transcript|ENST00000631376.1|lncRNA||||||||||rs868831437|C|C/T|4978|-1||HGNC|HGNC:33581|YES|||||||||||"

TEST_VCF = f"""##fileformat=VCFv4.2
##fileDate=20160824
##INFO=<ID=CSQ,Number=.,Type=String,Description="{CSQ_DESCRIPTION}">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr19	82664	.	C	T	50	PASS	CSQ={CSQ_1}
chr19	82829	my_var	T	A	50	PASS	CSQ={CSQ_2}
chrX	982829	.	G	C	.	50	PASS    .

"""

TEST_EMPTY_VCF = f"""##fileformat=VCFv4.2
##fileDate=20160824
##INFO=<ID=CSQ,Number=.,Type=String,Description="{CSQ_DESCRIPTION}">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

"""

TEST_PAGING_VCF = f"""##fileformat=VCFv4.2
##fileDate=20160824
##INFO=<ID=CSQ,Number=.,Type=String,Description="{CSQ_DESCRIPTION}">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr19	82664	id_01	C	T	50	PASS	CSQ={CSQ_1}
chr19	82829	id_02	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_03	T	A	50	PASS	.
chr19	82829	id_04	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_05	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_06	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_07	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_08	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_09	T	A	50	PASS	.
chr19	82829	id_10	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_11	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_12	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_13	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_14	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_15	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_16	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_17	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_18	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_19	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_20	T	A	50	PASS	CSQ={CSQ_2}
chr19	82829	id_21	T	A	50	PASS	CSQ={CSQ_2}
"""

VCF_PATH = FilePath("tests/test_vep.vcf")

def test_get_prediction_index_map():

    expected_index = {TARGET_COLUMNS[x]: x for x in range(0, len(TARGET_COLUMNS))}

    csq_header = f"""Consequence annotations from Ensembl VEP. Format: {'|'.join(TARGET_COLUMNS)}"""

    prediction_index_map = _get_prediction_index_map(csq_header)
    assert prediction_index_map == expected_index


def test_set_allele_type():

    outcomes = {
        "SNV": (True, True, True),
        "deletion": (True, False, False),
        "insertion": (False, True, False),
        "indel": (False, False, False),
        "substitution": (False, False, True),
    }

    for expected, args in outcomes.items():
        assert _set_allele_type(*args)[0] == expected


def test_get_csq_value():
    index_map = {
        "TEST_STR": 0,
        "TEST_NUM": 1,
        "TEST_BOOL": 2,
        "TEST_EMPTY": 3,
    }
    csq_values = ["foo", 2, True, ""]

    assert _get_csq_value(csq_values, "TEST_STR", "ERROR", index_map) == "foo"
    assert _get_csq_value(csq_values, "TEST_NUM", -1, index_map) == 2
    assert _get_csq_value(csq_values, "TEST_BOOL", False, index_map)
    assert _get_csq_value(csq_values, "TEST_MISSING", "ERROR", index_map) == "ERROR"
    assert _get_csq_value(csq_values, "TEST_EMPTY", None, index_map) == None


def test_get_alt_allele_details():
    index_map = _get_prediction_index_map(CSQ_DESCRIPTION)
    csq_list = [CSQ_1, CSQ_2, CSQ_NO_FREQ]

    results = _get_alt_allele_details("C", "T", csq_list, index_map)
    #assert type(results) == model.AlternativeVariantAllele
    assert results.allele_sequence == "T"
    assert results.allele_type == "SNV"
    assert results.representative_population_allele_frequency == 0.4860
    assert len(results.predicted_molecular_consequences) == 2
    #assert (
    #    results.predicted_molecular_consequences[0].feature_type
    #    == model.FeatureType.transcript
    #)
    assert results.predicted_molecular_consequences[0].biotype == "lncRNA"
    assert results.predicted_molecular_consequences[0].gene_symbol == "FAM138F"


def test_get_alt_allele_no_consequence():
    index_map = _get_prediction_index_map(CSQ_DESCRIPTION)

    csq_list = [CSQ_NO_CON]

    results = _get_alt_allele_details("C", "T", csq_list, index_map)

    #assert type(results) == model.AlternativeVariantAllele
    assert results.allele_sequence == "T"
    assert len(results.predicted_molecular_consequences) == 1
    assert results.predicted_molecular_consequences[0].consequences == []


def test_get_alt_allele_details_intergenic():
    index_map = _get_prediction_index_map(CSQ_DESCRIPTION)

    csq_list = [CSQ_2]

    # model.AlternativeVariantAllele
    results = _get_alt_allele_details("C", "A", csq_list, index_map)

    #assert type(results) == model.AlternativeVariantAllele
    assert results.allele_sequence == "A"
    assert results.allele_type == "SNV"
    assert len(results.predicted_molecular_consequences) == 1
    assert results.predicted_molecular_consequences[0].feature_type == None
    assert len(results.predicted_molecular_consequences[0].consequences) == 1
    assert (
        results.predicted_molecular_consequences[0].consequences[0]
        == "intergenic_variant"
    )

@pytest.mark.skip(reason="Unknown bug")
def test_get_results_from_stream():
    variant_count = 3
    results = get_results_from_stream(100, 1, variant_count, StringIO(TEST_VCF))

    print(results.variants)
    assert len(results.variants) == variant_count

    assert results.metadata.pagination.page == 1
    assert results.metadata.pagination.per_page == 100
    assert results.metadata.pagination.total == variant_count

    assert results.variants[0].name == "."
    assert results.variants[1].name == "my_var"
    assert results.variants[2].name == "."

    assert results.variants[0].reference_allele.allele_sequence == "C"
    assert results.variants[1].reference_allele.allele_sequence == "T"
    assert results.variants[2].reference_allele.allele_sequence == "G"

    assert results.variants[0].alternative_alleles[0].allele_sequence == "T"
    assert results.variants[0].alternative_alleles[0].allele_type == "SNV"

    assert (
        results.variants[0]
        .alternative_alleles[0]
        .representative_population_allele_frequency
        == 0.4860
    )
    assert (
        results.variants[1]
        .alternative_alleles[0]
        .representative_population_allele_frequency
        == None
    )

def test_paging():
    variant_count = 21
    results = get_results_from_path(5, 1, VCF_PATH)

    assert(results.metadata.pagination.page == 1)
    assert(results.metadata.pagination.per_page == 5)
    assert results.metadata.pagination.total == variant_count

    assert(results.variants[0].name == "id_01")
    assert(results.variants[-1].name == "id_05")

    results = get_results_from_path(5, 2, VCF_PATH)
    assert(results.variants[0].name == "id_06")
    assert(results.variants[-1].name == "id_10")

    results = get_results_from_path(5, 3, VCF_PATH)
    assert(results.variants[0].name == "id_11")
    assert(results.variants[-1].name == "id_15")

    results = get_results_from_path(5, 4, VCF_PATH)
    assert(results.variants[0].name == "id_16")
    assert(results.variants[-1].name == "id_20")

    #results = get_results_from_path(5, 5, VCF_PATH)
    #assert(results.variants[0].name == "id_21")
    #assert(len(results.variants) == 1)

def test_negative_paging():
    results = get_results_from_path(5, 6, VCF_PATH)
    assert(len(results.variants) == 0)
    assert(results.metadata.pagination.total == 21)


@pytest.mark.skip(reason="Used to test against a real VCF file")
def test_get_results_with_file_and_dump():

    vcf_path = (
        #"/Users/jon/Programming/vep-vcf-results/vep-output-phase1-options-plus-con.vcf"
        "/Users/jon/Programming/ensembl-web-tools-api/test_VEP.vcf.gz"
    )
    results = get_results_from_path(100, 1, vcf_path)

    with open("dump.json", "w") as test_dump:
        test_dump.write(results.json())

    #assert results.variants[0].name =="rs1405511870"
    assert results.metadata.pagination.total == 1
    assert len(results.variants) == 1


def test_gzip_text_stream_roundtrip():
    """The download compressor yields a valid gzip container that decompresses
    back to the concatenated input text."""
    chunks = ["Uploaded_variation\tLocation\n", "rs1\t19:100\n", "rs2\t19:200\n"]
    compressed = b"".join(gzip_text_stream(iter(chunks)))
    assert compressed[:2] == b"\x1f\x8b"  # gzip magic
    assert gzip.decompress(compressed).decode() == "".join(chunks)


def test_stream_vep_tsv_gzip_matches_plain(tmp_path):
    """Gzipping the TSV stream decompresses to the same bytes as the plain TSV."""
    vcf_path = tmp_path / "output.vcf.gz"
    with gzip.open(vcf_path, "wt") as f:
        f.write(TEST_VCF)

    plain = "".join(stream_vep_tsv(FilePath(vcf_path)))
    compressed = b"".join(gzip_text_stream(stream_vep_tsv(FilePath(vcf_path))))
    assert gzip.decompress(compressed).decode() == plain


def test_stream_vep_tsv(tmp_path):
    """The flattened TSV export: one row per CSQ entry (transcript AND
    intergenic), with location/ref columns plus every CSQ field."""
    vcf_path = tmp_path / "output.vcf.gz"
    with gzip.open(vcf_path, "wt") as f:
        f.write(TEST_VCF)

    rows = list(stream_vep_tsv(FilePath(vcf_path)))
    lines = [r.rstrip("\n").split("\t") for r in rows]

    header = lines[0]
    # Uploaded_variation, Location, Ref + every CSQ field.
    assert header[:4] == ["Uploaded_variation", "Location", "Ref", "Allele"]
    assert header == ["Uploaded_variation", "Location", "Ref"] + CSQ_DESCRIPTION.split("Format: ")[1].split("|")

    data = lines[1:]
    # Two records carry CSQ (CSQ_1, CSQ_2); the third has no CSQ and is skipped.
    assert len(data) == 2

    first = dict(zip(header, data[0]))
    assert first["Uploaded_variation"] == "."
    assert first["Location"] == "19:82664"  # "chr" prefix stripped
    assert first["Ref"] == "C"
    assert first["Allele"] == "T"
    assert first["Consequence"] == "upstream_gene_variant"

    # The intergenic entry is fully expanded too.
    intergenic = dict(zip(header, data[1]))
    assert intergenic["Consequence"] == "intergenic_variant"
    assert intergenic["Location"] == "19:82829"
