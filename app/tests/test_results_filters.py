"""Tests for server-side results filtering (app/vep/utils/results_filters.py and
the filtered scan path in vcf_results.get_results_from_path).

Filtered requests can't use the BGZF page index, so they scan the whole file
with gzip.open — meaning a plain gzip VCF fixture is enough here (no BGZF/page
index needed).
"""

import gzip

import pytest
from pydantic import FilePath

from app.vep.utils import results_filters as rf
from app.vep.utils.vcf_results import get_results_from_path

CSQ_DESC = (
    "Consequence annotations from Ensembl VEP. Format: "
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE"
)

# CSQ column -> index for the format above (Consequence is column 1).
INDEX_MAP = {
    name: i
    for i, name in enumerate(
        "Allele Consequence IMPACT SYMBOL Gene Feature_type Feature BIOTYPE".split()
    )
}

# A wider layout that also carries the canonical / MANE columns, for group tests.
GROUP_COLUMNS = (
    "Allele Consequence Feature CANONICAL MANE_SELECT MANE_PLUS_CLINICAL".split()
)
GROUP_INDEX_MAP = {name: i for i, name in enumerate(GROUP_COLUMNS)}


def _group_entry(feature: str, canonical: str, mane_select: str, mane_plus: str) -> str:
    return "|".join(["T", "missense_variant", feature, canonical, mane_select, mane_plus])


def _group_record(pos: int, entries: list[str]) -> str:
    return f"chr1\t{100 + pos}\tid_{pos:02d}\tC\tT\t.\t.\tCSQ={','.join(entries)}\n"


def _record(pos: int, consequences: list[str]) -> str:
    """One VCF data line whose CSQ carries one entry per given consequence."""
    entries = ",".join(
        f"T|{cons}|MODERATE|GENE{pos}|ENSG{pos}|Transcript|ENST{pos}|protein_coding"
        for cons in consequences
    )
    return f"chr1\t{100 + pos}\tid_{pos:02d}\tC\tT\t.\t.\tCSQ={entries}\n"


def _write_vcf(path, records: list[str]) -> str:
    text = (
        "##fileformat=VCFv4.2\n"
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="{CSQ_DESC}">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(records)
    )
    vcf_path = path / "results.vcf.gz"
    with gzip.open(vcf_path, "wt") as handle:
        handle.write(text)
    return str(vcf_path)


def _consequence_filter(*values: str) -> rf.ResultsFilter:
    return rf.ResultsFilter(
        field=rf.CONSEQUENCE_FIELD, operator=rf.OPERATOR_IN, values=list(values)
    )


def _transcript_record(pos: int, features: list[tuple[str, str]]) -> str:
    """A VCF data line with one CSQ entry per (transcript_feature, consequence)."""
    entries = ",".join(
        f"T|{cons}|MODERATE|GENE|ENSG|Transcript|{feature}|protein_coding"
        for feature, cons in features
    )
    return f"chr1\t{100 + pos}\tid_{pos:02d}\tC\tT\t.\t.\tCSQ={entries}\n"


def _transcript_filter(*values: str) -> rf.ResultsFilter:
    return rf.ResultsFilter(
        field=rf.TRANSCRIPT_FIELD, operator=rf.OPERATOR_IN, values=list(values)
    )


# --- parse_filters -----------------------------------------------------------


def test_parse_filters_empty_returns_none_list():
    assert rf.parse_filters(None) == []
    assert rf.parse_filters("") == []


def test_parse_filters_valid():
    parsed = rf.parse_filters(
        '[{"field": "consequence", "operator": "in", "values": ["missense_variant"]}]'
    )
    assert len(parsed) == 1
    assert parsed[0].field == "consequence"
    assert parsed[0].values == ["missense_variant"]


def test_parse_filters_bad_json_raises():
    with pytest.raises(rf.FilterError):
        rf.parse_filters("not json")


def test_parse_filters_non_array_raises():
    with pytest.raises(rf.FilterError):
        rf.parse_filters('{"field": "consequence"}')


# --- extract_csq_entries -----------------------------------------------------


def test_extract_csq_entries_splits_entries_and_subfields():
    line = _record(1, ["missense_variant", "synonymous_variant"])
    entries = rf.extract_csq_entries(line)
    assert len(entries) == 2
    assert entries[0][INDEX_MAP["Consequence"]] == "missense_variant"
    assert entries[1][INDEX_MAP["Consequence"]] == "synonymous_variant"


def test_extract_csq_entries_no_csq():
    assert rf.extract_csq_entries("chr1\t1\t.\tC\tT\t.\t.\tAC=1\n") == []


# --- compile_filters + pipeline ----------------------------------------------


def test_compile_rejects_unknown_field():
    with pytest.raises(rf.FilterError):
        rf.compile_filters(
            [rf.ResultsFilter(field="nonsense", operator="in", values=["x"])],
            INDEX_MAP,
        )


def test_compile_rejects_unknown_operator():
    with pytest.raises(rf.FilterError):
        rf.compile_filters(
            [rf.ResultsFilter(field="consequence", operator="gt", values=["x"])],
            INDEX_MAP,
        )


def test_pipeline_keeps_matching_and_counts_removed():
    lines = [
        _record(1, ["missense_variant"]),
        _record(2, ["synonymous_variant"]),
        _record(3, ["intron_variant"]),
        _record(4, ["missense_variant"]),
    ]
    compiled = rf.compile_filters([_consequence_filter("missense_variant")], INDEX_MAP)
    kept, stats = rf.apply_filter_pipeline(lines, compiled)

    assert len(kept) == 2
    assert kept == [lines[0], lines[3]]
    assert stats[0].field == "consequence"
    assert stats[0].removed == 2


def test_pipeline_prunes_nonmatching_entries():
    lines = [
        # kept, pruned to just the synonymous entry
        _record(1, ["missense_variant", "synonymous_variant"]),
        # kept intact: its single &-joined entry intersects the selection
        _record(2, ["splice_region_variant&intron_variant"]),
        # removed: neither term selected
        _record(3, ["upstream_gene_variant"]),
    ]
    compiled = rf.compile_filters(
        [_consequence_filter("synonymous_variant", "intron_variant")], INDEX_MAP
    )
    kept, stats = rf.apply_filter_pipeline(lines, compiled)

    assert len(kept) == 2
    ci = INDEX_MAP["Consequence"]
    # record 1 lost its non-matching missense entry
    assert [e[ci] for e in rf.extract_csq_entries(kept[0])] == ["synonymous_variant"]
    # record 2 kept as-is
    assert [e[ci] for e in rf.extract_csq_entries(kept[1])] == [
        "splice_region_variant&intron_variant"
    ]
    assert stats[0].removed == 1


# --- filter_records (streaming, page-bounded) --------------------------------


def test_filter_records_retains_only_the_page_window():
    lines = [_record(i, ["missense_variant"]) for i in range(1, 8)]  # 7 matches
    compiled = rf.compile_filters([_consequence_filter("missense_variant")], INDEX_MAP)

    outcome = rf.filter_records(lines, compiled, start=2, count=3)

    # Full counts are tallied regardless of which window is asked for...
    assert outcome.scanned_total == 7
    assert outcome.matched_total == 7
    # ...but only the [2, 5) slice of survivors is materialised.
    assert outcome.page == [lines[2], lines[3], lines[4]]
    assert outcome.stats[0].removed == 0


def test_filter_records_counts_all_scanned_including_dropped():
    lines = [
        _record(1, ["missense_variant"]),
        _record(2, ["synonymous_variant"]),  # dropped
        _record(3, ["missense_variant"]),
    ]
    compiled = rf.compile_filters([_consequence_filter("missense_variant")], INDEX_MAP)

    outcome = rf.filter_records(lines, compiled, start=0, count=10)

    assert outcome.scanned_total == 3
    assert outcome.matched_total == 2
    assert outcome.page == [lines[0], lines[2]]
    assert outcome.stats[0].removed == 1


def test_filter_records_page_stays_bounded_below_match_count():
    # 100 matches but a page of 5: only 5 lines are ever held, and a lazy
    # iterator source is consumed without materialising the input.
    lines = [_record(i, ["missense_variant"]) for i in range(1, 101)]
    compiled = rf.compile_filters([_consequence_filter("missense_variant")], INDEX_MAP)

    outcome = rf.filter_records(iter(lines), compiled, start=0, count=5)

    assert outcome.matched_total == 100
    assert len(outcome.page) == 5


def test_apply_filter_pipeline_wrapper_keeps_every_survivor():
    lines = [
        _record(1, ["missense_variant"]),
        _record(2, ["intron_variant"]),
        _record(3, ["missense_variant"]),
    ]
    compiled = rf.compile_filters([_consequence_filter("missense_variant")], INDEX_MAP)

    kept, stats = rf.apply_filter_pipeline(lines, compiled)

    assert kept == [lines[0], lines[2]]
    assert stats[0].removed == 1


# --- raw-line pre-filter (fast reject before splitting the CSQ) ---------------


def test_membership_filters_carry_a_line_prefilter_gene_symbol_does_not():
    cases = [
        (_consequence_filter("missense_variant"), True),
        (_transcript_filter("ENST00000341065"), True),
        (rf.ResultsFilter(field=rf.GENE_ID_FIELD, operator=rf.OPERATOR_IN, values=["ENSG001"]), True),
        # Gene symbol matching is case-insensitive, so a case-sensitive substring
        # test could false-negative — no prefilter.
        (rf.ResultsFilter(field=rf.GENE_SYMBOL_FIELD, operator=rf.OPERATOR_IN, values=["BRCA1"]), False),
    ]
    for filt, has_prefilter in cases:
        (cf,) = rf.compile_filters([filt], INDEX_MAP)
        assert (cf.line_prefilter is not None) is has_prefilter


def test_prefilter_substring_false_positive_is_still_excluded_by_exact_check():
    # The selected term appears in the line only as the SYMBOL, not as a
    # Consequence: the cheap substring prefilter admits the line, but the exact
    # per-entry check must still drop it (the prefilter is necessary, not sufficient).
    line = (
        "chr1\t100\tv\tC\tT\t.\t.\t"
        "CSQ=T|synonymous_variant|MODERATE|missense_variant|ENSG|Transcript|ENST|protein_coding\n"
    )
    compiled = rf.compile_filters([_consequence_filter("missense_variant")], INDEX_MAP)

    assert compiled[0].line_prefilter(line) is True  # prefilter alone would admit it
    outcome = rf.filter_records([line], compiled)
    assert outcome.matched_total == 0  # ...but the exact check drops it
    assert outcome.stats[0].removed == 1


def test_prefilter_rejects_without_error_when_no_token_present():
    lines = [_record(1, ["synonymous_variant"]), _record(2, ["intron_variant"])]
    compiled = rf.compile_filters([_consequence_filter("missense_variant")], INDEX_MAP)
    outcome = rf.filter_records(lines, compiled)
    assert outcome.matched_total == 0
    assert outcome.scanned_total == 2
    assert outcome.stats[0].removed == 2


def test_transcript_group_filter_has_no_prefilter():
    # Transcript-group tests read CANONICAL/MANE columns, not literal values.
    (cf,) = rf.compile_filters([_transcript_group_filter("canonical")], GROUP_INDEX_MAP)
    assert cf.line_prefilter is None


def test_transcript_filter_matches_ignoring_version():
    lines = [
        _transcript_record(
            1, [("ENST00000341065.8", "missense_variant"), ("ENST00000999.2", "intron_variant")]
        ),
        _transcript_record(2, [("ENST00000111.1", "missense_variant")]),
    ]
    fi = INDEX_MAP["Feature"]
    # user supplies the id without a version; the file has versioned ids
    compiled = rf.compile_filters([_transcript_filter("ENST00000341065")], INDEX_MAP)
    kept, stats = rf.apply_filter_pipeline(lines, compiled)

    assert len(kept) == 1
    kept_features = [e[fi] for e in rf.extract_csq_entries(kept[0])]
    # only the matching transcript survives; the other transcript is pruned
    assert kept_features == ["ENST00000341065.8"]
    assert stats[0].removed == 1


def test_transcript_filter_matches_with_version_supplied():
    lines = [_transcript_record(1, [("ENST00000341065.8", "missense_variant")])]
    # a versioned id supplied by the user still matches (version ignored)
    compiled = rf.compile_filters([_transcript_filter("ENST00000341065.3")], INDEX_MAP)
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert len(kept) == 1


def test_consequence_and_transcript_combined():
    # transcript A is missense, transcript B (same variant) is intron.
    lines = [
        _transcript_record(
            1,
            [
                ("ENST_A.1", "missense_variant"),
                ("ENST_B.1", "intron_variant"),
            ],
        )
    ]
    fi = INDEX_MAP["Feature"]
    # consequence=missense AND transcript in {A, B}: only A satisfies both
    compiled = rf.compile_filters(
        [_consequence_filter("missense_variant"), _transcript_filter("ENST_A", "ENST_B")],
        INDEX_MAP,
    )
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert len(kept) == 1
    assert [e[fi] for e in rf.extract_csq_entries(kept[0])] == ["ENST_A.1"]


def _gene_record(pos: int, genes: list[tuple[str, str]]) -> str:
    """A VCF data line with one CSQ entry per (SYMBOL, Gene) pair."""
    entries = ",".join(
        f"T|missense_variant|MODERATE|{symbol}|{gene_id}|Transcript|ENST_{i}|protein_coding"
        for i, (symbol, gene_id) in enumerate(genes)
    )
    return f"chr1\t{100 + pos}\tid_{pos:02d}\tC\tT\t.\t.\tCSQ={entries}\n"


def test_gene_symbol_filter_case_insensitive_and_prunes():
    lines = [
        _gene_record(1, [("TP53", "ENSG00000141510"), ("BRCA1", "ENSG00000012048")]),
        _gene_record(2, [("EGFR", "ENSG00000146648")]),
    ]
    si = INDEX_MAP["SYMBOL"]
    # lower-case query still matches TP53
    compiled = rf.compile_filters(
        [rf.ResultsFilter(field=rf.GENE_SYMBOL_FIELD, operator=rf.OPERATOR_IN, values=["tp53"])],
        INDEX_MAP,
    )
    kept, stats = rf.apply_filter_pipeline(lines, compiled)

    assert len(kept) == 1
    # only the TP53 entry survives; the BRCA1 entry on the same record is pruned
    assert [e[si] for e in rf.extract_csq_entries(kept[0])] == ["TP53"]
    assert stats[0].removed == 1


def test_gene_id_filter_version_insensitive():
    lines = [
        _gene_record(1, [("TP53", "ENSG00000141510.17")]),
        _gene_record(2, [("EGFR", "ENSG00000146648")]),
    ]
    compiled = rf.compile_filters(
        [rf.ResultsFilter(field=rf.GENE_ID_FIELD, operator=rf.OPERATOR_IN, values=["ENSG00000141510"])],
        INDEX_MAP,
    )
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert len(kept) == 1
    gi = INDEX_MAP["Gene"]
    assert [e[gi] for e in rf.extract_csq_entries(kept[0])] == ["ENSG00000141510.17"]


def _transcript_group_filter(*values: str) -> rf.ResultsFilter:
    return rf.ResultsFilter(
        field=rf.TRANSCRIPT_GROUP_FIELD, operator=rf.OPERATOR_IN, values=list(values)
    )


def test_transcript_group_canonical():
    lines = [
        _group_record(
            1,
            [
                _group_entry("ENST_A", "YES", "", ""),  # canonical
                _group_entry("ENST_B", "", "", ""),  # neither
            ],
        )
    ]
    fi = GROUP_INDEX_MAP["Feature"]
    compiled = rf.compile_filters([_transcript_group_filter("canonical")], GROUP_INDEX_MAP)
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert len(kept) == 1
    assert [e[fi] for e in rf.extract_csq_entries(kept[0])] == ["ENST_A"]


def test_transcript_group_mane_any_of():
    lines = [
        _group_record(
            1,
            [
                _group_entry("ENST_A", "", "NM_1.1", ""),  # MANE Select (refseq present)
                _group_entry("ENST_B", "YES", "", ""),  # canonical only
                _group_entry("ENST_C", "", "", "NM_2.1"),  # MANE Plus Clinical
            ],
        )
    ]
    fi = GROUP_INDEX_MAP["Feature"]
    # select the two MANE groups; canonical-only transcript B is pruned
    compiled = rf.compile_filters(
        [_transcript_group_filter("mane_select", "mane_plus_clinical")], GROUP_INDEX_MAP
    )
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert len(kept) == 1
    assert [e[fi] for e in rf.extract_csq_entries(kept[0])] == ["ENST_A", "ENST_C"]


def test_transcript_group_rejects_unknown_group():
    with pytest.raises(rf.FilterError):
        rf.compile_filters([_transcript_group_filter("nonsense")], GROUP_INDEX_MAP)


# --- allele frequency ---------------------------------------------------------

AF_COLUMNS = (
    "Allele Consequence Feature "
    "gnomAD_exomes_AF gnomAD_exomes_AF_nfe AoU_gvs_all_af AoU_gvs_max_subpop"
).split()
AF_INDEX_MAP = {name: i for i, name in enumerate(AF_COLUMNS)}


def _af_entry(exomes: str, nfe: str, aou: str) -> str:
    return "|".join(["T", "missense_variant", "ENST_1", exomes, nfe, aou, "eur"])


def _af_record(pos: int, entries: list[str]) -> str:
    return f"chr1\t{100 + pos}\tid_{pos:02d}\tC\tT\t.\t.\tCSQ={','.join(entries)}\n"


def _af_filter(operator, threshold, match="any", values=None) -> rf.ResultsFilter:
    return rf.ResultsFilter(
        field=rf.ALLELE_FREQUENCY_FIELD,
        operator=operator,
        values=values or [],
        threshold=threshold,
        match=match,
    )


def test_af_columns_discovery_excludes_subpop_label():
    assert rf.af_columns(AF_INDEX_MAP) == [
        "gnomAD_exomes_AF",
        "gnomAD_exomes_AF_nfe",
        "AoU_gvs_all_af",
    ]


def test_af_source_descriptor():
    assert rf.af_source_descriptor("gnomAD_exomes_AF") == {
        "key": "gnomAD_exomes_AF",
        "source": "gnomad_exomes",
        "population": "",
    }
    assert rf.af_source_descriptor("gnomAD_genomes_AF_grpmax") == {
        "key": "gnomAD_genomes_AF_grpmax",
        "source": "gnomad_genomes",
        "population": "grpmax",
    }
    assert rf.af_source_descriptor("AoU_gvs_all_af") == {
        "key": "AoU_gvs_all_af",
        "source": "all_of_us",
        "population": "",
    }
    assert rf.af_source_descriptor("AoU_gvs_afr_af")["population"] == "afr"
    assert rf.af_source_descriptor("SYMBOL") is None


def test_af_le_any_keeps_when_one_meets():
    lines = [
        _af_record(1, [_af_entry("0.3", "0.01", "0.5")]),  # nfe 0.01 <= 0.05 -> keep
        _af_record(2, [_af_entry("0.3", "0.2", "0.5")]),  # none <= 0.05 -> drop
    ]
    compiled = rf.compile_filters([_af_filter("le", 0.05, "any")], AF_INDEX_MAP)
    kept, stats = rf.apply_filter_pipeline(lines, compiled)
    assert len(kept) == 1
    assert kept[0] == lines[0]
    assert stats[0].removed == 1


def test_af_le_all_requires_every_value():
    lines = [
        _af_record(1, [_af_entry("0.3", "0.01", "0.5")]),  # not all <= 0.05 -> drop
        _af_record(2, [_af_entry("0.01", "0.02", "0.03")]),  # all <= 0.05 -> keep
    ]
    compiled = rf.compile_filters([_af_filter("le", 0.05, "all")], AF_INDEX_MAP)
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert len(kept) == 1
    assert kept[0] == lines[1]


def test_af_ge_and_eq():
    line = _af_record(1, [_af_entry("0.3", "0.01", "0.5")])
    kept_ge, _ = rf.apply_filter_pipeline(
        [line], rf.compile_filters([_af_filter("ge", 0.4, "any")], AF_INDEX_MAP)
    )
    assert len(kept_ge) == 1  # aou 0.5 >= 0.4
    kept_eq, _ = rf.apply_filter_pipeline(
        [line], rf.compile_filters([_af_filter("eq", 0.01, "any")], AF_INDEX_MAP)
    )
    assert len(kept_eq) == 1  # nfe == 0.01


def test_af_specific_columns_only():
    lines = [_af_record(1, [_af_entry("0.3", "0.3", "0.5")])]
    # test only the exomes overall column; 0.3 > 0.05 -> drop
    compiled = rf.compile_filters(
        [_af_filter("le", 0.05, "any", values=["gnomAD_exomes_AF"])], AF_INDEX_MAP
    )
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert kept == []


def test_af_no_data_excluded_for_now():
    # all AF columns empty -> no data -> currently dropped (revisit)
    lines = [_af_record(1, [_af_entry("", "", "")])]
    compiled = rf.compile_filters([_af_filter("le", 0.05, "any")], AF_INDEX_MAP)
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert kept == []


def test_af_ignores_missing_but_tests_present():
    # exomes empty (ignored), nfe 0.01 present and <= 0.05 -> keep
    lines = [_af_record(1, [_af_entry("", "0.01", "")])]
    compiled = rf.compile_filters([_af_filter("le", 0.05, "any")], AF_INDEX_MAP)
    kept, _ = rf.apply_filter_pipeline(lines, compiled)
    assert len(kept) == 1


def test_af_rejects_in_operator():
    with pytest.raises(rf.FilterError):
        rf.compile_filters(
            [rf.ResultsFilter(field=rf.ALLELE_FREQUENCY_FIELD, operator="in", values=[])],
            AF_INDEX_MAP,
        )


def test_empty_values_is_noop():
    lines = [_record(1, ["missense_variant"])]
    compiled = rf.compile_filters([_consequence_filter()], INDEX_MAP)
    assert compiled == []
    kept, stats = rf.apply_filter_pipeline(lines, compiled)
    assert kept == lines
    assert stats == []


# --- end to end via get_results_from_path ------------------------------------


def test_get_results_filtered_totals_and_metadata(tmp_path):
    records = [
        _record(1, ["missense_variant"]),
        _record(2, ["synonymous_variant"]),
        _record(3, ["missense_variant"]),
        _record(4, ["intron_variant"]),
        _record(5, ["missense_variant"]),
    ]
    vcf_path = _write_vcf(tmp_path, records)

    result = get_results_from_path(
        page_size=10,
        page=1,
        vcf_path=FilePath(vcf_path),
        filters=[_consequence_filter("missense_variant")],
    )

    # Three missense records survive; pagination total reflects the filtered set.
    assert result.metadata.pagination.total == 3
    assert len(result.variants) == 3
    assert result.metadata.filters is not None
    assert result.metadata.filters.unfiltered_total == 5
    assert result.metadata.filters.filtered_total == 3
    assert result.metadata.filters.stats[0].field == "consequence"
    assert result.metadata.filters.stats[0].removed == 2


def test_get_results_prunes_nonmatching_transcripts(tmp_path):
    # One variant with two transcripts: one missense, one upstream. Filtering on
    # missense must keep the variant but drop the upstream transcript.
    csq = (
        "T|missense_variant|MODERATE|GENE|ENSG|Transcript|ENST_A|protein_coding,"
        "T|upstream_gene_variant|MODIFIER|GENE|ENSG|Transcript|ENST_B|protein_coding"
    )
    records = [f"chr1\t200\tv1\tC\tT\t.\t.\tCSQ={csq}\n"]
    vcf_path = _write_vcf(tmp_path, records)

    result = get_results_from_path(
        page_size=10,
        page=1,
        vcf_path=FilePath(vcf_path),
        filters=[_consequence_filter("missense_variant")],
    )

    assert len(result.variants) == 1
    all_consequences = [
        consequence
        for allele in result.variants[0].alternative_alleles
        for prediction in allele.predicted_molecular_consequences
        for consequence in prediction.consequences
    ]
    assert "missense_variant" in all_consequences
    assert "upstream_gene_variant" not in all_consequences


def test_get_results_filtered_pagination_slices(tmp_path):
    records = [_record(i, ["missense_variant"]) for i in range(1, 8)]
    # add some non-matching noise interleaved
    records += [_record(i, ["intron_variant"]) for i in range(8, 11)]
    vcf_path = _write_vcf(tmp_path, records)

    page1 = get_results_from_path(
        page_size=5,
        page=1,
        vcf_path=FilePath(vcf_path),
        filters=[_consequence_filter("missense_variant")],
    )
    page2 = get_results_from_path(
        page_size=5,
        page=2,
        vcf_path=FilePath(vcf_path),
        filters=[_consequence_filter("missense_variant")],
    )

    assert page1.metadata.pagination.total == 7
    assert len(page1.variants) == 5
    assert len(page2.variants) == 2  # remainder of the 7 filtered records
