"""Tests for the runtime missing-expected-field check (design §6.2): at results
time, warn if the pipeline output header is missing a CSQ column the submitted
options required. Non-fatal — extras are ignored and a missing pin is a no-op.
"""

import gzip
import logging

from pydantic import FilePath

from app.vep.utils.spec_loader import write_expected_columns_sidecar
from app.vep.utils.vcf_results import _check_expected_columns, _read_csq_columns


def _write_vcf(path, csq_columns):
    """A minimal gzipped VCF whose CSQ header declares `csq_columns`."""
    description = (
        "Consequence annotations from Ensembl VEP. Format: " + "|".join(csq_columns)
    )
    lines = [
        "##fileformat=VCFv4.2",
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="{description}">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(lines) + "\n")


def _write_vcf_without_csq(path):
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


# --- _read_csq_columns ------------------------------------------------------


def test_read_csq_columns_returns_the_header_columns(tmp_path):
    vcf = tmp_path / "output.vcf.gz"
    _write_vcf(vcf, ["Allele", "REVEL", "gnomAD_exomes_AF"])
    assert _read_csq_columns(FilePath(vcf)) == {"Allele", "REVEL", "gnomAD_exomes_AF"}


def test_read_csq_columns_none_when_no_csq_header(tmp_path):
    vcf = tmp_path / "output.vcf.gz"
    _write_vcf_without_csq(vcf)
    assert _read_csq_columns(FilePath(vcf)) is None


# --- _check_expected_columns ------------------------------------------------


def test_present_columns_do_not_warn(tmp_path, caplog):
    vcf = tmp_path / "output.vcf.gz"
    _write_vcf(vcf, ["Allele", "REVEL", "ClinVar_CLNSIG"])
    write_expected_columns_sidecar(tmp_path, {"REVEL", "ClinVar_CLNSIG"})
    with caplog.at_level(logging.WARNING):
        _check_expected_columns(FilePath(vcf))
    assert caplog.text == ""


def test_missing_column_warns_and_names_it(tmp_path, caplog):
    vcf = tmp_path / "output.vcf.gz"
    _write_vcf(vcf, ["Allele", "REVEL"])  # ClinVar_CLNSIG dropped by the pipeline
    write_expected_columns_sidecar(tmp_path, {"REVEL", "ClinVar_CLNSIG"})
    with caplog.at_level(logging.WARNING):
        _check_expected_columns(FilePath(vcf))
    assert "ClinVar_CLNSIG" in caplog.text
    assert "REVEL" not in caplog.text.replace("ClinVar", "")  # only the missing one


def test_extra_columns_are_ignored(tmp_path, caplog):
    vcf = tmp_path / "output.vcf.gz"
    _write_vcf(vcf, ["Allele", "REVEL", "SOMETHING_EXTRA"])
    write_expected_columns_sidecar(tmp_path, {"REVEL"})
    with caplog.at_level(logging.WARNING):
        _check_expected_columns(FilePath(vcf))
    assert caplog.text == ""


def test_missing_sidecar_is_a_noop(tmp_path, caplog):
    vcf = tmp_path / "output.vcf.gz"
    _write_vcf(vcf, ["Allele", "REVEL"])  # no expected_columns.json written
    with caplog.at_level(logging.WARNING):
        _check_expected_columns(FilePath(vcf))
    assert caplog.text == ""


def test_no_csq_header_warns(tmp_path, caplog):
    vcf = tmp_path / "output.vcf.gz"
    _write_vcf_without_csq(vcf)
    write_expected_columns_sidecar(tmp_path, {"REVEL"})
    with caplog.at_level(logging.WARNING):
        _check_expected_columns(FilePath(vcf))
    assert "No CSQ header" in caplog.text
