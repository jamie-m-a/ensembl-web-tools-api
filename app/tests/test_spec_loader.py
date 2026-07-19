"""Tests for spec_loader: content-digest versioning, assembly -> merged-spec
resolution, and the sidecar that pins the merged document to a job.
"""

import json

import pytest
from pydantic import FilePath

from app.vep.utils.spec_loader import (
    EXPECTED_COLUMNS_SIDECAR_FILE,
    SPEC_SIDECAR_FILE,
    _content_digest,
    load_expected_columns_sidecar,
    load_merged_spec,
    load_merged_spec_file,
    load_spec_sidecar,
    resolve_merged_spec,
    write_expected_columns_sidecar,
    write_spec_sidecar,
)
from app.vep.utils.vcf_results import _load_pinned_spec

SAMPLE = {
    "genome": {"assembly": "GRCh38"},
    "config": {"entries": []},
    "parsing": {"plugins": []},
}


# --- content digest -----------------------------------------------------


def test_digest_is_independent_of_key_order():
    reordered = {
        "parsing": {"plugins": []},
        "config": {"entries": []},
        "genome": {"assembly": "GRCh38"},
    }
    assert _content_digest(SAMPLE) == _content_digest(reordered)


def test_digest_ignores_any_spec_version_already_present():
    """spec_version can't affect its own value -- it must be excluded before
    hashing, or the digest would depend on whatever was there before it."""
    with_version = {**SAMPLE, "spec_version": "sha256:whatever"}
    assert _content_digest(SAMPLE) == _content_digest(with_version)


def test_digest_changes_with_real_content():
    changed = {**SAMPLE, "genome": {"assembly": "GRCh37"}}
    assert _content_digest(SAMPLE) != _content_digest(changed)


# --- load_merged_spec_file: version is computed, not authored ------------


def test_load_merged_spec_file_computes_version_ignoring_file_placeholder(tmp_path):
    path = tmp_path / "spec.json"
    path.write_text(json.dumps({**SAMPLE, "spec_version": "sha256:0000"}))
    spec = load_merged_spec_file(path)
    assert spec.spec_version.startswith("sha256:")
    assert spec.spec_version != "sha256:0000"
    # deterministic: the same content computes the same digest on a second load
    assert load_merged_spec_file(path).spec_version == spec.spec_version
    # and the digest is mirrored onto the nested parsing view
    assert spec.parsing.spec_version == spec.spec_version


def test_load_merged_spec_file_works_with_no_version_in_the_file(tmp_path):
    """The bundled spec files don't carry a spec_version at all -- it is purely
    computed at load time."""
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(SAMPLE))
    spec = load_merged_spec_file(path)
    assert spec.spec_version.startswith("sha256:")


def test_bundled_human_grch38_spec_loads_and_has_a_real_digest():
    spec = load_merged_spec("human_grch38")
    assert spec.spec_version.startswith("sha256:")
    assert spec.spec_version != "sha256:" + "0" * 64
    assert len(spec.config_entries()) > 0
    assert len(spec.parse_plugins()) > 0


# --- resolve_merged_spec -------------------------------------------------


def test_resolve_grch38_returns_the_human_spec():
    resolved = resolve_merged_spec("GRCh38.p14")
    bundled = load_merged_spec("human_grch38")
    assert resolved.spec_version == bundled.spec_version
    assert {p.plugin for p in resolved.parse_plugins()} == {
        p.plugin for p in bundled.parse_plugins()
    }


def test_resolve_matches_by_prefix_not_exact_string():
    """Real assembly_name values carry a patch suffix, e.g. "GRCh38.p14" -- an
    exact-match lookup would never resolve anything."""
    assert (
        resolve_merged_spec("GRCh38.p14").spec_version
        == load_merged_spec("human_grch38").spec_version
    )


def test_resolve_unknown_assembly_raises():
    """No spec exists yet for other assemblies. This must fail loudly rather
    than silently falling back to human_grch38 for, say, a mouse submission."""
    with pytest.raises(ValueError, match="GRCh37"):
        resolve_merged_spec("GRCh37.p13")


def test_resolve_empty_assembly_raises():
    with pytest.raises(ValueError):
        resolve_merged_spec("")


# --- sidecar ---------------------------------------------------------------


def test_write_and_load_spec_sidecar_round_trip(tmp_path):
    """The digest must survive a write -> reload of the pinned document. Both
    spec_version fields are stamped when written; the loader must exclude them
    again before hashing or the reloaded digest would differ from the original
    (the bug this test guards)."""
    spec = load_merged_spec("human_grch38")
    written_path = write_spec_sidecar(tmp_path, spec)
    assert written_path == tmp_path / SPEC_SIDECAR_FILE
    assert written_path.exists()

    (tmp_path / "output.vcf.gz").write_bytes(b"")
    loaded = load_spec_sidecar(FilePath(tmp_path / "output.vcf.gz"))
    assert loaded is not None
    assert loaded.spec_version == spec.spec_version
    assert len(loaded.parse_plugins()) == len(spec.parse_plugins())
    assert len(loaded.config_entries()) == len(spec.config_entries())


def test_load_spec_sidecar_missing_is_none(tmp_path):
    (tmp_path / "output.vcf.gz").write_bytes(b"")
    assert load_spec_sidecar(FilePath(tmp_path / "output.vcf.gz")) is None


# --- expected-columns sidecar (the per-job missing-field check) --------------


def test_expected_columns_sidecar_round_trip(tmp_path):
    columns = {"REVEL", "ClinVar_CLNSIG", "gnomAD_exomes_AF"}
    written = write_expected_columns_sidecar(tmp_path, columns)
    assert written == tmp_path / EXPECTED_COLUMNS_SIDECAR_FILE
    (tmp_path / "output.vcf.gz").write_bytes(b"")
    loaded = load_expected_columns_sidecar(FilePath(tmp_path / "output.vcf.gz"))
    assert loaded == columns


def test_load_expected_columns_sidecar_missing_is_none(tmp_path):
    (tmp_path / "output.vcf.gz").write_bytes(b"")
    assert load_expected_columns_sidecar(FilePath(tmp_path / "output.vcf.gz")) is None


def test_write_spec_sidecar_overwrites_the_previous_one(tmp_path):
    """Matches the DUMP_INI dev harness: one job in flight at a time, so a new
    submission's spec replaces the last one rather than accumulating."""
    write_spec_sidecar(tmp_path, load_merged_spec("human_grch38"))
    write_spec_sidecar(tmp_path, load_merged_spec("human_grch38"))
    assert (tmp_path / SPEC_SIDECAR_FILE).exists()
    # still exactly one sidecar file, not two
    assert len(list(tmp_path.glob("*spec*"))) == 1


# --- _load_pinned_spec: the results-time seam (vcf_results) -----------------
# The defensive wrapper get_results_from_path uses to load the pinned spec at
# results time. It must never let a missing or corrupt pin break parsing, and it
# returns the parsing half of the merged document.


def test_load_pinned_spec_returns_the_sidecar_parsing_when_present(tmp_path):
    write_spec_sidecar(tmp_path, load_merged_spec("human_grch38"))
    (tmp_path / "output.vcf.gz").write_bytes(b"")
    spec = _load_pinned_spec(FilePath(tmp_path / "output.vcf.gz"))
    assert spec is not None
    assert spec.spec_version == load_merged_spec("human_grch38").spec_version
    assert len(spec.plugins) > 0


def test_load_pinned_spec_missing_sidecar_is_none(tmp_path):
    (tmp_path / "output.vcf.gz").write_bytes(b"")
    assert _load_pinned_spec(FilePath(tmp_path / "output.vcf.gz")) is None


def test_load_pinned_spec_unreadable_sidecar_is_none_not_raised(tmp_path):
    """A corrupt pin must fall back, not 500 the results endpoint."""
    (tmp_path / SPEC_SIDECAR_FILE).write_text("{ not valid json")
    (tmp_path / "output.vcf.gz").write_bytes(b"")
    assert _load_pinned_spec(FilePath(tmp_path / "output.vcf.gz")) is None
