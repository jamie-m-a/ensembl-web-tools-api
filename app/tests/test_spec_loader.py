"""Tests for spec_loader: content-digest versioning, assembly -> spec
resolution, and the sidecar that pins a spec to a job.
"""

import json

import pytest
from pydantic import FilePath

from app.vep.utils.spec_loader import (
    SPEC_DIR,
    SPEC_SIDECAR_FILE,
    _content_digest,
    load_spec,
    load_spec_file,
    load_spec_sidecar,
    resolve_spec,
    write_spec_sidecar,
)

SAMPLE = {"genome": {"assembly": "GRCh38"}, "plugins": []}


# --- content digest -----------------------------------------------------


def test_digest_is_independent_of_key_order():
    reordered = {"plugins": [], "genome": {"assembly": "GRCh38"}}
    assert _content_digest(SAMPLE) == _content_digest(reordered)


def test_digest_ignores_any_spec_version_already_present():
    """spec_version can't affect its own value -- it must be excluded before
    hashing, or the digest would depend on whatever was there before it."""
    with_version = {**SAMPLE, "spec_version": "sha256:whatever"}
    assert _content_digest(SAMPLE) == _content_digest(with_version)


def test_digest_changes_with_real_content():
    changed = {"genome": {"assembly": "GRCh37"}, "plugins": []}
    assert _content_digest(SAMPLE) != _content_digest(changed)


# --- load_spec_file: version is computed, not authored -------------------


def test_load_spec_file_computes_version_ignoring_file_placeholder(tmp_path):
    path = tmp_path / "spec.json"
    path.write_text(json.dumps({**SAMPLE, "spec_version": "sha256:0000"}))
    spec = load_spec_file(path)
    assert spec.spec_version == _content_digest(SAMPLE)
    assert spec.spec_version != "sha256:0000"


def test_load_spec_file_works_with_no_version_in_the_file(tmp_path):
    """The bundled spec files don't carry a spec_version at all -- it is
    purely computed at load time."""
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(SAMPLE))
    spec = load_spec_file(path)
    assert spec.spec_version == _content_digest(SAMPLE)


def test_bundled_human_grch38_spec_loads_and_has_a_real_digest():
    spec = load_spec("human_grch38")
    assert spec.spec_version.startswith("sha256:")
    assert spec.spec_version != "sha256:" + "0" * 64
    assert len(spec.plugins) > 0


# --- resolve_spec ---------------------------------------------------------


def test_resolve_spec_grch38_returns_the_human_spec():
    resolved = resolve_spec("GRCh38.p14")
    bundled = load_spec("human_grch38")
    assert resolved.spec_version == bundled.spec_version
    assert {p.plugin for p in resolved.plugins} == {p.plugin for p in bundled.plugins}


def test_resolve_spec_matches_by_prefix_not_exact_string():
    """Real assembly_name values carry a patch suffix, e.g. "GRCh38.p14" -- an
    exact-match lookup would never resolve anything."""
    assert resolve_spec("GRCh38.p14").spec_version == load_spec("human_grch38").spec_version


def test_resolve_spec_unknown_assembly_raises():
    """No spec exists yet for other assemblies. This must fail loudly rather
    than silently falling back to human_grch38 for, say, a mouse submission."""
    with pytest.raises(ValueError, match="GRCh37"):
        resolve_spec("GRCh37.p13")


def test_resolve_spec_empty_assembly_raises():
    with pytest.raises(ValueError):
        resolve_spec("")


# --- sidecar ---------------------------------------------------------------


def test_write_and_load_spec_sidecar_round_trip(tmp_path):
    spec = load_spec("human_grch38")
    written_path = write_spec_sidecar(tmp_path, spec)
    assert written_path == tmp_path / SPEC_SIDECAR_FILE
    assert written_path.exists()

    (tmp_path / "output.vcf.gz").write_bytes(b"")
    loaded = load_spec_sidecar(FilePath(tmp_path / "output.vcf.gz"))
    assert loaded is not None
    assert loaded.spec_version == spec.spec_version
    assert len(loaded.plugins) == len(spec.plugins)


def test_load_spec_sidecar_missing_is_none(tmp_path):
    (tmp_path / "output.vcf.gz").write_bytes(b"")
    assert load_spec_sidecar(FilePath(tmp_path / "output.vcf.gz")) is None


def test_write_spec_sidecar_overwrites_the_previous_one(tmp_path):
    """Matches the DUMP_INI dev harness: one job in flight at a time, so a new
    submission's spec replaces the last one rather than accumulating."""
    write_spec_sidecar(tmp_path, load_spec("human_grch38"))
    first_mtime = (tmp_path / SPEC_SIDECAR_FILE).stat().st_mtime_ns
    write_spec_sidecar(tmp_path, load_spec("human_grch38"))
    assert (tmp_path / SPEC_SIDECAR_FILE).exists()
    # still exactly one sidecar file, not two
    assert len(list(tmp_path.glob("*spec*"))) == 1
