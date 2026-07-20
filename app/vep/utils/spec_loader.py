"""Where the merged annotation spec comes from, how it is pinned to a job, and
which one applies to a given submission.

This is the seam that keeps "local JSON file" vs "annotation API" from being a
decision we have to make yet. Today `load_merged_spec` reads a JSON document
shipped in `vep/specs/`; when the API exists, only its body changes — an HTTP GET
plus a cache keyed on the spec's content digest. Everything downstream takes a
validated `MergedSpec` either way (its `.config` half drives config.ini
generation at submission, its `.parsing` half parses the results).

The file is not a mock of the API: it is the same document the API will serve,
and the same one pinned alongside a job at submission time.
"""

import hashlib
import json
from pathlib import Path

from pydantic import FilePath

from vep.models.merged_spec_model import MergedSpec

SPEC_DIR = Path(__file__).resolve().parent.parent / "specs"

# Written alongside a job's config at submission time, so the spec used to
# generate its options is the one used to parse its results, even if the
# bundled spec changes in between (see resolve_merged_spec / write_spec_sidecar).
# The name is retained (results-meta / page-index sidecars sit beside it); its
# content is now the whole merged document, not just the parsing half.
SPEC_SIDECAR_FILE = "parsing_spec.json"

# The per-job CSQ columns the submitted options require, pinned beside the spec
# at submission and checked against the pipeline output header at results time
# (the runtime missing-expected-field check). Job-specific, so kept separate from
# the (assembly-generic) merged spec document.
EXPECTED_COLUMNS_SIDECAR_FILE = "expected_columns.json"


def _content_digest(payload: dict) -> str:
    """A stable digest of a spec's meaning, ignoring its own `spec_version`.

    Independent of key order and of whitespace. `payload` must be a *validated
    model's* dump (see load_merged_spec_file), not raw user-authored JSON:
    hand-written specs use aliases (`from`, `as`) and omit fields at their
    default, while a round-tripped `model_dump()` uses field names and fills in
    every default. Hashing the raw file directly would make the digest depend on
    which of those wrote it — exactly the instability version pinning exists to
    avoid.
    """
    content = {key: value for key, value in payload.items() if key != "spec_version"}
    canonical = json.dumps(content, sort_keys=True, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def load_merged_spec_file(path: Path) -> MergedSpec:
    """Parse and validate a merged spec document, raising if it does not conform
    (including the config↔parsing consistency check in MergedSpec).

    `spec_version` is computed here, not read from the file: the version is a
    property of the content. The digest is taken from the *validated model's*
    canonical dump, and the nested `parsing.spec_version` is excluded before
    hashing (the top-level one is stripped by `_content_digest`) so a round-trip
    through the sidecar — where both are stamped — hashes to the same value. The
    computed digest is stamped on both the document and its `.parsing` view, so
    the pinned parse half reports the unified id.
    """
    payload = json.loads(path.read_text())
    payload.setdefault("spec_version", "")  # satisfy the field before it's computed
    spec = MergedSpec.model_validate(payload)
    canonical_dump = spec.model_dump(
        mode="json", by_alias=True, exclude={"parsing": {"spec_version": True}}
    )
    digest = _content_digest(canonical_dump)
    spec.spec_version = digest
    spec.parsing.spec_version = digest
    return spec


def load_merged_spec(name: str) -> MergedSpec:
    """The named merged spec from the bundled spec directory, e.g. "human_grch38"."""
    return load_merged_spec_file(SPEC_DIR / f"{name}.json")


# Assembly-name prefixes, mirroring ConfigIniParams' own is_human_grch38 /
# is_human_grch37 / is_mouse_reference checks (pipeline_model.py) so a spec is
# picked using the same notion of "which genome is this" as the ini builder.
# Only GRCh38 has a spec today; a submission for any other assembly fails
# loudly here rather than being silently parsed with the wrong one.
_ASSEMBLY_SPECS = {
    "GRCh38": "human_grch38",
}


def resolve_merged_spec(assembly_name: str) -> MergedSpec:
    """The merged spec for a submission's assembly.

    Only `assembly_name` is available at submission time (it is a field on
    ConfigIniParams already); `species_taxonomy_id` is not — that is only ever
    sent to /form_config today. Real per-species branching (as opposed to
    per-assembly) would need that added to the submission contract.
    """
    for prefix, spec_name in _ASSEMBLY_SPECS.items():
        if (assembly_name or "").startswith(prefix):
            return load_merged_spec(spec_name)
    raise ValueError(f"No spec available for assembly {assembly_name!r}")


def write_spec_sidecar(directory: str | Path, spec: MergedSpec) -> Path:
    """Pin `spec` to a job by writing the whole merged document into the job's
    directory.

    In the real pipeline, `directory` is the job's own outdir, alongside its
    config.ini and (eventually) its output VCF, so `load_spec_sidecar` finds it
    from the results path with no other bookkeeping needed.

    In the DUMP_INI dev harness there is no per-job outdir — results are read
    from a fixed LOCAL_RESULTS_VCF path, decoupled from any submission_id — so
    `directory` there is DUMP_INI_DIR, the same directory the config dump goes
    into. A submission there overwrites the previous sidecar, which matches how
    that harness already works: one manually-run job at a time.
    """
    path = Path(directory) / SPEC_SIDECAR_FILE
    path.write_text(spec.model_dump_json())
    return path


def load_spec_sidecar(vcf_path: FilePath) -> MergedSpec | None:
    """The merged spec pinned alongside `vcf_path`'s directory, or None if there
    isn't one (e.g. output from before this existed). Keyed off the VCF path the
    same way results_meta.json and the page-index sidecar are, via `.with_name()`."""
    sidecar_path = vcf_path.with_name(SPEC_SIDECAR_FILE)
    if not sidecar_path.exists():
        return None
    return load_merged_spec_file(sidecar_path)


def write_expected_columns_sidecar(
    directory: str | Path, columns: set[str]
) -> Path:
    """Pin the CSQ columns this job's options require, beside its spec sidecar,
    for the results-time missing-expected-field check. Sorted for a stable file."""
    path = Path(directory) / EXPECTED_COLUMNS_SIDECAR_FILE
    path.write_text(json.dumps(sorted(columns)))
    return path


def load_expected_columns_sidecar(vcf_path: FilePath) -> set[str] | None:
    """The expected CSQ columns pinned alongside `vcf_path`, or None if there is
    no sidecar (output from before this existed). Keyed off the VCF path via
    `.with_name()`, like the spec and page-index sidecars."""
    sidecar_path = vcf_path.with_name(EXPECTED_COLUMNS_SIDECAR_FILE)
    if not sidecar_path.exists():
        return None
    return set(json.loads(sidecar_path.read_text()))
