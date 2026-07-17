"""Where a parsing spec comes from, how it is pinned to a job, and which one
applies to a given submission.

This is the seam that keeps "local JSON file" vs "annotation API" from being a
decision we have to make yet. Today `load_spec` reads a JSON document shipped in
`vep/parsing_specs/`; when the API exists, only its body changes — an HTTP GET
plus a cache keyed on the spec's content digest. Everything downstream takes a
validated `ParsingSpec` either way.

The file is not a mock of the API: it is the same document the API will serve,
and the same one pinned alongside a job at submission time.
"""

import hashlib
import json
from pathlib import Path

from pydantic import FilePath

from vep.models.parsing_spec_model import ParsingSpec

SPEC_DIR = Path(__file__).resolve().parent.parent / "parsing_specs"

# Written alongside a job's config at submission time, so the spec used to
# generate its options is the one used to parse its results, even if the
# bundled spec changes in between (see resolve_spec / write_spec_sidecar).
SPEC_SIDECAR_FILE = "parsing_spec.json"


def _content_digest(payload: dict) -> str:
    """A stable digest of a spec's meaning, ignoring its own `spec_version`.

    Independent of key order and of whitespace. `payload` must be a *validated
    model's* dump (see load_spec_file), not raw user-authored JSON: hand-written
    specs use aliases (`from`, `as`) and omit fields at their default, while a
    round-tripped `model_dump()` uses field names and fills in every default.
    Hashing the raw file directly would make the digest depend on which of
    those wrote it — exactly the instability version pinning exists to avoid.
    """
    content = {key: value for key, value in payload.items() if key != "spec_version"}
    canonical = json.dumps(content, sort_keys=True, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def load_spec_file(path: Path) -> ParsingSpec:
    """Parse and validate a spec document, raising if it does not conform.

    `spec_version` is computed here, not read from the file: the version is a
    property of the content, not something to hand-author and risk going stale.
    Any `spec_version` present in the file is ignored. The digest is taken from
    the *validated model's* canonical dump, not the raw file, so it is the same
    whether the spec was hand-written (aliases, defaults omitted) or came back
    from a round-trip through model_dump_json (field names, defaults filled in)
    — both produce the same spec and must produce the same digest.
    """
    payload = json.loads(path.read_text())
    payload.setdefault("spec_version", "")  # satisfy the required field for now
    spec = ParsingSpec.model_validate(payload)
    canonical_dump = spec.model_dump(mode="json", by_alias=True)
    spec.spec_version = _content_digest(canonical_dump)
    return spec


def load_spec(name: str) -> ParsingSpec:
    """The named spec from the bundled spec directory, e.g. "human_grch38"."""
    return load_spec_file(SPEC_DIR / f"{name}.json")


# Assembly-name prefixes, mirroring ConfigIniParams' own is_human_grch38 /
# is_human_grch37 / is_mouse_reference checks (pipeline_model.py) so a spec is
# picked using the same notion of "which genome is this" as the ini builder.
# Only GRCh38 has a spec today; a submission for any other assembly fails
# loudly here rather than being silently parsed with the wrong one.
_ASSEMBLY_SPECS = {
    "GRCh38": "human_grch38",
}


def resolve_spec(assembly_name: str) -> ParsingSpec:
    """The parsing spec for a submission's assembly.

    Only `assembly_name` is available at submission time (it is a field on
    ConfigIniParams already); `species_taxonomy_id` is not — that is only ever
    sent to /form_config today. Real per-species branching (as opposed to
    per-assembly) would need that added to the submission contract.
    """
    for prefix, spec_name in _ASSEMBLY_SPECS.items():
        if (assembly_name or "").startswith(prefix):
            return load_spec(spec_name)
    raise ValueError(f"No parsing spec available for assembly {assembly_name!r}")


def write_spec_sidecar(directory: str | Path, spec: ParsingSpec) -> Path:
    """Pin `spec` to a job by writing it into the job's directory.

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


def load_spec_sidecar(vcf_path: FilePath) -> ParsingSpec | None:
    """The spec pinned alongside `vcf_path`'s directory, or None if there isn't
    one (e.g. output from before this existed). Keyed off the VCF path the same
    way results_meta.json and the page-index sidecar are, via `.with_name()`."""
    sidecar_path = vcf_path.with_name(SPEC_SIDECAR_FILE)
    if not sidecar_path.exists():
        return None
    return load_spec_file(sidecar_path)
