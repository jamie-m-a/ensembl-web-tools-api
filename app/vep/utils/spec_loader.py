"""Where a parsing spec comes from.

This is the seam that keeps "local JSON file" vs "annotation API" from being a
decision we have to make yet. Today `load_spec` reads a JSON document shipped in
`vep/parsing_specs/`; when the API exists, only this function's body changes — an
HTTP GET plus a cache keyed on the spec's content digest. Everything downstream
takes a validated `ParsingSpec` either way.

The file is not a mock of the API: it is the same document the API will serve,
and the same one pinned alongside a job at submission time.
"""

import json
from pathlib import Path

from vep.models.parsing_spec_model import ParsingSpec

SPEC_DIR = Path(__file__).resolve().parent.parent / "parsing_specs"


def load_spec_file(path: Path) -> ParsingSpec:
    """Parse and validate a spec document. Raises if it does not conform."""
    return ParsingSpec.model_validate_json(path.read_text())


def load_spec(name: str) -> ParsingSpec:
    """The named spec from the bundled spec directory, e.g. "human_grch38"."""
    return load_spec_file(SPEC_DIR / f"{name}.json")
