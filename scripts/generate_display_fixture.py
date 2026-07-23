"""Regenerate the frontend display-spec fixture from the bundled merged spec.

The standalone-web-vep tests render against the *real* `display` payload the
tools API serves (`MergedSpec.display_payload`), captured in
`displaySpec.fixture.ts`. That file must stay byte-equal to the served payload,
so run this after any change to a genome spec's `display` section (and after any
change to the display-spec models' serialisation).

It rewrites only the JSON body of the fixture, preserving the file's licence
header, import and doc comment (everything up to `= `).

Usage (from the repo root, with a sibling standalone-web-vep checkout):

    PYTHONPATH=app .venv/bin/python scripts/generate_display_fixture.py

Pass an explicit fixture path as the first argument to override the default
(which assumes standalone-web-vep sits beside this repo).
"""

import json
import sys
from pathlib import Path

from vep.utils.spec_loader import load_merged_spec

GENOME = "human_grch38"
MARKER = "export const displaySpecFixture: DisplaySpec = "

# .../vep/ensembl-web-tools-api/scripts/this.py -> parents[2] == .../vep
_REPOS_DIR = Path(__file__).resolve().parents[2]
DEFAULT_FIXTURE = (
    _REPOS_DIR
    / "standalone-web-vep/src/content/app/tools/vep/views/vep-submission-results"
    / "components/vep-results-annotation-detail/displaySpec.fixture.ts"
)


def regenerate(fixture: Path) -> None:
    text = fixture.read_text()
    if MARKER not in text:
        raise SystemExit(f"marker {MARKER!r} not found in {fixture}")

    payload = load_merged_spec(GENOME).display_payload()
    if payload is None:
        raise SystemExit(f"{GENOME} spec has no display section")

    body = json.dumps(
        payload.model_dump(mode="json", by_alias=True),
        indent=2,
        ensure_ascii=False,
    )
    prefix = text[: text.index(MARKER) + len(MARKER)]
    fixture.write_text(prefix + body + ";\n")


def main() -> None:
    fixture = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_FIXTURE
    if not fixture.exists():
        raise SystemExit(
            f"fixture not found: {fixture}\n"
            "Pass the path to displaySpec.fixture.ts as the first argument."
        )
    regenerate(fixture)
    print(f"regenerated {fixture}")


if __name__ == "__main__":
    main()
