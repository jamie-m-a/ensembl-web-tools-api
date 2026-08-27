"""Regenerate the frontend display-spec fixture from the bundled merged spec.

LOCAL DEV HARNESS (fork only) — upstream deleted this in 68898ea ("Remove dev
mode, cleanup"); see the LOCAL DEV HARNESS block in ``core/config.py``.

The frontend tests render against the *real* `display` payload the tools API
serves (`MergedSpec.display_payload`), captured in `displaySpec.fixture.ts`.
That file must stay equal to the served payload, so run this after any change
to a genome spec's `display` section (and after any change to the display-spec
models' serialisation).

It rewrites only the JSON body of the fixture, preserving the file's licence
header, import and doc comment (everything up to `= `).

This emits raw JSON, which is *not* how the fixture is committed — run Prettier
over it afterwards or the diff is unreadable and CI's format check fails:

    PYTHONPATH=app ../.venv/bin/python app/vep/scripts/generate_display_fixture.py
    cd ../ensembl-client && npx prettier --write <fixture path>

Pass an explicit fixture path as the first argument to override the default
(which assumes an ensembl-client checkout sits beside this repo).
"""

import json
import sys
from pathlib import Path

from vep.utils.spec_loader import load_merged_spec

GENOME = "human_grch38"
MARKER = "export const displaySpecFixture: DisplaySpec = "

# .../vep/ensembl-web-tools-api/app/vep/scripts/this.py -> parents[4] == .../vep
_REPOS_DIR = Path(__file__).resolve().parents[4]
DEFAULT_FIXTURE = (
    _REPOS_DIR
    / "ensembl-client/src/content/app/tools/vep/views/vep-submission-results"
    / "components/vep-results-annotation-detail/displaySpec.fixture.ts"
)


def regenerate(fixture: Path) -> None:
    text = fixture.read_text()
    if MARKER not in text:
        raise SystemExit(f"marker {MARKER!r} not found in {fixture}")

    payload = load_merged_spec(GENOME).display_payload()
    if payload is None:
        raise SystemExit(f"{GENOME} spec has no display section")

    # `exclude_none` so the fixture carries what the spec *says*. Without it
    # every value wrote out every field it could have, so a three-key cell
    # became eleven keys of mostly null — and adding a field to a model churned
    # the fixture everywhere, a frontend diff for a change no frontend can see.
    #
    # Not `exclude_defaults`, which also drops defaults that are *not* None and
    # that the frontend has no way to reconstruct: the house truncation
    # (`visible_count: 3`) is a backend default, and dropping it silently
    # stopped every long list truncating. An omitted key and a null read the
    # same to the frontend; an omitted key and a 3 do not.
    body = json.dumps(
        payload.model_dump(mode="json", by_alias=True, exclude_none=True),
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
