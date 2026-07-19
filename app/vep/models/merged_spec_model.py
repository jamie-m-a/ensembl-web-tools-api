"""The merged annotation-spec document: config + parsing for one genome.

One document, one content digest, pinned per job (spec_loader.py). It joins the
two halves the annotation API will serve — the option→`config.ini` rules
(`config_spec_model.py`) and the CSQ parsing rules (`parsing_spec_model.py`) —
under a single `spec_version`, so a job's options and the parsing of its results
are provably the same ruleset (design §8).

The two halves live as sibling sections rather than one per-plugin entry: the
config-set and parse-set only partly overlap and do not align 1:1 (`eve` config
feeds both the `eve` and `popeve` parsers; `hgvs`+`hgvsg` feed one `hgvs` parser;
10 config options have no parser at all). The explicit config→parse relation is
carried on each config entry's `parsed_as`, and this model's `model_validator`
is the load-time **consistency check** (design §6.1) that guards it.

See docs/design/merged-annotation-spec.md.
"""

import logging

from pydantic import BaseModel, ConfigDict, model_validator

from vep.models.config_spec_model import (
    ConfigEntry,
    ConfigSpec,
    CustomEmitter,
    LiteralFields,
)
from vep.models.parsing_spec_model import ParsingSpec, PluginSpec


class MergedSpec(BaseModel):
    """Config + parsing for one genome, under one content digest."""

    model_config = ConfigDict(extra="forbid")

    # Computed by spec_loader from the document's content, not authored; mirrored
    # onto `parsing.spec_version` so the pinned parse view carries the same id.
    spec_version: str = ""
    genome: dict | None = None
    config: ConfigSpec
    parsing: ParsingSpec

    def config_entries(self) -> list[ConfigEntry]:
        return self.config.entries

    def parse_plugins(self) -> list[PluginSpec]:
        return self.parsing.plugins

    @model_validator(mode="after")
    def _config_parsing_consistent(self) -> "MergedSpec":
        """Config↔parsing consistency check (design §6.1), run at load time.

        - every `parsed_as` id must resolve to a real parse plugin (error);
        - a `custom` emitter's derived columns must line up with its mapped parse
          plugin's `csq_fields` — exact for literal fields (ClinVar), prefix-only
          for the combinatorial gnomAD/AoU builders whose per-ancestry columns
          are discovered by the parser's `from_pattern` (error);
        - `plugin`/`flag` emitters are presence-checked only, since VEP derives
          their CSQ column names internally and the config line never states them;
        - a parse plugin that no config entry points at is a soft warning (it can
          never run), not a failure.
        """
        parse_ids = {p.plugin for p in self.parsing.plugins}
        referenced: set[str] = set()
        errors: list[str] = []

        for entry in self.config.entries:
            for parse_id in entry.parsed_as:
                referenced.add(parse_id)
                if parse_id not in parse_ids:
                    errors.append(
                        f"config entry {entry.id!r} references unknown parse "
                        f"plugin {parse_id!r}"
                    )
            if isinstance(entry.config, CustomEmitter) and entry.parsed_as:
                errors += self._check_custom_columns(entry, entry.config)

        if errors:
            raise ValueError("config/parsing inconsistency: " + "; ".join(errors))

        orphans = parse_ids - referenced
        if orphans:
            logging.warning(
                "parse plugins with no config entry enabling them: %s",
                sorted(orphans),
            )
        return self

    def _check_custom_columns(
        self, entry: ConfigEntry, emitter: CustomEmitter
    ) -> list[str]:
        short_name = emitter.params.get("short_name")
        if not isinstance(short_name, str):
            # A non-literal short_name (by_assembly / from_option) can't be
            # resolved to column names statically; nothing to check.
            return []

        mapped = [p for p in self.parsing.plugins if p.plugin in entry.parsed_as]
        csq_fields = {field for plugin in mapped for field in plugin.csq_fields}

        if isinstance(emitter.fields, LiteralFields):
            return [
                f"custom entry {entry.id!r} emits column "
                f"{short_name}_{field!s} that no mapped parse plugin "
                f"{sorted(entry.parsed_as)} declares"
                for field in emitter.fields.literal
                if f"{short_name}_{field}" not in csq_fields
            ]

        # Builder-based (gnomAD / All of Us): the combinatorial per-ancestry
        # columns are discovered by the parser's `from_pattern`, not listed, so
        # only require that the short_name prefix aligns with a declared column.
        if not any(field.startswith(f"{short_name}_") for field in csq_fields):
            return [
                f"custom entry {entry.id!r} short_name {short_name!r} matches no "
                f"CSQ column of its mapped parse plugin(s) {sorted(entry.parsed_as)}"
            ]
        return []
