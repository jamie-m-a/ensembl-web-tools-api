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
    FlagEmitter,
    FromOption,
    LiteralFields,
    PluginEmitter,
)
from vep.models.display_spec_model import (
    DisplayListBlock,
    DisplayPayload,
    DisplaySpec,
)
from vep.models.parsing_spec_model import ParsingSpec, PluginSpec
from vep.utils.config_interpreter import build_fields


def _is_simple_plugin(emitter: PluginEmitter) -> bool:
    """A plugin whose CSQ columns all appear in the header whenever it runs — no
    sub-option gates one away. That means no variadic `flags` (IntAct) and no
    `from_option` params (ProtVar/mutfunc/DosageSensitivity sub-flags), so all of
    its `csq_fields` are safe to *require*. Sub-flagged plugins are excluded from
    the expected set: turning a sub-flag off legitimately drops its column, and
    requiring it anyway would false-positive."""
    return emitter.flags is None and not any(
        isinstance(value, FromOption) for value in emitter.params.values()
    )


class MergedSpec(BaseModel):
    """Config + parsing for one genome, under one content digest."""

    model_config = ConfigDict(extra="forbid")

    # Computed by spec_loader from the document's content, not authored; mirrored
    # onto `parsing.spec_version` so the pinned parse view carries the same id.
    spec_version: str = ""
    genome: dict | None = None
    config: ConfigSpec
    parsing: ParsingSpec
    # How the parsed annotations are laid out in the results detail. Optional,
    # and the default matters: every spec pinned to a job before this section
    # existed has no `display` key and must still load (the results path then
    # falls back to the current genome's display spec).
    display: DisplaySpec | None = None

    def config_entries(self) -> list[ConfigEntry]:
        return self.config.entries

    def plugin_scopes(self) -> dict[str, str]:
        """plugin id -> "allele" | "transcript", derived from `parsing`.

        The display spec's rows name a plugin but deliberately do not say which
        entity it hangs off; that is a property of the parser and is stated
        exactly once, here. Authoring it a second time in `display` would create
        precisely the hand-synced seam the merged document exists to remove.
        """
        return {plugin.plugin: plugin.scope for plugin in self.parsing.plugins}

    def display_payload(self) -> DisplayPayload | None:
        """The display spec plus its derived scopes, as served on the results
        response. None when this document has no display section."""
        if self.display is None:
            return None
        return DisplayPayload(
            options=self.display.options, plugin_scopes=self.plugin_scopes()
        )

    def parse_plugins(self) -> list[PluginSpec]:
        return self.parsing.plugins

    def expected_csq_columns(self, options: dict) -> set[str]:
        """The CSQ columns a job with these selected `options` must have in its
        output header (design §6.2) — the per-job basis for the runtime
        missing-expected-field check. The union, over *enabled* config entries, of:

          * custom emitters → the exact columns the emitted `fields=` names
            (`<short_name>_<field>`), including the combinatorial gnomAD/AoU set,
            derived from the *same* `build_fields` that wrote the config line;
          * simple plugin emitters (no column-gating sub-flags) → their mapped
            parse plugins' `csq_fields`;
          * flag emitters → only the *allele-scoped* mapped parse plugins'
            `csq_fields`. A flag can emit conditional columns (HGVSc/HGVSp exist
            only where a variant has transcript context), so transcript-scoped
            flag columns are not required; but an allele-scoped one is present for
            every variant (HGVSg whenever `--hgvsg` is on, SPDI whenever `--spdi`
            is), so it is safe to require.

        Sub-flagged plugins and transcript-scoped flag columns are deliberately
        excluded — a sub-option (or the absence of a transcript) can legitimately
        drop one of their columns (see `_is_simple_plugin`). Extras are never
        required. gnomAD/AoU with nothing selected emit no line and so contribute
        nothing, matching the config.
        """
        by_plugin = {p.plugin: p for p in self.parsing.plugins}
        expected: set[str] = set()
        for entry in self.config.entries:
            if not options.get(entry.id):
                continue
            emitter = entry.config
            if isinstance(emitter, CustomEmitter):
                short_name = emitter.params.get("short_name")
                if isinstance(short_name, str):
                    for field in build_fields(emitter.fields, options):
                        expected.add(f"{short_name}_{field}")
            elif isinstance(emitter, PluginEmitter) and _is_simple_plugin(emitter):
                for parse_id in entry.parsed_as:
                    plugin = by_plugin.get(parse_id)
                    if plugin is not None:
                        expected.update(plugin.csq_fields)
            elif isinstance(emitter, FlagEmitter):
                for parse_id in entry.parsed_as:
                    plugin = by_plugin.get(parse_id)
                    if plugin is not None and plugin.scope == "allele":
                        expected.update(plugin.csq_fields)
        return expected

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
          never run), not a failure;
        - every display row's `from`/`compose` field reference must resolve to a
          real parse plugin and one of its declared target fields (error) — the
          display-side analogue of the above, and the main guard against the
          laid-out labels drifting from what the parsers actually produce.
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

        errors += self._check_display_refs()

        if errors:
            raise ValueError("config/parsing inconsistency: " + "; ".join(errors))

        orphans = parse_ids - referenced
        if orphans:
            logging.warning(
                "parse plugins with no config entry enabling them: %s",
                sorted(orphans),
            )
        return self

    def _check_display_refs(self) -> list[str]:
        """Display↔parsing consistency: resolve every field a display option
        reads against the parsing plugins and their declared targets — a fixed
        row's `<plugin>.<field>`, a list block's `<plugin>.<listField>`, and each
        list cell's item-relative field against that list target's `item_fields`
        — plus every block's `requires`."""
        if self.display is None:
            return []

        targets_by_plugin = {
            plugin.plugin: {t.field: t for t in plugin.targets}
            for plugin in self.parsing.plugins
        }
        errors: list[str] = []

        def field_error(option_id: str, plugin: str, field: str) -> str | None:
            if plugin not in targets_by_plugin:
                return (
                    f"display option {option_id!r} references unknown parse "
                    f"plugin {plugin!r}"
                )
            if field not in targets_by_plugin[plugin]:
                return (
                    f"display option {option_id!r} references field {field!r} "
                    f"that parse plugin {plugin!r} does not produce"
                )
            return None

        for option in self.display.options:
            oid = option.option_id
            for block in option.blocks:
                if block.requires and block.requires not in targets_by_plugin:
                    errors.append(
                        f"display option {oid!r} requires unknown parse plugin "
                        f"{block.requires!r}"
                    )
                if not isinstance(block, DisplayListBlock):
                    continue
                # The list field itself must be a target; then each cell's
                # item-relative refs must be in that target's item_fields.
                plugin, list_field = block.list_ref()
                err = field_error(oid, plugin, list_field)
                if err:
                    errors.append(err)
                    continue
                item_fields = set(
                    targets_by_plugin[plugin][list_field].item_fields or []
                )
                for cell in block.item.cells:
                    for item_field in cell.item_field_refs():
                        if item_field not in item_fields:
                            errors.append(
                                f"display option {oid!r} list "
                                f"{plugin}.{list_field} references item field "
                                f"{item_field!r} not in its target's item_fields"
                            )
            for plugin, field in option.scalar_field_refs():
                err = field_error(oid, plugin, field)
                if err:
                    errors.append(err)
        return errors

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
