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

See app/vep/docs/design/spec-and-extension-guide.md.
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
    DisplayGroupBlock,
    DisplayListBlock,
    DisplayMapRowsBlock,
    DisplayPayload,
    DisplayRowsBlock,
    DisplaySpec,
    DisplayTableBlock,
)
from vep.models.option_help_model import HelpSpec
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


# --- Display format <-> parsing type compatibility -------------------------- #
#
# Each display `format` assumes a value of a particular *shape*, and applying it
# to the wrong shape crashes the frontend renderer: `num` calls `.toPrecision`
# (throws on a string), `join`/`humanize_join` call `.join`/`.map` (throw on a
# string / non-list). The parsing target a display field reads fixes that shape,
# so the mismatch is caught here at load time instead of at render. A shape is a
# small tag: ("scalar", "num"|"string"), ("list", "num"|"string"|"object"),
# ("dict",) or ("object",).

_NUMERIC_TYPES = frozenset({"float", "int"})


def _scalar_shape(value_type: str) -> tuple[str, str]:
    """A scalar value's shape from its parsing `type` (`raw` is source text, so
    string-like)."""
    return ("scalar", "num" if value_type in _NUMERIC_TYPES else "string")


def _target_shape(target) -> tuple[str, ...]:
    """The shape a display field gets when it reads a whole parsing target,
    derived from the target's transform (and, for scalars/lists, element type)."""
    transform = target.transform
    if transform in ("scalar", "first"):
        return _scalar_shape(target.type)
    if transform == "list":
        return ("list", _scalar_shape(target.type)[1])
    if transform in ("zip", "chunk", "records", "stack"):
        return ("list", "object")
    if transform == "positional":
        return ("list", "object") if target.wrap == "list" else ("object",)
    if transform == "regex":
        return ("list", "object") if target.each else ("object",)
    if transform in ("pattern_map", "key_value"):
        return ("dict",)
    return ("object",)  # unreachable for the current Transform set; be safe


def _element_shapes(list_target) -> dict[str, tuple[str, ...]]:
    """The shape of each field a list target's elements declare.

    A `stack` declares its fields per group rather than once: the same field
    name in every group, plus each group's constant tags (always strings)."""
    shapes = {
        field.field: _scalar_shape(field.type)
        for field in list_target.as_fields or []
    }
    for group in list_target.of or []:
        for field in group.as_fields:
            shapes[field.field] = _scalar_shape(field.type)
        for tag in group.const:
            shapes[tag] = ("scalar", "string")
    return shapes


def _format_suits_shape(fmt: str, shape: tuple[str, ...]) -> bool:
    """Whether `format` can be applied to `shape` without crashing / misreading.
    `text` (and any unlisted format) only stringifies, so it suits anything."""
    if fmt == "num":
        return shape == ("scalar", "num")
    if fmt in ("humanize", "phenotype", "humanize_terms"):
        return shape == ("scalar", "string")
    if fmt == "join":
        return shape[0] == "list" and shape[1] in ("string", "num")
    if fmt == "humanize_join":
        return shape == ("list", "string")
    if fmt == "count":
        return shape[0] == "list" or shape == ("scalar", "string")
    return True


_FORMAT_NEEDS = {
    "num": "a numeric field",
    "humanize": "a string field",
    "phenotype": "a string field",
    "humanize_terms": "a string field",
    "join": "a list of scalars",
    "humanize_join": "a list of strings",
    "count": "a list, or a delimited string",
}

_SHAPE_DESCRIPTIONS = {
    ("scalar", "num"): "a numeric field",
    ("scalar", "string"): "a string field",
    ("list", "num"): "a list of numbers",
    ("list", "string"): "a list of strings",
    ("list", "object"): "a list of objects",
    ("dict",): "a map",
    ("object",): "an object",
}


def _describe_shape(shape: tuple[str, ...]) -> str:
    return _SHAPE_DESCRIPTIONS.get(shape, str(shape))


def _compose_errors(oid: str, compose, target_of) -> list[str]:
    """A `with_score` value renders `num(score) (humanize(classification))`, so
    the score must be numeric and the classification a string — either wrong
    crashes exactly as a bad row `format` would."""
    errors: list[str] = []
    for ref, needed in (
        (compose.score, ("scalar", "num")),
        (compose.classification, ("scalar", "string")),
    ):
        target = target_of(ref)
        if target is not None and _target_shape(target) != needed:
            errors.append(
                f"display option {oid!r} uses {ref!r} in a with_score value, "
                f"which needs {_describe_shape(needed)}, but it is "
                f"{_describe_shape(_target_shape(target))}"
            )
    return errors


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
    # The help behind each option's (?) button. Optional for the same reason as
    # `display`: every spec pinned to a job before this section existed has no
    # `help` key and must still load.
    help: HelpSpec | None = None

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
            options=self.display.options,
            plugin_scopes=self.plugin_scopes(),
            rating_scales=self.display.rating_scales,
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
        # A forced option's columns are expected too — its config line is emitted,
        # so its output must be there (ClinVar rides in on Phenotypes this way).
        options = self.config.effective_options(options)
        expected: set[str] = set()
        for entry in self.config.entries:
            if not options.get(entry.id):
                continue
            if not entry.requirements_met(options):
                continue  # a sub-option entry contributes nothing without its parent
            emitter = entry.config
            if isinstance(emitter, CustomEmitter):
                short_name = emitter.params.get("short_name")
                if isinstance(short_name, str):
                    # VEP always emits the bare `short_name` match column for a
                    # custom (the matched record / overlap), so it is expected
                    # whenever the option is on — independent of `fields`. This is
                    # also what keeps an identity field read from it (gnomAD SV's
                    # id) from being nulled by the AF-column gate.
                    expected.add(short_name)
                    # A fields-less custom's other columns are source-derived
                    # (gff/bed overlap), not statically known, so it requires no
                    # more than the base.
                    if emitter.fields is not None:
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

        errors += self._check_display()
        errors += self._check_stars_scales()

        if errors:
            raise ValueError("config/parsing inconsistency: " + "; ".join(errors))

        orphans = parse_ids - referenced
        if orphans:
            logging.warning(
                "parse plugins with no config entry enabling them: %s",
                sorted(orphans),
            )
        return self

    def _list_element_fields(
        self,
    ) -> dict[str, dict[tuple[str, ...], dict[str, tuple[str, ...] | None]]]:
        """Per plugin, what an element of each of its lists carries — each field
        with its shape where that is known — keyed by the path to that list.

        A path rather than a name because lists nest: a ClinVar condition holds
        the classifications its submitters gave, and each of those holds the
        submissions it counted. `("conditions",)` is the target itself,
        `("conditions", "classifications", "submitters")` a submission three
        levels in. Only the top level used to be known, so a display ref into a
        nested list had nothing to check against and a typo rendered an empty
        cell.

        What a joined-in row carries is *derived*, not declared: a join attaches
        the source list's own rows unchanged (see `_apply_joins`), so declaring
        their fields would be a second copy to keep in step -- and the two
        declarations that existed had already drifted from the targets they
        described.

        Shape is None for a name the target lists but does not type -- the two
        come from different places, `item_fields` naming what an element carries
        and `as`/`of` typing what the transform builds.
        """
        by_plugin: dict[
            str, dict[tuple[str, ...], dict[str, tuple[str, ...] | None]]
        ] = {}
        for plugin in self.parsing.plugins:
            paths: dict[tuple[str, ...], dict[str, tuple[str, ...] | None]] = {}
            for target in plugin.targets:
                shapes = _element_shapes(target)
                paths[(target.field,)] = {
                    name: shapes.get(name) for name in target.item_fields or []
                }
            # In declaration order: a join may draw from a list an earlier one
            # enriched, and then the attached rows carry that enrichment too.
            for join in plugin.joins or []:
                into = paths.get((join.into,))
                if into is None:
                    continue
                if join.count_into:
                    into[join.count_into] = ("scalar", "num")
                if not join.as_field:
                    continue
                into[join.as_field] = ("list", "object")
                source = dict(paths.get((join.source,), {}))
                attached = (join.into, join.as_field)
                if not join.count_by:
                    paths[attached] = source
                    continue
                # A counted join attaches groups, not rows: each is the value
                # grouped on plus its count, and -- with `nest_as` -- the rows
                # behind it.
                paths[attached] = {
                    join.count_by: source.get(join.count_by),
                    "count": ("scalar", "num"),
                }
                if join.nest_as:
                    paths[attached][join.nest_as] = ("list", "object")
                    paths[attached + (join.nest_as,)] = source
            # A `collapse` post-op moves fields out of the row and into a
            # nested list of its own, so the display resolves them a level down
            # from where the target declared them.
            for operation in plugin.post_joins or []:
                if operation.op != "collapse":
                    continue
                row = paths.get((operation.target,))
                if row is None:
                    continue
                target = next(
                    (t for t in plugin.targets if t.field == operation.target),
                    None,
                )
                shapes = _element_shapes(target) if target is not None else {}
                gathered = operation.fields or []
                for field in gathered:
                    row.pop(field, None)
                row[operation.into] = ("list", "object")
                paths[(operation.target, operation.into)] = {
                    field: shapes.get(field) for field in gathered
                }
            by_plugin[plugin.plugin] = paths
        return by_plugin

    def _check_display(self) -> list[str]:
        """Display↔parsing consistency, in one walk of the block tree.

        Two things are checked per site, because they were two walks and a
        construct kept being added to one and forgotten in the other -- which is
        how a stacked row's cells came to be ref-checked but never type-checked.
        Anything reachable here is now reachable by both.

        **That the field exists**: resolve every field a display option
        reads against the parsing plugins and their declared targets — a fixed
        row's `<plugin>.<field>`, a block's `when` field, a list block's
        `<plugin>.<listField>`, and each list element's item-relative refs (label
        and cells) against that list's element fields — plus every block's
        `requires`. Groups are flattened by `iter_blocks`, so their sub-blocks and
        their own `when` are checked the same way.

        Element refs resolve by *path*, so a column's `items` and their `expand`
        are checked against the nested lists they actually read rather than
        against the table's own row (see `_list_element_fields`).

        **And that its `format` suits it**: a format assumes a shape, and
        applying it to the wrong one crashes the renderer (`num` ->
        `.toPrecision`, `join`/`humanize_join` -> `.join`/`.map`). Only refs
        that resolve are shape-checked; an unresolved one is already reported
        above rather than complained about twice."""
        if self.display is None:
            return []

        targets_by_plugin = {
            plugin.plugin: {t.field: t for t in plugin.targets}
            for plugin in self.parsing.plugins
        }
        # Fields a plugin builds from the variant rather than from a CSQ
        # column. They have no target, so a display row naming one has to be
        # recognised here or the consistency check reads it as a typo.
        variant_link_fields = {
            plugin.plugin: set(plugin.variant_links or {})
            for plugin in self.parsing.plugins
        }
        paths_by_plugin = self._list_element_fields()
        errors: list[str] = []

        def target_of(ref: str):
            plugin, _, field = ref.partition(".")
            return targets_by_plugin.get(plugin, {}).get(field)

        def field_error(option_id: str, plugin: str, field: str) -> str | None:
            if plugin not in targets_by_plugin:
                return (
                    f"display option {option_id!r} references unknown parse "
                    f"plugin {plugin!r}"
                )
            if field in variant_link_fields.get(plugin, set()):
                return None
            if field not in targets_by_plugin[plugin]:
                return (
                    f"display option {option_id!r} references field {field!r} "
                    f"that parse plugin {plugin!r} does not produce"
                )
            if targets_by_plugin[plugin][field].join_source:
                return (
                    f"display option {option_id!r} references {plugin}.{field}, "
                    "which is a join source and is dropped from the output; read "
                    "it through the list it was joined into"
                )
            return None

        def scalar_ref_error(option_id: str, ref: str) -> str | None:
            plugin, _, field = ref.partition(".")
            return field_error(option_id, plugin, field)

        def item_errors(
            oid: str, plugin: str, path: tuple[str, ...], refs
        ) -> list[str]:
            where = f"{plugin}.{'.'.join(path)}"
            fields = paths_by_plugin.get(plugin, {}).get(path)
            if fields is None:
                return [
                    f"display option {oid!r} reads items of {where}, which the "
                    "parse does not produce as a list"
                ]
            return [
                f"display option {oid!r} references item field {ref!r} not in "
                f"{where}"
                for ref in refs
                if ref not in fields
            ]

        def column_errors(
            oid: str, plugin: str, path: tuple[str, ...], items
        ) -> list[str]:
            """A column's nested `items`/`expand` refs, one list level at a
            time. Recursive because an expanded line may itself expand."""
            found = item_errors(oid, plugin, path, items.item_field_refs())
            if not items.expand or found:
                # A bad `expand.from` is already reported; descending through it
                # would only add a second complaint about the same typo.
                return found
            nested = path + (items.expand.source,)
            for cell in items.expand.cells:
                found += column_errors(oid, plugin, nested, cell)
            if items.expand.emphasis:
                found += item_errors(
                    oid, plugin, nested, [items.expand.emphasis.field]
                )
            return found

        def check(oid: str, ref: str, fmt: str, shape: tuple[str, ...]) -> None:
            if not _format_suits_shape(fmt, shape):
                errors.append(
                    f"display option {oid!r} formats {ref!r} as {fmt!r}, but "
                    f"that is {_describe_shape(shape)}; {fmt!r} needs "
                    f"{_FORMAT_NEEDS.get(fmt, 'a compatible value')}"
                )

        def check_at(
            oid: str, plugin: str, path: tuple[str, ...], field: str | None, fmt: str
        ) -> None:
            """The same, for a field of a list reached by path. A format written
            one or two lists in was checked nowhere before, so `humanize_terms`
            was type-checked nowhere at all."""
            if field is None:
                return
            shape = paths_by_plugin.get(plugin, {}).get(path, {}).get(field)
            if shape is not None:  # an unresolved ref is item_errors' complaint
                check(oid, f"{plugin}.{'.'.join(path)}.{field}", fmt, shape)

        def check_items(
            oid: str, plugin: str, path: tuple[str, ...], items
        ) -> None:
            """A column's `items`, and the detail they expand onto."""
            if items.format:
                check_at(oid, plugin, path, items.source, items.format)
            if not items.expand:
                return
            nested = path + (items.expand.source,)
            for cell in items.expand.cells:
                check_items(oid, plugin, nested, cell)

        for option in self.display.options:
            oid = option.option_id
            for block in option.iter_blocks():
                # `when` reads a scalar `<plugin>.<field>`, like a row's `from`.
                if block.when:
                    err = scalar_ref_error(oid, block.when.field_ref)
                    if err:
                        errors.append(err)
                # A group only carries `when`; its children are visited too.
                if isinstance(block, DisplayGroupBlock):
                    continue
                if block.requires and block.requires not in targets_by_plugin:
                    errors.append(
                        f"display option {oid!r} requires unknown parse plugin "
                        f"{block.requires!r}"
                    )
                if isinstance(block, DisplayListBlock):
                    # The list field itself must be a target; then each element's
                    # item-relative refs must be in that target's item_fields.
                    plugin, list_field = block.list_ref()
                    err = field_error(oid, plugin, list_field)
                    if err:
                        errors.append(err)
                        continue
                    refs = list(block.item.item_field_refs())
                    if block.group_by:
                        # The field the items group on is an item field too, and
                        # a typo there would silently collapse every item into
                        # one unnamed section rather than fail.
                        refs.append(block.group_by.field)
                    errors += item_errors(oid, plugin, (list_field,), refs)
                    label = block.item.label
                    if label and label.format:
                        check_at(oid, plugin, (list_field,), label.source, label.format)
                    for part in (block.item.cells or []) + (block.item.rows or []):
                        if part.format:
                            check_at(
                                oid, plugin, (list_field,), part.source, part.format
                            )
                elif isinstance(block, DisplayTableBlock):
                    if block.rows is not None:
                        # Fixed matrix: each row value is a scalar `plugin.field`,
                        # checked like a row's `from`.
                        for ref in block.matrix_value_refs():
                            err = scalar_ref_error(oid, ref)
                            if err:
                                errors.append(err)
                        # A value column's format applies to each row value it
                        # holds -- checked against that scalar field.
                        for ref, fmt in block.value_column_formats():
                            target = target_of(ref)
                            if target is not None:
                                check(oid, ref, fmt, _target_shape(target))
                        continue
                    # List mode: iterates a list target like a list block; each
                    # column reads one of that target's item_fields.
                    plugin, list_field = block.list_ref()
                    err = field_error(oid, plugin, list_field)
                    if err:
                        errors.append(err)
                        continue
                    errors += item_errors(
                        oid, plugin, (list_field,), block.column_field_refs()
                    )
                    for column in block.columns:
                        if column.format:
                            check_at(
                                oid,
                                plugin,
                                (list_field,),
                                column.source,
                                column.format,
                            )
                        if column.items and column.source:
                            nested = (list_field, column.source)
                            errors += column_errors(oid, plugin, nested, column.items)
                            check_items(oid, plugin, nested, column.items)
                elif isinstance(block, DisplayMapRowsBlock):
                    # Every ref is a scalar-or-dict `plugin.field`; the rows
                    # themselves come from a vocabulary the response ships, so
                    # there is nothing here to check them against — a wrong
                    # `vocabulary`/`scope` yields no rows, which the equivalence
                    # check catches rather than the spec loader.
                    for ref in block.field_refs():
                        err = scalar_ref_error(oid, ref)
                        if err:
                            errors.append(err)
                elif isinstance(block, DisplayRowsBlock):
                    for row in block.rows:
                        for ref in row.field_refs():
                            err = scalar_ref_error(oid, ref)
                            if err:
                                errors.append(err)
                        if row.source and row.format:
                            target = target_of(row.source)
                            if target is not None:
                                check(
                                    oid, row.source, row.format, _target_shape(target)
                                )
                        if row.compose:
                            errors += _compose_errors(oid, row.compose, target_of)
                        # A row that stacks a list reads that list's element
                        # fields, exactly as a list block's item does.
                        list_ref = row.list_ref()
                        if list_ref is None:
                            continue
                        plugin, list_field = list_ref
                        if field_error(oid, plugin, list_field):
                            continue  # already reported by field_refs above
                        refs = list(row.item.item_field_refs())
                        if row.where:
                            # The table's identical `where` was checked and the
                            # row's was not; a typo here keeps every element
                            # rather than the ones meant, so a list split
                            # between two places shows up whole in both.
                            refs.append(row.where.field)
                        errors += item_errors(oid, plugin, (list_field,), refs)
                        for cell in row.item.cells or []:
                            if cell.format:
                                check_at(
                                    oid,
                                    plugin,
                                    (list_field,),
                                    cell.source,
                                    cell.format,
                                )
        return errors

    def _check_stars_scales(self) -> list[str]:
        """A `stars_from` cell names a *field* whose values name scales, so the
        display alone cannot tell whether those scales exist.

        Where the field is filled from a stack's `const` the values are stated
        in the parse spec, so they can be checked after all: ClinVar's aggregate
        classification is rated `clinvar_aggregate` for a germline one and
        `clinvar_somatic` for a somatic one, and neither name appears anywhere
        the display-side check can see. A field filled from the data is left
        alone -- an unknown scale there is a source adding wording, which
        correctly shows no stars.
        """
        if self.display is None:
            return []
        targets = {
            (plugin.plugin, target.field): target
            for plugin in self.parsing.plugins
            for target in plugin.targets
        }

        def const_values(target, field: str) -> set[str] | None:
            groups = getattr(target, "of", None) or []
            values = set()
            for group in groups:
                value = (group.const or {}).get(field)
                if not isinstance(value, str):
                    return None  # data, not a constant
                values.add(value)
            return values or None

        errors: list[str] = []
        for option in self.display.options:
            for block in option.iter_blocks():
                for row in getattr(block, "rows", None) or []:
                    list_ref = getattr(row, "list_ref", None)
                    list_ref = list_ref() if callable(list_ref) else None
                    if list_ref is None or row.item is None:
                        continue
                    target = targets.get(list_ref)
                    if target is None:
                        continue
                    for cell in row.item.cells:
                        if not cell.stars_from:
                            continue
                        values = const_values(target, cell.stars_from)
                        for value in sorted(values or ()):
                            if value not in self.display.rating_scales:
                                errors.append(
                                    f"display option {option.option_id!r} rates "
                                    f"on scale {value!r}, which "
                                    f"{'.'.join(list_ref)}.{cell.stars_from} "
                                    "states but `rating_scales` does not define"
                                )
        return errors

    def _check_custom_columns(
        self, entry: ConfigEntry, emitter: CustomEmitter
    ) -> list[str]:
        short_name = emitter.params.get("short_name")
        if not isinstance(short_name, str):
            # A non-literal short_name (by_assembly / from_option) can't be
            # resolved to column names statically; nothing to check.
            return []
        if emitter.fields is None:
            # A fields-less overlap custom emits source-derived columns VEP names
            # itself; nothing to check statically.
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
