"""Turn selected options + a ConfigSpec into VEP config.ini lines.

The declarative counterpart to
`pipeline_model.ConfigIniParams.create_config_ini_file`: given the options a
submission selected, this emits the same `plugin …`, `custom …` and flag lines
the hardcoded builder does. Additive for now — built alongside the old path and
proved equal by `tests/test_config_interpreter.py` (a differential test over
option combinations) before it replaces it.

NOT here: the always-on base config (`force_overwrite`, `numbers`, `symbol`, …
and the assembly-gated `mane`/`assembly`). Those are VEP-invocation invariants
that stay in the backend next to the per-genome `gff`/`fasta` resolution — this
module only turns *selected options* into lines. See docs/design/.
"""

from vep.models.config_spec_model import (
    AllofusPopulationFields,
    ByAssembly,
    ConfigSpec,
    CustomEmitter,
    FlagEmitter,
    FromOption,
    GnomadAncestrySexFields,
    LiteralFields,
    PluginEmitter,
)


# Returned by _param_value for an assembly-conditional param that has no value
# on this assembly (SpliceAI's snv_ensembl on GRCh37): the kwarg is dropped.
_SKIP = object()


def _interpolate(text: str, context: dict[str, str]) -> str:
    """Substitute the backend-provided `{path}` / `{gff}` tokens. Anything else
    (notably the pipeline's `###CHR###` per-chromosome placeholder) is left
    untouched."""
    for token, value in context.items():
        text = text.replace("{" + token + "}", value)
    return text


def _param_value(param, options: dict, assembly: str, context: dict) -> str:
    if isinstance(param, str):
        return _interpolate(param, context)
    if isinstance(param, ByAssembly):
        if assembly in param.by_assembly:
            return _interpolate(param.by_assembly[assembly], context)
        if param.omit_if_absent:
            return _SKIP
        return _interpolate(param.by_assembly["GRCh38"], context)
    if isinstance(param, FromOption):
        value = options.get(param.from_option)
        if param.equals is not None:
            return "1" if value == param.equals else "0"
        if param.as_type == "value":
            return str(value)
        return str(int(bool(value)))
    raise ValueError(f"unknown param value: {param!r}")


def _params_str(params: dict, options, assembly, context) -> list[str]:
    parts = []
    for key, value in params.items():
        resolved = _param_value(value, options, assembly, context)
        if resolved is not _SKIP:
            parts.append(f"{key}={resolved}")
    return parts


def _variadic_suffix(flags, options: dict) -> str:
    """IntAct's selected sub-flags: `,all=1` when every one is on, else
    `,<kw>=1` for each selected, else nothing."""
    selected = [f for f in flags.options if options.get(f.option)]
    if flags.all_shortcut and len(selected) == len(flags.options):
        return f",{flags.all_shortcut}=1"
    if selected:
        return "".join(f",{f.keyword}=1" for f in selected)
    return ""


def build_fields(fields, options: dict) -> list[str]:
    if isinstance(fields, LiteralFields):
        return list(fields.literal)

    if isinstance(fields, GnomadAncestrySexFields):
        # non_ukb is inserted after base when the UK-Biobank toggle is off
        # (exomes only; genomes has no include_ukb_option).
        non_ukb = bool(fields.include_ukb_option) and not options.get(
            fields.include_ukb_option
        )

        def _code(anc_code: str, sex_code: str) -> str:
            parts = [fields.base]
            if non_ukb:
                parts.append("non_ukb")
            if anc_code:
                parts.append(anc_code)
            if sex_code:
                parts.append(sex_code)
            return "_".join(parts)

        result: list[str] = []
        for ancestry in fields.ancestries:
            if not options.get(ancestry.option):
                continue
            if not ancestry.sex_split:  # grpmax: one field, no XX/XY
                result.append(_code(ancestry.code, ""))
                continue
            for sex in fields.sexes:
                if options.get(f"{ancestry.option}_{sex.suffix}"):
                    result.append(_code(ancestry.code, sex.code))
        return result

    if isinstance(fields, AllofusPopulationFields):
        result = []
        for population in fields.populations:
            if options.get(population.option):
                result.extend(population.codes)
        return result

    raise ValueError(f"unknown field builder: {fields!r}")


def _emit_entry(entry, options, assembly, context) -> str | None:
    emitter = entry.config

    if isinstance(emitter, FlagEmitter):
        # Flags are always written (as 0 or 1), like the base flag block.
        return f"{emitter.keyword} {int(bool(options.get(entry.id)))}"

    if not options.get(entry.id):
        return None  # plugin / custom lines only appear when the option is on

    if isinstance(emitter, PluginEmitter):
        parts = _params_str(emitter.params, options, assembly, context)
        line = f"plugin {emitter.name}"
        if parts:
            line += "," + ",".join(parts)
        if emitter.flags:
            line += _variadic_suffix(emitter.flags, options)
        return line

    if isinstance(emitter, CustomEmitter):
        # A fields-less custom (gff/bed overlap) writes no `fields=` clause at
        # all — VEP emits the source's attributes itself.
        field_list = build_fields(emitter.fields, options) if emitter.fields else []
        if emitter.omit_if_no_fields and not field_list:
            return None
        join = getattr(emitter.fields, "join", "%")
        parts: list[str] = []
        for key, value in emitter.params.items():
            resolved = _param_value(value, options, assembly, context)
            if resolved is not _SKIP:
                parts.append(f"{key}={resolved}")
            if emitter.fields is not None and key == emitter.fields_after:
                parts.append(f"fields={join.join(field_list)}")
        return "custom " + ",".join(parts)

    raise ValueError(f"unknown emitter: {emitter!r}")


def emit_config_lines(
    spec: ConfigSpec,
    options: dict,
    *,
    assembly: str,
    plugin_path: str,
    gff: str,
) -> list[str]:
    """The option-driven config.ini lines for `options`, in entry `order`.

    `options` is a flat {option_id: value} map (a `ConfigIniParams.model_dump()`).
    `plugin_path`/`gff` fill the `{path}`/`{gff}` tokens.
    """
    context = {"path": plugin_path, "gff": gff}
    # A selected option can force other options on for config emission only
    # (ProtVar needs HGVSg computed to build its link). This is confined to the
    # config lines — it never touches the options the results view gates display
    # on — so a forced flag is computed without adding its row.
    effective = dict(options)
    for entry in spec.entries:
        if options.get(entry.id):
            for forced_id in entry.forces_on:
                effective[forced_id] = True
    lines: list[str] = []
    for entry in sorted(spec.entries, key=lambda e: e.order):
        line = _emit_entry(entry, effective, assembly, context)
        if line is not None:
            lines.append(line)
    return lines
