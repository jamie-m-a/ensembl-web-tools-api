"""Static, strongly-typed model of the *config-generation* half of the merged
annotation spec.

Sibling to `parsing_spec_model.py`. Where the parsing spec says how to read a
plugin's CSQ columns, this says how a *selected option* becomes a line in the VEP
`config.ini` — replacing the hardcoded `PLUGIN_CONFIG_LINES` /
`PLUGIN_CONFIG_LINES_BY_ASSEMBLY` maps and the `create_config_ini_file` body in
`pipeline_model.py`. It is *data* (the `config` section of the merged JSON under
`specs/`, later served by the annotation API), so it is validated hard on arrival
(`extra="forbid"`).

The emitters are a **small closed set** — `{flag, plugin, custom}` — derived by
enumerating what `create_config_ini_file` actually emits, not invented, exactly
as the parsing transforms were. The always-on base config (`force_overwrite`,
`numbers`, `symbol`, … and the assembly-gated `mane`/`assembly`) is deliberately
NOT here: it is a VEP-invocation invariant that stays in the backend, next to the
per-genome `gff`/`fasta` resolution. See docs/design/merged-annotation-spec.md.
"""

from typing import Annotated, Literal, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator


# --------------------------------------------------------------------------- #
# Param values — the right-hand side of one plugin/custom kwarg.               #
# A value is a literal string, or one of these small typed forms.             #
# --------------------------------------------------------------------------- #

class ByAssembly(BaseModel):
    """A value chosen by the submission's assembly, e.g. an assembly-specific
    data file. Keys are assembly prefixes (`GRCh38`/`GRCh37`); the interpreter
    falls back to `GRCh38` when the assembly isn't listed, matching the ini
    builder's `by_assembly.get(assembly, by_assembly["GRCh38"])`.
    """

    model_config = ConfigDict(extra="forbid")

    by_assembly: dict[str, str]
    # When the submission's assembly isn't a key: fall back to `GRCh38` (mirrors
    # the ini builder's `by_assembly.get(assembly, by_assembly["GRCh38"])`)
    # unless this is set — then the whole param is dropped. For a param that
    # exists on some assemblies only (SpliceAI's `snv_ensembl`, GRCh38 only).
    omit_if_absent: bool = False


class FromOption(BaseModel):
    """A flag derived from another option. Without `equals`: `int(bool(option))`
    (ProtVar's stability/pocket/int, dosage's cover, mutfunc's sub-flags).
    With `equals`: 1 when a *select* option equals the value, else 0 (the
    TSSDistance direction radio → three upstream/downstream/both flags).
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    from_option: str
    equals: str | None = None
    # `as` is a Python keyword.
    #   int    -> int(bool(option)), or 1 when `equals` matches (the 0/1 flags)
    #   value  -> the option's own value verbatim (NearestExonJB's max_range)
    as_type: Literal["int", "value"] = Field(alias="as")


# A literal wins the union unambiguously; the two dict forms discriminate on
# their distinct required keys (`by_assembly` vs `from_option`).
ParamValue = Union[str, ByAssembly, FromOption]


# --------------------------------------------------------------------------- #
# Variadic sub-flags (IntAct) and the genome gate.                            #
# --------------------------------------------------------------------------- #

class VariadicFlag(BaseModel):
    model_config = ConfigDict(extra="forbid")

    option: str   # the boolean sub-option id
    keyword: str  # the ini flag it appends when on, e.g. `feature_ac`


class VariadicFlags(BaseModel):
    """IntAct-style: append `,<keyword>=1` for each selected sub-option, or a
    single `,<all_shortcut>=1` when *all* of them are selected. None selected
    leaves the base line untouched.
    """

    model_config = ConfigDict(extra="forbid")

    options: list[VariadicFlag]
    all_shortcut: str | None = None


class GenomeGate(BaseModel):
    """Emit only for these assemblies (prefix match, like the ini builder's
    `is_human_grch38` etc.). Mirrors the parser's CSQ `when`, on the genome."""

    model_config = ConfigDict(extra="forbid")

    assembly: list[str]


# --------------------------------------------------------------------------- #
# `fields=` for custom lines — a small closed set of named builders over the   #
# open field-code data that lives on the option definitions.                   #
# --------------------------------------------------------------------------- #

class LiteralFields(BaseModel):
    """A fixed field list, e.g. ClinVar's `CLNSIG%CLNSIGCONF`."""

    model_config = ConfigDict(extra="forbid")

    literal: list[str]


class AncestryCode(BaseModel):
    """One gnomAD ancestry option and its field-code component. `code` is empty
    for "all". `sex_split=False` marks a code that takes no XX/XY suffix and is a
    plain toggle (genomes' `grpmax`)."""

    model_config = ConfigDict(extra="forbid")

    option: str
    code: str
    sex_split: bool = True


class SexCode(BaseModel):
    """A sex sub-option suffix and its field-code component (both="", female=XX,
    male=XY). The sub-option id is `<ancestry.option>_<suffix>`."""

    model_config = ConfigDict(extra="forbid")

    suffix: str
    code: str


class PopulationCode(BaseModel):
    """One All of Us population option and the field code(s) it contributes
    ("max" contributes two)."""

    model_config = ConfigDict(extra="forbid")

    option: str
    codes: list[str]


class GnomadAncestrySexFields(BaseModel):
    """The gnomAD exomes/genomes grammar: for each selected ancestry, for each
    selected sex-of-that-ancestry, emit `<base>[_non_ukb][_<anc>][_<XX|XY>]`.

    The builder is only the combinatorial algorithm; the ancestry/sex codes are
    open data. TODO (at merge): move `ancestries`/`sexes` onto the option
    definitions (`field_code` per Q1) and reference them here rather than
    inlining — carried here for now so the config interpreter is self-contained.
    """

    model_config = ConfigDict(extra="forbid")

    builder: Literal["gnomad_ancestry_sex"]
    base: str = "AF"
    # A boolean option that, when *false*, inserts `non_ukb` after `base`
    # (exomes only; genomes has no UK Biobank subset so it omits this).
    include_ukb_option: str | None = None
    join: str = "%"
    ancestries: list[AncestryCode]
    sexes: list[SexCode]


class AllofusPopulationFields(BaseModel):
    """All of Us: concatenate the codes of each selected population (no sex
    split; "max" contributes two). Same TODO as above — codes move to the option
    defs at merge."""

    model_config = ConfigDict(extra="forbid")

    builder: Literal["allofus_populations"]
    join: str = "%"
    populations: list[PopulationCode]


FieldsSpec = Union[LiteralFields, GnomadAncestrySexFields, AllofusPopulationFields]


# --------------------------------------------------------------------------- #
# The three emitters.                                                          #
# --------------------------------------------------------------------------- #

class FlagEmitter(BaseModel):
    """`<keyword> {0|1}` from the entry's own boolean option (hgvs, hgvsg, spdi,
    protein)."""

    model_config = ConfigDict(extra="forbid")

    emit: Literal["flag"]
    keyword: str


class PluginEmitter(BaseModel):
    """`plugin <name>,<k>=<v>,…` when the entry's option is on. Static params,
    assembly-keyed files, and sub-flag interpolation are all `ParamValue`s;
    IntAct's variadic sub-flags use `flags`.
    """

    model_config = ConfigDict(extra="forbid")

    emit: Literal["plugin"]
    name: str
    params: dict[str, ParamValue] = {}
    flags: VariadicFlags | None = None
    when: GenomeGate | None = None


class CustomEmitter(BaseModel):
    """`custom file=…,short_name=…,fields=…,format=…` — gnomAD/AoU/ClinVar. When
    `omit_if_no_fields`, the whole line is dropped if `fields` resolves empty
    (nothing selected)."""

    model_config = ConfigDict(extra="forbid")

    emit: Literal["custom"]
    params: dict[str, ParamValue] = {}
    fields: FieldsSpec
    # `fields=` is emitted immediately after this param, matching the arg order
    # VEP's `custom` lines use (…,short_name=…,fields=…,format=…).
    fields_after: str = "short_name"
    omit_if_no_fields: bool = False
    when: GenomeGate | None = None


ConfigEmitter = Annotated[
    Union[FlagEmitter, PluginEmitter, CustomEmitter], Field(discriminator="emit")
]


# --------------------------------------------------------------------------- #
# Document.                                                                    #
# --------------------------------------------------------------------------- #

class ConfigEntry(BaseModel):
    """One option's config rule. `id` matches the ConfigIniParams field / form
    option id, so a selected option finds its emitter. `order` is the position
    the line takes in the generated ini (the current builder's emission order is
    load-bearing for the golden-file tests, so it is explicit here)."""

    model_config = ConfigDict(extra="forbid")

    id: str
    order: int
    # The parse-plugin id(s) this option's output is parsed by — the explicit
    # config→parsing link the consistency check uses (merged_spec_model.py).
    # Empty for config-only options (flags like spdi/protein, and loeuf,
    # geno2mp, the nearest-* plugins). The
    # relation is not 1:1: one config may feed several parse entries
    # (eve → eve + popeve) and several configs may feed one ({hgvs, hgvsg} →
    # hgvs). Kept on the config side so the parsing specs stay untouched; it is
    # also the seed of a future per-entry merge (design §3).
    parsed_as: list[str] = []
    # Other option ids to treat as on for config emission whenever this option
    # is selected — a config-only dependency. ProtVar reads HGVSg to build its
    # link, so `protvar` forces `hgvsg` to be computed; this never touches the
    # user's own HGVSg selection, which is what the results view gates the HGVSg
    # row's display on, so the value is computed without showing the row.
    forces_on: list[str] = []
    config: ConfigEmitter


class ConfigSpec(BaseModel):
    """The config half of the merged document, for one genome."""

    model_config = ConfigDict(extra="forbid")

    genome: dict | None = None
    entries: list[ConfigEntry]

    @model_validator(mode="after")
    def _unique_ids(self) -> "ConfigSpec":
        ids = [entry.id for entry in self.entries]
        dupes = {i for i in ids if ids.count(i) > 1}
        if dupes:
            raise ValueError(f"duplicate config entry ids: {sorted(dupes)}")
        return self

    def entry(self, option_id: str) -> ConfigEntry | None:
        return next((e for e in self.entries if e.id == option_id), None)
