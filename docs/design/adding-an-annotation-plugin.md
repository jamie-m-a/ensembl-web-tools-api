# Enabling a new annotation plugin — handover

One annotation is described by **one merged JSON document per genome** — three
sibling sections that own *generating its config*, *parsing its CSQ columns*, and
*laying out its results*, under a single content digest pinned to every job.

This is the practical companion to [`merged-annotation-spec.md`](./merged-annotation-spec.md):
what to load into that document, the grammar available to each section, and the
few seams that are still code. The grammar below is drawn verbatim from the
models — `config_spec_model.py`, `parsing_spec_model.py`, `display_spec_model.py`,
`merged_spec_model.py`, `spec_loader.py`, and the frontend `vepDisplaySpec.ts`.

---

## 1. What is data, what is still code

Adding a plugin used to mean editing four hand-synced places. Now the
data-driven parts collapse into the spec JSON the API serves; a small residue
stays in the backend because it is genuine logic or an invocation invariant.

**JSON — loaded to the API** (`app/vep/specs/<genome>.json`):

- **`config`** — how the option becomes a `config.ini` line (flag / plugin / custom).
- **`parsing`** — how its CSQ columns become structured data (9 transforms).
- **`display`** — how that data is laid out in the results detail (rows / blocks / formats).

**Backend — still code:**

- **Panel activation + the option's form definition** — `form_panels.py`
  (species/assembly gating is logic, not data).
- **A boolean field per option** — `ConfigIniParams` in `pipeline_model.py`.
- **The plugin data file** behind the `{path}` token.
- **A frontend override** — *only* if the display is interactive/derived.

The one document is resolved by assembly prefix in `spec_loader._ASSEMBLY_SPECS`.
It is served at submission to build the config, pinned as a `parsing_spec.json`
sidecar, and re-served at results to parse + lay out. **Change any section and
the content digest moves** — so any pinned `dev-data` sidecar must be
regenerated.

---

## 2. The per-plugin skeleton

A new option adds *up to* three entries — one per section. A config-only flag
(e.g. `spdi`) needs only the first; an unparsed option (e.g. `geno2mp`) skips
parsing and display. Fill the blanks:

```jsonc
// app/vep/specs/<genome>.json — the three sibling sections
{
  "config": { "entries": [
    {
      "id": "my_option",            // == the ConfigIniParams field / form option id
      "order": 30,                  // position in the emitted config.ini
      "parsed_as": ["my_plugin"],   // → parsing plugin id(s); [] if unparsed
      "config": { "emit": "plugin", "name": "MyPlugin",
                  "params": { "file": "{path}/my_data.vcf.gz" } }
    }
  ]},
  "parsing": { "plugins": [
    {
      "plugin": "my_plugin", "scope": "transcript", "output": "my_plugin",
      "csq_fields": ["MyPlugin_score"],   // header cols it owns; none present ⇒ didn't run
      "targets": [
        { "field": "score", "from": "MyPlugin_score", "transform": "scalar", "type": "float" }
      ]
    }
  ]},
  "display": { "options": [
    {
      "option_id": "my_option",
      "blocks": [ { "kind": "rows", "heading": "My plugin", "rows": [
        { "label": "Score", "from": "my_plugin.score", "format": "num" }
      ] } ]
    }
  ]}
}
```

**Two joins are load-time checked** (`MergedSpec` model validator):
`config.parsed_as` must name a real `parsing` plugin, and every display
`from`/`compose` must resolve to a real plugin **and** one of its declared target
fields. `scope` (allele vs transcript) is stated **once** on the parser and
derived for display (`plugin_scopes`) — never authored twice.

---

## 3. `config` — emitting the config line

A closed set of three emitters (`config.entries[].config`), no open escape hatch.

| emit | shape | produces / when to use |
|------|-------|------------------------|
| **`flag`** | `{ keyword }` | → `keyword 1\|0`. A VEP base-output switch with no data file — `hgvs`, `hgvsg`, `spdi`, `protein`. |
| **`plugin`** | `{ name, params, flags?, when? }` | → `plugin Name,k=v,…`. The common case: a VEP plugin with a data file and/or sub-flags. |
| **`custom`** | `{ params, fields, fields_after?, omit_if_no_fields?, when? }` | → `custom file=…,short_name=…,fields=…,format=…`. A VCF/BED overlay — gnomAD, All of Us, ClinVar. |

**Param values** (plugin & custom `params`):

| value | shape | meaning |
|-------|-------|---------|
| literal | `"…"` | A string. Tokens `{path}` (plugin-data root) and `{gff}` (per-genome gff) are interpolated by the backend. |
| by_assembly | `{ by_assembly, omit_if_absent? }` | Assembly-keyed file (`GRCh38`/`GRCh37`). Falls back to GRCh38 unless `omit_if_absent` drops the whole param. |
| from_option | `{ from_option, as, equals? }` | `as:"int"` → 1/0 from a boolean sub-option; `+ equals` → 1 when a *select* equals the value; `as:"value"` → the sub-option's value verbatim. |

- **Variadic sub-flags — `flags`** (IntAct-style): append `,keyword=1` per
  selected sub-option, or a single `,all_shortcut=1` when every one is on.
- **Assembly gate — `when`**: `{ assembly: ["GRCh38"] }` emits the line only for
  matching assemblies (prefix match), inside a single multi-assembly spec.

**Field builders** (custom `fields`):

| builder | for | notes |
|---------|-----|-------|
| `literal` | A fixed list | `{ literal: ["CLNSIG", …] }` — ClinVar. |
| `gnomad_ancestry_sex` | Combinatorial AF | ancestry × sex × UK-Biobank subset → `AF[_non_ukb][_<anc>][_XX\|_XY]`. Codes live on the option defs. |
| `allofus_populations` | Flat population AF | Concatenate each selected population's field code(s); no sex split. |

`omit_if_no_fields` drops the whole custom line when nothing is selected.
`fields_after` (default `short_name`) controls where `fields=` lands in the arg
order.

---

## 4. `parsing` — reading the CSQ columns

A `PluginSpec` declares `scope` (`allele` | `transcript`), `output`, the
`csq_fields` it owns, and its `targets`. Two "nothing here" guards mirror the old
hand-written parsers: `csq_fields` absent from the header ⇒ the plugin didn't
run; `require_any_input` / `require_any_output` ⇒ present but empty ⇒ no
annotation.

Each target (`targets[].transform`) is a small transform over column(s):

| transform | key params | column(s) → value |
|-----------|-----------|-------------------|
| **`scalar`** | `from, type` | One column → one value. |
| **`list`** | `from` | One column → `&`-split list; empties and `NA` dropped. |
| **`first`** | `from` | One column → first real item of a `&`-split list. |
| **`zip`** | `from[], as[], align, drop_when?, post?` | N aligned `&`-lists → list of objects, positions preserved. `align: max\|min`. |
| **`regex`** | `from, pattern, as[], each?` | Named groups → object(s). `each` applies per `&`-item; non-matches skipped. |
| **`pattern_map`** | `from_pattern, exclude?` | Columns matching `"X_{ph}"` → dict keyed by the wildcard. Discovered from the header — the field set need not be known up front (gnomAD per-ancestry AF). |
| **`chunk`** | `from, size, as[]` | Take `size` `&`-items per object → list (ProtVar partner & score pairs). |
| **`positional`** | `from, as[], wrap?` | Assign `&`-items to `as` strictly by index; extras ignored, missing null. `wrap:"list"` for a single-element list. |
| **`key_value`** | `from, pair_delimiter, kv_delimiter` | Split into a dict — order-independent (UTRAnnotator's unstable key order). |

**Shared modifiers** on any target:

- `type` = `string \| float \| int \| raw` (`raw` keeps source text verbatim).
- `when` `{ field, includes }` gates the build on another column's `&`-membership.
- `drop_when` `{ all_null \| null:<field> }` discards produced elements.
- `post` `[{ op: dedup \| sort, by, desc, nulls }]` reshapes the list.
- Sub-fields (`as[]` / regex groups) are `FieldSpec`: `{ field, type, replace?, strip? }`.

---

## 5. `display` — laying out the result

Two mutually-exclusive routes, chosen by the shape of the output:

**Declarative DSL — zero frontend code.** Use when the output is *a heading and
some label/value rows*. Authored entirely in JSON; the generic renderer walks it.

- a **row**: `label` + `from: "plugin.field"` *or* `compose`
- a **block**: `heading?`, `requires?`, `rows[]`
- an **option**: a *sequence* of blocks (EVE = a bare row plus a sibling popEVE block)

**Frontend override — a renderer case.** Use when the output is *conditional,
derived or interactive* — no declarative row can express it: shape-flips on the
data (ClinVar), derived counts (IntAct), popups / external link-outs / truncated
lists / imported summaries (ProtVar, OpenTargets, GO, phenotypes…). Register the
id in `OVERRIDE_OPTION_IDS` and add a `case` in `VepResultsAnnotationDetail`.

Display row fields:

| field | values | note |
|-------|--------|------|
| `from` | `"plugin.field"` | The value's source — mutually exclusive with `compose`. |
| `compose` | `{ format:"with_score", classification, score }` | One string from two fields — "Likely benign (0.07)". Drops the row if the classification is absent. |
| `format` | `text · num · humanize · phenotype · join` | Which frontend formatter renders the value. `text` is the default. |
| `placeholder` | `"—"` | Unset ⇒ an absent value **drops the row**; set ⇒ keeps it and shows this (SpliceAI's eight deltas always read as a set). |
| `mono` / `help` / `key` | `bool / str / str` | Monospace value · a (?) tooltip beside the label · an explicit React key. |

A block's `requires: "plugin"` makes the whole block vanish when that plugin
produced nothing — needed only where placeholder rows would otherwise render as a
wall of dashes.

---

## 6. Enablement checklist

An ordered sequence — the joins fail loudly at load time if a step is skipped.

1. **Place the plugin data file** behind the `{path}` token.
   *(Note: `PLUGIN_PATH` in `pipeline_model.py` is a pre-production placeholder —
   real per-genome resolution is still to be wired.)*
2. **Define the form option** — add it to a panel with its label / default /
   sub-options and species-assembly gating. `form_panels.py`
3. **Add the `ConfigIniParams` field** — a boolean per option (plus any
   sub-option fields), default matching the form. A test asserts every panel
   option id is a field (`test_form_panels.py`). `pipeline_model.py`
4. **Author the `config` entry** — the emitter (flag / plugin / custom) and
   `parsed_as`. `specs/<genome>.json`
5. **Author the `parsing` plugin** — scope, `csq_fields`, targets. Only if the
   plugin emits columns to read.
6. **Author the `display` — or add an override** — declarative blocks for the
   simple case; otherwise register the id and add a renderer `case` (+ any typed
   accessor).
7. **Load & check** — loading the spec runs the config↔parsing and
   display→parsing consistency checks. Run both suites (backend pytest, frontend
   vitest + tsc).
8. **Regenerate the dev-data sidecar** — the content digest moved, so any pinned
   `parsing_spec.json` sidecar is stale; regenerate before testing the results path.
9. **Differential-test against a real VCF** — the validation that consistently
   pays off: run the spec over a real dev-data VCF carrying the columns, not just
   fixtures.

---

**A note on the asymmetry.** The config-set and parse-set only partly overlap: a
config option can have no parser (`parsed_as: []` — flags like `spdi`/`protein`,
and `geno2mp`, `gnomad_mt`, the nearest-* options), and a parser can attach where
no display row reads it. That asymmetry is the normal case; the consistency
checks in `merged_spec_model.py` police it (a `parsed_as` pointing at a missing
plugin is an error; a parse plugin no config entry enables is a soft warning).
