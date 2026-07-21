# Merged annotation spec — design

> **Status: IMPLEMENTED & MERGED (2026-07-21).** This started as a design draft
> (2026-07-18); everything it proposes has since shipped to `main` in both VEP
> repos across a series of merged PRs — the config→JSON migration, the go-flat
> results cutover, the display DSL (results-panel pinning + the per-option
> display spec), and the popup-link templates. It is kept as the original design
> rationale; where it and the shipped code disagree, **the code is
> authoritative**. A few forward-looking sections below (notably §9) have been
> refreshed in place. Genuinely-open follow-ups that remain: the AF
> `kind: "allele_frequency"` marker, and converting the still-unspecced
> SIFT/PolyPhen / uniprot / protein-matches tail (blocked on sample data).

Covers both VEP repos: `ensembl-web-tools-api` (backend) and `standalone-web-vep`
(frontend).

This document is the design for collapsing the drift-prone, hand-synced chain
that spanned four places into a single data-driven document served by the new
annotation API.

---

## 1. Problem & goal

Adding or changing one annotation plugin/option today means editing four places,
kept in sync by hand:

1. `app/vep/form_panels.py::get_visible_panels` — which panels/options show per
   species (served by `GET /form_config/{genome_id}`).
2. `app/vep/models/pipeline_model.py::ConfigIniParams` — option ids → the VEP
   `config.ini` lines → which plugins the pipeline runs (built at
   `POST /submissions`).
3. `app/vep/utils/vcf_results.py` `_parse_*` bank → parsed JSON models (run at
   `GET /submissions/{id}/results`); now largely superseded by the declarative
   parsing spec in `app/vep/parsing_specs/human_grch38.json` +
   `app/vep/utils/spec_interpreter.py`.
4. Frontend display — `VepResultsAnnotationDetail.tsx` (+ shared summary
   components) render each plugin bespoke-ly.

**Goal:** one **merged JSON document**, served by the new API and pinned per job,
that owns the option→config translation, the parsing rules, and the display
rules. The backend keeps only what is genuinely logic or invocation-invariant
(panel *activation*, the always-on base config, per-genome file resolution).

**What this buys us:** adding a plugin becomes (mostly) one JSON entry; the
front/back contract is enforced by content-pinning + a load-time consistency
check; and unknown plugins degrade gracefully instead of being invisible.

---

## 2. Decisions (settled in review, 2026-07-18)

- **The new API owns the option→`config.ini` mapping**, not just parse+display.
- **One merged document**, one content-derived `spec_version` digest — the
  existing sidecar/pinning machinery (`spec_loader.py`) covers both halves. A
  **config↔parsing consistency check** guards their coherence.
- **Panel *activation* stays in the backend** (`get_visible_panels`); the option
  **definitions** (labels/defaults/sub-options) move to the JSON — data in JSON,
  conditional logic in code.
- **Always-on config lines stay in the backend** — they are VEP-invocation
  invariants, not options.
- **The missing-expected-field check is greenlit**, sequenced after the config
  migration (it needs the API's per-job enabled-plugin knowledge).
- **Display does not fully close like parsing.** The realistic model is a set of
  generic renderer primitives the JSON composes, plus a small shared interactive
  kit — not a per-plugin override registry.

The rest of this document specifies each.

---

## 3. The merged document

Top level is unchanged in shape from today's `ParsingSpec`
(`app/vep/models/parsing_spec_model.py`); the per-entry model grows a `config`,
an `option`, and a `display` section:

```jsonc
{
  "spec_version": "<computed digest, mechanism unchanged>",
  "genome": { "species_taxonomy_id": "9606", "assembly": "GRCh38" },
  "plugins": [
    {
      "id": "protvar",          // option id == ConfigIniParams param name
      "scope": "transcript",
      "output": "protvar",

      "option":  { /* §5 — form option definition (data only) */ },
      "config":  { /* §4 — how the option becomes a config.ini line */ },
      "csq_fields": [ /* … */ ], "require_any_input": [ /* … */ ],
      "targets":  [ /* existing parse half, unchanged */ ],
      "display": { /* §7 — which generic renderer + config, or an override */ }
    }
  ]
}
```

**Config-set ≠ parse-set (important).** The set of config options and the set of
parseable plugins only partly overlap:

- `hgvs`, `spdi`, `protein` are **config flags** that change VEP's base output —
  no dedicated plugin/parser.
- `loeuf`, `geno2mp`, `enformer`, `maxentscan`, `gnomad_mt`, `tss_distance`,
  `nearest_gene`, `nearest_exon_jb` are **config options with no parsing entry**
  in today's 21-plugin spec.

So an entry must allow `config`-without-parse and parse-without-`config`. That
asymmetry is the normal case, and policing it is exactly the consistency check's
job (§6).

---

## 4. Config half — the emitter grammar

Mirror the parsing side's philosophy: a **small closed set of named emitters**
with typed params, no open escape hatch. Four emitters cover everything currently
in `ConfigIniParams.create_config_ini_file`.

### 4.1 `flag`
`hgvs`, `hgvsg`, `spdi`, `protein`:
```jsonc
"config": { "emit": "flag", "keyword": "hgvs" }   // -> "hgvs 1" / "hgvs 0"
```

### 4.2 `plugin`
Static, assembly-keyed files, and sub-flag interpolation all expressed as typed
param values:
```jsonc
// static
"config": { "emit": "plugin", "name": "MaveDB",
            "params": { "file": "{path}/MaveDB_variants.tsv.gz" } }

// assembly-keyed data file (alphamissense, cadd, revel, spliceai, …)
"config": { "emit": "plugin", "name": "AlphaMissense",
            "params": { "file": { "by_assembly": {
              "GRCh38": "{path}/AlphaMissense_hg38.tsv.gz",
              "GRCh37": "{path}/AlphaMissense_hg19.tsv.gz" } } } }

// sub-flag interpolation (protvar, dosage_sensitivity, mutfunc)
"config": { "emit": "plugin", "name": "ProtVar",
            "params": { "db": "{path}/ProtVar_data.db",
              "stability": { "from_option": "protvar_stability", "as": "int" },
              "pocket":    { "from_option": "protvar_pocket",    "as": "int" },
              "int":       { "from_option": "protvar_int",       "as": "int" } } }

// select radio -> flags (TSSDistance direction)
"config": { "emit": "plugin", "name": "TSSDistance",
            "params": {
              "upstream":   { "from_option": "tss_distance_direction", "equals": "upstream",   "as": "int" },
              "downstream": { "from_option": "tss_distance_direction", "equals": "downstream", "as": "int" },
              "both":       { "from_option": "tss_distance_direction", "equals": "both",       "as": "int" },
              "gff3": "{gff}" } }
```

A **param value** is one of:
- a literal string,
- `{ "by_assembly": { "GRCh38": …, "GRCh37": … } }`,
- `{ "from_option": <id>, "as": "int" }` — 1/0 from a boolean sub-option,
- `{ "from_option": <id>, "equals": <value>, "as": "int" }` — 1 when a select
  equals the value.

`{path}` (deploy plugin-data root) and `{gff}` (per-genome resolved gff) are
**backend-provided interpolation tokens**. Like `gff`/`fasta`, they are runtime
values, not static data — they stay backend.

**Variadic sub-flags (IntAct).** IntAct appends one `,<keyword>=1` per selected
sub-option, or `,all=1` when all are selected. Modeled with a `flags` sub-spec:
```jsonc
"config": { "emit": "plugin", "name": "IntAct",
  "params": { "mutation_file": "{path}/mutations.tsv",
              "mapping_file":  "{path}/mutation_gc_map.txt.gz" },
  "flags": { "all_shortcut": "all",
             "options": [
               { "option": "intact_feature_ac",              "keyword": "feature_ac" },
               { "option": "intact_feature_short_label",      "keyword": "feature_short_label" },
               { "option": "intact_feature_annotation",       "keyword": "feature_annotation" },
               { "option": "intact_ap_ac",                    "keyword": "ap_ac" },
               { "option": "intact_interaction_participants",  "keyword": "interaction_participants" },
               { "option": "intact_pmid",                     "keyword": "pmid" } ] } }
```

### 4.3 `custom`
The gnomAD / All-of-Us / ClinVar `custom …` lines.

```jsonc
// ClinVar — static fields, fully declarative
"config": { "emit": "custom",
  "params": { "file": { "by_assembly": {
                 "GRCh38": "{path}/clinvar_GRCh38.vcf.gz",
                 "GRCh37": "{path}/clinvar_GRCh37.vcf.gz" } },
              "short_name": "ClinVar", "format": "vcf", "type": "exact" },
  "fields": { "literal": ["CLNSIG", "CLNSIGCONF"] } }

// gnomAD exomes — combinatorial fields via a NAMED, closed-set builder
"config": { "emit": "custom", "omit_if_no_fields": true,
  "params": { "file": "{path}/gnomad.exomes.v4.1.sites.chr###CHR###.vcf.bgz",
              "short_name": "gnomAD_exomes", "format": "vcf" },
  "fields": { "builder": "gnomad_ancestry_sex", "base": "AF",
              "include_ukb_option": "gnomad_exomes_include_ukb", "join": "%" } }
```

**Named field-builders over open data (the one deliberate non-flat call).**
The gnomAD `fields=` grammar ("for each selected ancestry, for each selected
sex-of-that-ancestry, emit `AF[_non_ukb][_<code>][_<XX|XY>]`") stays a named
builder in a **small closed registry** — `gnomad_ancestry_sex`,
`allofus_populations`, `literal`. This is the exact stance the parsing side took
(a closed transform set, no escape hatch). Crucially the builder is a **closed
algorithm over *open data***: the ancestry/sex/population **labels and field
codes live in the JSON option definitions** (§5), not in the builder. The builder
only combines them.

So "where do the labels/codes come from" → the same place they come from today
(`form_panels.py`'s `_GNOMAD_*_ANCESTRIES` and `pipeline_model.py`'s
`GNOMAD_*_ANCESTRIES` / `GNOMAD_SEXES` / `ALLOFUS_POPULATIONS` tuples), relocated
onto the option sub-options. Adding an ancestry is then JSON-only; a genuinely
new *shape* of grammar is a new named builder (a reviewed code change). This also
decouples the param id from the field code, killing the implicit "`all` → empty
string" special-case.

### 4.4 Emitter gating
Each emitter may carry a genome gate, the config-side analog of the parser's CSQ
`when`:
```jsonc
"config": { "emit": "plugin", "name": "GO", "when": { "assembly": ["GRCh38"] },
            "params": { "file": "{path}/GO.pm_homo_sapiens_116_GRCh38.gff.gz" } }
```

### 4.5 Always-on base config — stays backend
These are written unconditionally (not options, not per-selection). Keep them in
the backend as one explicit base-config definition, next to the per-genome
`gff`/`fasta` resolution (which is a runtime lookup and cannot be static data):

```
force_overwrite 1
numbers 1
symbol 1
biotype 1
transcript_version 1
canonical 1
database 0            # NEW — not currently emitted anywhere
```
Assembly-conditional, also backend:
```
mane 1                # GRCh38 AND mouse GRCm39 (current code is correct)
assembly GRCh38|GRCh37
```

Rationale: these are VEP-invocation invariants, not annotation options; the
payload difference is negligible, so decide on ownership, and ownership is
backend. It does not cost the data-driven goal — adding a *plugin* still touches
only the JSON.

> Cleanup already in flight: the stale `symbol`/`biotype` form checkboxes
> (`submission_form.py` `FormConfig` + the dead `Checkbox` class; frontend
> `vepFormConfig.ts` types) — the form display was removed but the models were
> not. `symbol 1`/`biotype 1` are unconditional, so those toggles were
> meaningless. Tracked as its own cross-repo task.

---

## 5. Option definitions vs panel activation

Split **data** from **logic**:

- **Option definitions** (id, label, type, default, category, panel, sub-options,
  and the field codes for AF builders) move into the JSON `option` block — this
  is the static data currently inside `form_panels.py`'s dicts.
- **Panel activation** (the `is_human_grch38` gating and species-tier layering)
  **stays in `get_visible_panels`**, which becomes "given the JSON option catalog
  + this genome's attributes, decide which to show."

Worked gnomAD ancestry sub-option, carrying both label and field code so both the
form and the config builder read from one place:
```jsonc
{ "id": "gnomad_exomes_afr", "label": "African & African-American",
  "field_code": "afr", "default": false, "sex_split": true,
  "sub_options": [
    { "id": "gnomad_exomes_afr_both",   "label": "Both",   "field_code": "",   "default": true  },
    { "id": "gnomad_exomes_afr_female", "label": "Female", "field_code": "XX", "default": false },
    { "id": "gnomad_exomes_afr_male",   "label": "Male",   "field_code": "XY", "default": false } ] }
```
(`all` → `field_code: ""`; genomes `grpmax` → `{ "field_code": "grpmax",
"sex_split": false }`; All-of-Us populations carry `"field_codes": [...]`, e.g.
`max` → `["gvs_max_af", "gvs_max_subpop"]`.)

Adding a plugin is then one JSON entry. The backend changes only for a genuinely
new *visibility rule* or a new *named builder/transform* — both rare, both
deliberately reviewed.

**Endpoint 1 is relayed, not called directly by the frontend** (decided in
review). The frontend requests the input form from the backend
(`GET /form_config`, keyed on the selected species/assembly, fired on species
selection), and the backend fetches the option catalog from the new API's
*endpoint 1* and runs activation. Reason: the frontend cannot be relied on to
hold the genome metadata needed to resolve options for every species — GRCh38 is
simple, others get more complex — so activation stays server-side. See §10 for
the full request flow.

---

## 6. The two checks

### 6.1 Config↔parsing consistency (static, load-time)
A `model_validator` on the merged model, run when the document loads:

- **Plugin-presence pairing** — surface a `config` entry that produces parseable
  output but has no `targets`, and (softly) a `targets` entry with no `config`.
- **Column-level, exact for `custom` emitters** — a `custom` line's `short_name`
  + `fields` literally determine the CSQ column names (`gnomAD_exomes` + `AF` →
  `gnomAD_exomes_AF`), which must match the parse half's
  `csq_fields`/`from_pattern`.
- **Presence-level only for `plugin` emitters** — VEP derives a plugin's CSQ
  column names internally, so the config line does not state them; the check can
  only verify the plugin/parser pairing, not the exact columns.

### 6.2 Missing-expected-field (runtime, per-job) — greenlit
Given the job's **enabled** options: expected columns = union of
(custom emitter → derived column names; plugin emitter → its `csq_fields`). Fail
loud if a required column is absent from the VCF header; keep silently ignoring
extras (the tolerance is intentionally one-directional). gnomAD's expected set
comes from the **same** `gnomad_ancestry_sex` builder that wrote its config line —
one reason the builder is shared code, not per-side duplication.

Sequencing: this depends on the config→JSON migration (it needs the enabled-plugin
knowledge), so it ships as part of / after that work, not standalone.

**On a miss (operational, decided in review):** in production the backend reruns
the pipeline, which is responsible for adding the required headers; retries are
**capped at 3**, after which it fails with a logged detail. In dev it only warns
(the manual loop has no rerun). Missing headers should be rare — the headers are
required by the parser — so the cap is a safety valve, not an expected path.

---

## 7. Display half

### 7.1 What is already data-driven
Layout is **already** driven by the `panels` contract: `renderPanel` /
`renderOption` in `VepResultsAnnotationDetail.tsx` walk panel → category → option
generically. What is hardcoded is only the **per-option value node** — the big
`switch (optionId)` in `optionContent`.

### 7.2 Display does not fully close — and that is fine
Unlike parsing (finite structural transforms), presentation is open-ended.
Reading the whole renderer, the per-option nodes split three ways:

- **Category 0/2 — a generic renderer vocabulary** (the majority). Primitives the
  default renderer needs, which the JSON composes per plugin:
  - `labelled-rows` — `{ heading?, rows: [{ path, label, format, when_present }] }`
    with a fixed formatter set (`num`, `score`, `humanize`, `mono`, `join`).
    Covers revel, cadd, loeuf, alphamissense, hgvs, spliceai,
    dosage_sensitivity, utrannotator, riboseqorfs, mutfunc's default view.
  - `repeat` — map an array → a row each (mavedb assays, protvar
    pockets/interfaces, go terms, phenotypes, opentargets associations,
    frequency populations).
  - `show-more` — the truncated "+N more" list, currently reimplemented four
    times (mavedb, go, phenotypes, gwas); one primitive.
  - `link-template` — single-field external URL, `{ template: "…/{id}" }`
    (go AmiGO, opentargets target/disease, mavedb score-sets). Data, not logic.
  - `decode` — code → label via the option catalog. The frequency population
    labels (`frequencyPopulationLabels.ts`) are today a **hand-synced copy** of
    `form_panels.py` labels (per its own comment) — the exact drift this project
    removes. It is the inverse of the `gnomad_ancestry_sex` builder: a closed
    algorithm over the same open label data.
  - Sub-option enumeration (Show-all) is already generic (`renderRunSubOptions` +
    `subOptionRan`), parameterized by the sub-option list that moves to the JSON.

  With these, mavedb / go / phenotypes become fully data-driven and mutfunc is
  near-trivial.

- **Category 3 — judgment calls** (a few small cross-field conditionals). Either
  add conditional/derived operators to the vocabulary, or keep as overrides:
  clinvar (two shapes on `conflicting_breakdown` empty/non-empty), intact
  (derived count + pluralization), all-of-us "max" subpopulation bracket,
  opentargets (GWAS/QTL compound, multi-link rows).

- **Category 1 — genuinely frontend-only**, now a *small shared interactive kit*,
  not a per-plugin registry:
  - the **view-in-app popup** widget (`ViewInAppPopup` + `PointerBox`), shared
    and generic, fed template URLs (§7.3);
  - the **show-more** toggle (also a generic primitive);
  - **temporarily**, ProtVar's `buildProtvarUrl` variant minimisation (§7.3).

The JSON's `display` block is therefore either
`{ "renderer": "default", … }` (composing the primitives above) or
`{ "renderer": "override:opentargets" }` for the shrinking Category-3 tail.

### 7.3 In-app popup links move off the resolver to templates
The standalone frontend **does not have the host-router resolver** the integrated
Ensembl app has, so it cannot use resolved root-relative feature-explorer paths.
All four popups build their links from a **simple template** instead
(link-construction logic to be supplied):

| Popup | Call site | Links (intended) |
|---|---|---|
| Location | `VepResultsLocation.tsx:60` | genome browser |
| Gene | `VepResultsGene.tsx:57` | genome browser + feature browser |
| Transcript | `VepSubmissionResults.tsx:729` | feature explorer **+ add genome browser** |
| Protein | `VepResultsAnnotationDetail.tsx:188` | feature explorer |

Covers **both** link builders — `featureExplorerUrls` and `urlFor.browser`
(genome browser). The `protein` popup should share the `transcript` popup's code.
Once templated, the popup is a shared generic primitive fed data URLs; ProtVar's
algorithmic URL stays temporarily and disappears when upstream ProtVar exposes
template-constructible links.

---

## 8. Version pinning (unchanged, now covers both halves)

`spec_loader.py` already computes `spec_version` as a content digest of the
*validated model's canonical dump* and sidecars the pinned document beside the
job (`parsing_spec.json`). Merging config+parse into one document means one digest
pins the whole front↔back contract for a job; a step1↔step3 mismatch 409s. No
mechanism change — the merge is why Q1 ("join them") is the less brittle choice.

---

## 9. Sequencing & open items

**Workstreams (independent seams, can proceed in parallel):**
1. **Go-flat results** — switch the response to
   `spec_interpreter.apply_plugin_spec`; generic `{plugin, scope, data}` envelope;
   frontend generic types. (The results-time seam already loads the pinned spec.)
2. **Config→JSON migration** — lift `pipeline_model.py`'s generation into the
   emitter grammar + a thin runtime interpreter; the base config stays backend.
   Bundles the consistency check (§6.1) and the missing-field check (§6.2).
3. **Display DSL** — the generic renderer primitives (§7.2) + the popup-template
   work (§7.3) + the small override registry.

**Resolved (were awaiting user input):**
- The **link-construction templates** for the four popups — supplied and
  **shipped**; both link builders now emit templated `beta.ensembl.org` URLs.
- The future **ProtVar template** that would retire `buildProtvarUrl` — still
  pending upstream; `buildProtvarUrl` stays algorithmic for now.

**Parked (lack of data / out of scope for now):**
- SIFT/PolyPhen (relocated to typed consequence fields at the go-flat cutover —
  no longer inside a `_parse_pathogenicity`, which was deleted), uniprot,
  protein matches/DOMAINS — still unspecced for lack of a sample VCF carrying
  the columns.
- AF `kind: "allele_frequency"` marker on `PluginSpec` (still not added).

**Cross-repo layout:** backend owns the merged JSON, its model/interpreter, the
two checks, and the always-on base; frontend owns the renderer primitives +
override registry + popup templates. Use the same branch name in both repos when
implementation starts (this design branch is `feature/config-display-spec-design`
in `ensembl-web-tools-api`).

---

## 10. End-to-end request flow (verified against current code)

The new API exposes **two endpoints**:
- **Endpoint 1** — the available options + how to render them on the input form,
  keyed on the selected species/assembly.
- **Endpoint 2** — the config lines, parsing spec, and results-display spec for a
  specific set of *selected* options.

One submission, start to finish (order confirmed in both repos):

1. **Input form.** On species selection the frontend calls the backend
   `GET /form_config` (`vepApiSlice.ts`; fired via `skip: !selectedSpecies` in
   `VepFormOptionsSection.tsx`). The backend fetches endpoint 1, applies
   activation (`get_visible_panels`), and returns panels + options + render info.
   The call is **relayed, not direct** (§5). The frontend renders the form with
   generic renderers.
2. **Submit.** The user picks options, adds a VCF, and the frontend `POST`s
   `/submissions` (`vep_resources.py::submit_vep`). The backend sends the selected
   options to endpoint 2 and receives config + parsing + display.
3. **Check + pin.** The backend runs the consistency check over the returned
   document (config ⇄ parsing ⇄ display); on failure it retries the fetch, max 3,
   then fails and logs. It merges the config with the always-on base (§4.5), emits
   `config.ini`, and pins the parsing + display spec beside the job
   (`spec_loader.write_spec_sidecar`).
4. **Run.** Prod: launch a Nextflow run via Seqera (`nextflow.launch_workflow`),
   output mounted for the backend. Dev: dump `config.ini` to `dev-data`; the
   output VCF is hand-placed there (the manual loop).
5. **Poll.** The frontend polls `GET /status` every 15s
   (`vepSubmissionStatusPolling.ts`, `POLLING_INTERVAL`). Prod polls Seqera
   (`get_workflow_status`); dev reports SUCCEEDED immediately.
6. **Results.** The frontend calls `GET /results` (`get_results_from_path` +
   `_load_pinned_spec`). The backend runs the missing-header check (§6.2), parses
   with the **pinned** parsing spec, and returns the annotations plus the
   **pinned** display spec. The frontend renders with generic renderers (+ the
   small custom kit). NB today the results view re-fetches *live* form-config
   panels for layout — pinning the display spec fixes that.
7. **Later.** Filtering and download stay server-side (`parse_filters`,
   `stream_vep_tsv`).

A rendered version of this flow (dev/prod branches inline) lives beside this doc
as `dataflow-diagram.html` — a self-contained, theme-aware SVG (no runtime). It is
generated by `dataflow-diagram.py` (edit the `EVENTS` list and re-run); keep this
section and the diagram in step.
