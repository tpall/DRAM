# VFDB / pathogenesis annotation

Dev note for the VFDB integration. Place at `docs/dev/vfdb-pathogenesis.md`.

## Status & decision (2026-07-24)

**Shipped as opt-in raw annotation only — no distillate.** Touch-point 1 (the mmseqs
search) is implemented and smoke-validated on branch `feat/vfdb-pathogenesis`:
`--use_vfdb` / `--anno_dbs vfdb` adds `vfdb_id` / `vfdb_bitScore` / `vfdb_description`
to `raw-annotations.tsv`. Touch-points 2 (distill) and 3 (`--add_db`) are **NOT built**,
and the rule_parser `ID_EXPR_DICT` `vfdb_id` one-liner is **not committed** — it is only
needed if the distillate is ever built (rule_parser is the upstream
`WrightonLabCSU/Rule-Parser` submodule, so that would require a fork + re-pointing).

**Why no distillate:** a stricter-threshold pass on real data (3 extraves gut MAGs vs
VFDB setA, 4759 proteins) showed VFDB-presence is mostly *remote homology*, not
virulence. Of 846 hits, median amino-acid identity was **32%**; at a real presence bar
(bit ≥ 120 & ≥ 80% identity) only **10 (~1%)** survived, and those were
enzyme/metabolic/immune-modulation — **zero exotoxins, adherence or invasion factors**.
The dominant VF classes were Nutritional/Metabolic and Immune modulation. So on
commensal gut cohorts a "pathogenesis" distillate would over-report virulence.

**If revived** (e.g. for pathogen-suspect isolates): gate on **%identity + coverage** (a
code change — the mmseqs module filters bit score only), keep it opt-in, label the column
"VF homolog present", and re-add the `ID_EXPR_DICT` line via a Rule-Parser fork. The rest
of this note is the original three-touch-point design, kept for that scenario.

## Scope

Add **virulence-factor gene presence** to DRAM2 via VFDB, as an mmseqs2 search
that flows through the normal annotate → distill path. Deliberately
**gene-presence only** — no resistance point mutations (those would need
AMRFinderPlus `--organism`; a homology search can't call them). Because it's
presence-only, everything stays inside DRAM2's own mmseqs sweep; no external
tool (ABRicate / AMRFinderPlus) is invoked.

DRAM2's native DB set is metabolic (KEGG, UniRef, PFAM, dbCAN, MEROPS, VOGDB,
viral), so pathogenesis is a genuine gap, not a re-enable.

## The one invariant that governs everything

The VFG id must be **identical** in three places:

1. the mmseqs **target DB** headers (first whitespace token of each FASTA record),
2. `vfdb_descriptions.tsv` column 1 (`params.vfdb_list`), and
3. `distill_pathogenesis.tsv` `gene_id`.

If these drift, hits annotate but the distillate silently comes back empty.
The safest guarantee is to generate the description TSV and the distill sheet
from the **same** VFDB FASTA in one pass (`build_vfdb_sheets.R`), and build the
mmseqs DB from that same FASTA:

```
mmseqs createdb VFDB_setA_pro.fas databases/vfdb/mmseqs/vfdb.mmsdb
```

## Touch-point 1 — annotate (the search)

Modelled on the MEROPS/CAMPER mmseqs half (list-file descriptions, no SQL
sidecar). `use_vfdb` is threaded through the same five layers every other DB
uses:

- `nextflow.config` — `use_vfdb`, `vfdb_db`, `vfdb_list`.
- `workflows/dram.nf` — `params.use_vfdb`, `getDBFlag(..., 'vfdb', ...)`,
  and the `ANNOTATE(...)` call.
- `subworkflows/local/annotate.nf` — `take:` + the `DB_SEARCH(...)` call.
- `subworkflows/local/db_search.nf` — the `MMSEQS_SEARCH_VFDB` include, the
  `if (use_vfdb)` search block (mixed into `formattedOutputchannels`), and the
  `DB_channel_SETUP` empties / `if` / emit.

Toggle with `--use_vfdb` or `--anno_dbs vfdb`. Not wired into `--pool_searches`,
so it runs per-genome (like CAMPER); pooling it later is the same prefix-split
pattern in `search-bundling.md`.

### Confidence tiers, for free

`mmseqs_add_descriptions.py` reads optional `A_rank` / `B_rank` bit-score
columns from the description TSV and stamps a `vfdb_rank` (A/B/none) per hit.
That is the AMRFinderPlus-style per-gene cutoff idea with **no extra code** —
put curated cutoffs in `vfdb_descriptions.tsv`. Use a stricter score/e-value
for virulence than for metabolism: VF domains are shared with commensal
housekeeping genes, so loose thresholds inflate false positives.

## Touch-point 2 — distill (surfacing it)

- New sheet `bin/assets/forms/distill_sheets/distill_pathogenesis.tsv`
  (schema per `distill_metals.tsv`), grouping VFG ids into pathogenesis
  categories (adhesion, T3SS/T4SS/T6SS, toxins, iron acquisition, immune
  evasion, biofilm).
- One clause in `bin/distill.py`, mirroring the CAMPER auto-include:

  ```python
  if "vfdb_id" in annotations and ("default" in distil_topics or "pathogenesis" in distil_topics):
      distil_sheets_names.append(DISTILL_DIR / "distill_pathogenesis.tsv")
  ```

**No Rule-Parser edit needed.** `prepare_present_map_df` auto-discovers any
`*_id` column and routes `vfdb_id` (absent from `ID_EXPR_DICT`) through the
`DEFAULT` handler — a bare single best-hit id, which is exactly right for
top-hit gene presence. (Contrast `vogdb_id`, which had to be registered only
because it needs regex extraction.)

Run with `--distill_topic pathogenesis`. Drop the `"default" in distil_topics or`
clause if you want it strictly opt-in rather than appearing on default runs.

## Touch-point 3 — retrospective add (`--add_db`)

Adds VFDB to a run that already finished **and whose Nextflow work dir was
deleted** (so `-resume` is impossible). It reads only *published* outputs under
`<outdir>/ANNOTATE/`, so the work dir is irrelevant.

`subworkflows/local/add_db.nf`: published `*_called_genes.faa` → `GENE_LOCS`
(regenerate locs) → `MMSEQS_INDEX` → `MMSEQS_SEARCH` (VFDB) →
`JOIN_DB_ANNOTATIONS` → `SUMMARIZE`. Nothing expensive is recomputed.

```
nextflow run . --add_db \
  --input_genes  <outdir>/ANNOTATE/PRODIGAL \
  --annotations  <outdir>/.../raw-annotations.tsv \
  --distill_topic pathogenesis --outdir add_db_out
```

Two design decisions, both forced by the existing code:

- **Join keys on `[query_id, input_fasta]`, not the 5-key `ADD_ANNOTATIONS`
  merge.** `generate_faa_gene_loc_tsv.py` emits no strandedness, so the
  proteins-only path can't satisfy a key that includes it — a 5-key outer merge
  would double rows instead of adding columns. And megahit `k141_*` ids collide
  across bins, so `query_id` alone isn't unique; the genome tag is required.
- **Only `query_id` + `vfdb_*` columns cross the join**, so it strictly adds
  columns and cannot perturb existing annotations.

The genome tag is recovered from the mmseqs output filename and normalised by
stripping the `_called_genes` suffix so it matches `input_fasta` in the old
table.

## Smoke tests

- **Annotate:** diff the sorted `raw-annotations.tsv` against a non-VFDB run —
  the `vfdb_*` columns should be the only delta.
- **Distill:** `metabolism_summary.xlsx` / Product gains a `pathogenesis`
  group; everything else byte-identical (gated on the `vfdb_id` column).
- **add_db:** `raw-annotations-added.tsv` has the **same row count** as the
  original `raw-annotations.tsv` (left join, no inflation) plus `vfdb_*`
  columns. Row growth ⇒ genome-tag normalisation mismatch — check the
  `_called_genes` strip yields tags matching `input_fasta`.

## Gotchas

- Verify the FASTA-header regexes in `build_vfdb_sheets.R` against your VFDB
  release; setA/setB header formats drift.
- `--add_db` assumes `--input_genes` proteins came from the **same** Prodigal
  run as `--annotations`; a different calling silently under-matches on
  `query_id`.
- setA = experimentally verified (default); setB = full/predicted (higher
  recall, more noise — expose as opt-in).
