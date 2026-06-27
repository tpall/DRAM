# Search bundling & workflow-cache efficiency — design note

Status: **Phase-2 implemented + validated for mmseqs (merops/viral/methyl)** on branch
`perf/search-bundling` (off `dev`). Author: investigation 2026-06-27.

## Validated result (extraves, 56 bins, 2026-06-27)
`--pool_searches` for merops + viral + methyl is **byte-identical** to the per-genome
path (136,240 rows; every merops/viral/methyl column matches, order-insensitive). The
per-genome `MMSEQS_INDEX` is fully collapsed (0 tasks under pooling). mmseqs
index+search+split tasks: **224 → 7** for 3 DBs at 56 bins; on a 4621-bin cohort the
per-genome path would be ~18k mmseqs tasks vs ~7. Earlier single-DB run measured
6.1 → 1.5 CPU-h and ~12× less aggregate queue+exec time.

Key bug the validation caught: Prodigal gene ids are `{scaffold}_{gene}` and unique
only WITHIN a genome (megahit `k141_*` scaffolds collide across bins). Fixed by
`PREFIX_GENES_FOR_POOL` (`<genome>___<id>` on faa headers + gene-locs), stripped in
`SPLIT_POOLED_HITS`.

## Remaining
- **SQL descriptions still run per-genome** (SQL_MEROPS/SQL_VIRAL = 56 tasks each) — they
  consume the per-genome split CSVs. Cheap (sqlite lookups, no DB load) but poolable by
  running SQL on the pooled CSV before the split (next easy win).
- **kegg/pfam/camper/canthyd/uniref** mmseqs searches not yet pooled (guarded against
  `--pool_searches`); same ~6-line pattern, not in the eluring DB config so unvalidated here.
- **HMM searches** (kofam/sulfur/etc.) and the single-task COMBINE — Phase-2/3.

## Problem

Annotation DB searches run **per genome**. The genome is the atomic unit end to end:

- `CALL_GENES` is per-fasta — "1-by-1" (`subworkflows/local/call.nf:27`); emits
  `ch_called_proteins = [name, faa]` per genome.
- Every search consumes that per-genome channel (`subworkflows/local/db_search.nf:56`),
  so **N genomes × M databases = N×M search tasks**, each re-staging/loading its DB.
- `MMSEQS_SEARCH` is `process_huge`, `tag { input_fasta }`, one task per genome, symlinking
  the target DB every time (`modules/local/annotate/mmseqs_search.nf:9,30`).
- `MMSEQS_INDEX` is per genome with no `storeDir` (`modules/local/annotate/mmseqs_index.nf`)
  → re-indexes every run.
- `COMBINE_ANNOTATIONS` is a **single serial task** merging all N genomes
  (`subworkflows/local/db_search.nf:398-402`).

Concrete cost: a 4621-bin cohort (curated_metagenomic) launches ~4621 tasks **each** for
mmseqs-merops / viral / methyl / (kegg) / pfam / uniref, kofam-HMM, dbcan, and mmseqs-index —
tens of thousands of `process_huge` tasks. The search compute per small bin is trivial; the
**per-task overhead** (SLURM scheduling, container spin-up, DB staging) × N×M dominates.

### Why this is also an operational problem (observed 2026-06)
- **Array+container clone race / exit-127**: driven by huge SLURM array-submission bursts.
  Fewer, larger tasks shrink the burst.
- **Fragile resume** (wrong-session, `.nextflow/history` truncation): the cache holds
  N_genomes×M entries (tens of thousands) → slow to load, easy to corrupt. Pooling cuts it ~50×.
- **Disk-quota exhaustion**: thousands of per-task work dirs.

## Enabler: result attribution is already genome-safe

Pooling is only viable because hits can be split back to the originating genome losslessly:

- Prodigal gene IDs are `{genome}_{scaffold}_{gene_number}` (`call_genes_prodigal.nf`).
- Output filenames embed the genome via a `___` delimiter
  (`hmmsearch.nf`, `mmseqs_search.nf:24-25`).
- `combine_annotations.py` recovers the genome from both
  (`input_fasta_from_filepath()` split on `___`, lines 26-27, 108) and tags each row with
  `FASTA_COLUMN`, merging on `[query_id, FASTA_COLUMN]`.

So a pooled search's `convertalis`/hmmsearch output already carries genome identity in every
query id — splitting back is a `groupBy(genome-prefix)`.

## Dormant machinery to reuse

A `splitFasta → search → groupTuple → CONCAT_HMM_HITS` chunk path already exists for kofam/vog
(`db_search.nf:160-198, 343-381`) but is (a) **within a single genome**, (b) wired for only 2 of
~13 searches, (c) **off by default** (`kofam_chunk_size=0`, `vog_chunk_size=0`,
`nextflow.config:118-119`). The same pattern, lifted to operate **across genomes**, is the core
of the fix.

## Proposed direction — pool-then-chunk

1. Concatenate all genomes' `.faa` (already `.collect()`ed into `ch_collected_faa` for QUAST,
   `call.nf:35-38`).
2. `splitFasta(by: search_chunk_size, file: true)` into fixed-size chunks (~100k genes).
3. Index + search **each chunk once per DB** (chunk count ≪ genome count for many small MAGs).
4. Split hits back to genomes by the `{genome}_...` query-id prefix and `groupTuple`, then feed
   the existing per-genome `COMBINE_ANNOTATIONS` (or a parallel merge — see Phase 1).

Effect: MEROPS for curated drops from **4621 tasks → ~90 chunks**; DB loaded ~90× instead of
4621×. Cache becomes chunk-granular (editing one genome re-runs only its chunk).

### Exception: dbcan
`RUNDBCAN_EASYSUBSTRATE` (run_dbcan `easy_substrate`) needs each genome's GFF and is inherently
per-genome. It stays per-genome unless run_dbcan exposes a batch mode (Phase 3).

## Phasing

- **Phase 1 (low-risk):**
  - `storeDir` on `MMSEQS_INDEX` (deterministic per protein set) — opt-in via a path param so
    default behaviour is unchanged; caveat: storeDir keys on output *name*, so a renamed-but-
    changed genome could reuse a stale index (fine for stable cohort inputs).
  - Parallelize `COMBINE_ANNOTATIONS` (tree/partial merge instead of one serial N-genome task).
  - Surface the dormant chunk-size knobs.
- **Phase 2 (the win):** pool-then-chunk across genomes for mmseqs + hmm searches, with split-back.
- **Phase 3:** run_dbcan batch-mode investigation.

## Validation gate

Bundling is purely an **execution-layout** change — outputs must be identical. Run **extraves
(56 bins)** old-way vs new-way and `diff` `ANNOTATE/raw-annotations.tsv` (sorted): must match
modulo row order. Same smoke-test pattern already used to gate the dbcan3 rollout.

## Scaffolding added on this branch
- `params.pool_searches` (bool, default `false`) and `params.search_chunk_size` (int, default `0`)
  in `nextflow.config` — Phase-2 entry points, defaulting to current behaviour (no-op until wired).
- This note.
