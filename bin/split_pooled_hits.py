#!/usr/bin/env python
"""Split a pooled mmseqs formatted-hits CSV back into per-genome CSVs.

Phase-2 search-bundling (docs/dev/search-bundling.md): when `pool_searches` is set,
PREFIX_GENES_FOR_POOL rewrites each genome's gene ids to `<genome>___<orig_id>` so
they are globally unique (Prodigal ids are only unique within a genome — megahit
k141_* scaffold names collide across bins). All genomes are then searched in one
task, so the formatted-hits CSV contains every genome's hits, each `query_id`
carrying its `<genome>___` prefix.

Downstream (SQL descriptions + COMBINE_ANNOTATIONS) consumes ONE CSV per genome,
attributing rows to the genome in the `___`-delimited filename. So we recover the
genome from each row's prefix, strip the prefix to restore the original `query_id`,
and re-emit `<genome>___mmseqs_<db>_formatted.csv` — identical to what the
per-genome path produces (rows partitioned, values unchanged; final row order is
irrelevant because COMBINE_ANNOTATIONS sorts).

Usage: split_pooled_hits.py <pooled_formatted.csv> <db_name>
"""
import sys
import pandas as pd

DELIM = "___"


def main(pooled_csv, db_name):
    hits = pd.read_csv(pooled_csv)
    if "query_id" not in hits.columns:
        hits = hits.rename(columns={hits.columns[0]: "query_id"})

    # query_id == "<genome>___<original_id>"; split on the FIRST delimiter only.
    split = hits["query_id"].astype(str).str.split(DELIM, n=1, expand=True)
    if split.shape[1] < 2 or split[1].isna().any():
        sys.exit(
            f"ERROR: {int(split[1].isna().sum()) if split.shape[1] > 1 else len(hits)} "
            f"query_ids lack a '{DELIM}' genome prefix — was PREFIX_GENES_FOR_POOL run?"
        )
    hits = hits.assign(__genome=split[0], query_id=split[1])

    n = 0
    for genome, sub in hits.groupby("__genome", sort=True):
        sub = sub.drop(columns="__genome")
        sub.to_csv(f"{genome}{DELIM}mmseqs_{db_name}_formatted.csv", index=False)
        n += 1
    sys.stderr.write(f"split pooled {db_name} hits into {n} per-genome CSVs\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
