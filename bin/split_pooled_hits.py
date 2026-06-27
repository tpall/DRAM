#!/usr/bin/env python
"""Split a pooled mmseqs formatted-hits CSV back into per-genome CSVs.

Phase-2 search-bundling (docs/dev/search-bundling.md): when `pool_searches` is set,
all genomes' called proteins are concatenated and searched against a DB in one task.
The resulting formatted-hits CSV therefore contains hits for every genome, identified
by the genome-prefixed `query_id` ({genome}_{scaffold}_{gene}).

Downstream (SQL descriptions + COMBINE_ANNOTATIONS) consumes ONE CSV per genome,
attributing rows to the genome encoded in the `___`-delimited filename. So we
partition the pooled CSV by genome and re-emit `<genome>___mmseqs_<db>_formatted.csv`,
exactly as the per-genome search path would have produced (rows are partitioned, not
re-formatted, so columns/values are identical; final row order is irrelevant because
COMBINE_ANNOTATIONS sorts).

Genome membership is taken from the per-genome gene-locs tables
(`<genome>_called_genes_table.tsv`), each of which lists that genome's `query_id`s —
the same source of truth COMBINE uses, so attribution matches exactly.

Usage: split_pooled_hits.py <pooled_formatted.csv> <gene_locs_dir> <db_name>
"""
import sys
import os
import glob
import pandas as pd

LOCS_SUFFIX = "_called_genes_table.tsv"


def main(pooled_csv, locs_dir, db_name):
    hits = pd.read_csv(pooled_csv)
    if "query_id" not in hits.columns:
        hits = hits.rename(columns={hits.columns[0]: "query_id"})

    # query_id -> genome, from each genome's gene-locs table
    q2g = {}
    locs_files = glob.glob(os.path.join(locs_dir, "*" + LOCS_SUFFIX))
    if not locs_files:
        sys.exit(f"ERROR: no *{LOCS_SUFFIX} files in {locs_dir}")
    for f in locs_files:
        genome = os.path.basename(f)[: -len(LOCS_SUFFIX)]
        locs = pd.read_csv(f, sep="\t", usecols=["query_id"])
        for q in locs["query_id"].tolist():
            q2g[q] = genome

    hits["__genome"] = hits["query_id"].map(q2g)
    unmapped = int(hits["__genome"].isna().sum())
    if unmapped:
        sys.stderr.write(
            f"WARNING: {unmapped} pooled hits had no genome in gene-locs; dropping them\n"
        )
        hits = hits.dropna(subset=["__genome"])

    n = 0
    for genome, sub in hits.groupby("__genome", sort=True):
        sub = sub.drop(columns="__genome")
        sub.to_csv(f"{genome}___mmseqs_{db_name}_formatted.csv", index=False)
        n += 1
    sys.stderr.write(f"split pooled {db_name} hits into {n} per-genome CSVs\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
