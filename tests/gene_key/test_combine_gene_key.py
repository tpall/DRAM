"""Run with python3 -m unittest discover -s tests/gene_key -v (needs pandas, scikit-bio
and click, as in the annotate module environment).

query_id is not unique across genomes: per-sample assemblers reuse contig names
(MEGAHIT k141_*), so two bins in one run can both contain a gene "k141_100_1".
combine_annotations.py used to build its gene table keyed on query_id alone, so the
genome read last overwrote the other's entry and the other genome's rows came out with
no scaffold or coordinates. These tests give two genomes a colliding gene id and check
that each keeps its own coordinates and hits, and that a genome's rows do not change
when the other genome is added to the run.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
SHARED_ID = "k141_100_1"
GENOMES = {
    "binA": {"genes": {SHARED_ID: (1, 300, 1)}, "kegg": {SHARED_ID: "K01206"}},
    "binB": {"genes": {SHARED_ID: (10, 450, -1), "k141_200_1": (5, 200, 1)},
             "kegg": {SHARED_ID: "K01186"}},
}


def combine(genomes):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        ann, genes = tmp / "annotations", tmp / "genes"
        ann.mkdir(); genes.mkdir()
        for g in genomes:
            d = GENOMES[g]
            (genes / f"{g}_called_genes.faa").write_text("".join(
                f">{gid} # {a} # {b} # {s} # ID=1\nMCAACHKL\n" for gid, (a, b, s) in d["genes"].items()))
            (ann / f"{g}___kofam_annotations.csv").write_text(
                "query_id,kofam_id,kofam_bitScore\n" + "".join(f"{gid},{ko},400\n" for gid, ko in d["kegg"].items()))
        out = tmp / "raw-annotations.tsv"
        subprocess.run([sys.executable, str(ROOT / "bin/combine_annotations.py"),
                        "--annotations_dir", str(ann), "--genes_dir", str(genes),
                        "--output", str(out), "--threads", "1"],
                       check=True, cwd=ROOT / "bin", capture_output=True)
        return pd.read_csv(out, sep="\t")


class TestGeneKey(unittest.TestCase):
    def test_colliding_ids_keep_their_coordinates(self):
        df = combine(["binA", "binB"]).set_index(["input_fasta", "query_id"])
        self.assertEqual(len(df), 3)
        self.assertFalse(df["scaffold"].isna().any())
        self.assertEqual(df.loc[("binA", SHARED_ID), "stop_position"], 300)
        self.assertEqual(df.loc[("binB", SHARED_ID), "stop_position"], 450)
        self.assertEqual(df.loc[("binA", SHARED_ID), "kofam_id"], "K01206")
        self.assertEqual(df.loc[("binB", SHARED_ID), "kofam_id"], "K01186")

    def test_genome_rows_independent_of_run(self):
        alone = combine(["binA"]).set_index(["input_fasta", "query_id"])
        both = combine(["binA", "binB"]).set_index(["input_fasta", "query_id"])
        pd.testing.assert_frame_equal(alone.loc[["binA"]].dropna(axis=1, how="all"),
                                      both.loc[["binA"]].dropna(axis=1, how="all"),
                                      check_like=True, check_dtype=False)


if __name__ == "__main__":
    unittest.main()
