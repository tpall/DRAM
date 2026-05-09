"""Unit tests for combine_annotations.py --dbcan_dir parsing.

Exercises the full click entry point against fabricated run_dbcan output
to verify that the dbcan3 wiring produces the expected raw-annotations.tsv
schema without needing nextflow, containers, or a real DBCAN DB.

Run with:  pytest tests/unit/test_combine_annotations_dbcan.py
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
COMBINE_SCRIPT = ROOT / "bin" / "combine_annotations.py"

GENES_FAA = """\
>sampleA_scaffold1_1 # 1 # 300 # 1 # ID=1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIE
>sampleA_scaffold1_2 # 350 # 700 # -1 # ID=2
MELVKSTPDPLAQRQQLAGFETPYAARTHHL
>sampleA_scaffold2_1 # 1 # 250 # 1 # ID=3
MQLNAAEKLAEPHLPCRYALGAESAGPLLQR
"""

KOFAM_CSV = """\
query_id,kofam_id,kofam_bitScore,kofam_description
sampleA_scaffold1_1,K00001,180.5,alcohol dehydrogenase
sampleA_scaffold2_1,K01902,250.0,succinyl-CoA synthetase alpha subunit
"""

DBCAN_HMM_TSV = (
    "Target Name\tHMM Name\ti-Evalue\tCoverage\n"
    "sampleA_scaffold1_1\tGH18.hmm\t1.2e-30\t0.95\n"
    "sampleA_scaffold1_2\tCBM50.hmm\t3.4e-15\t0.88\n"
)

DBCAN_SUB_TSV = (
    "Target Name\tSubfam Name\tSubfam Composition\tSubfam EC\tSubstrate\ti-Evalue\n"
    "sampleA_scaffold1_1\tGH18_e1\tGH18\t3.2.1.14\tchitin\t5.6e-25\n"
    "sampleA_scaffold2_1\tGH5_e3\tGH5|GH5_4\t3.2.1.4\tcellulose\t1.1e-18\n"
)


@pytest.fixture
def fixture_dirs(tmp_path):
    """Build annotations/, genes/, dbcan/ directories with fabricated content
    and return a (annotations_dir, genes_dir, dbcan_dir, output_tsv) tuple."""
    annotations = tmp_path / "annotations"
    genes = tmp_path / "genes"
    dbcan = tmp_path / "dbcan"
    annotations.mkdir()
    genes.mkdir()
    dbcan.mkdir()

    (genes / "sampleA_called_genes.faa").write_text(GENES_FAA)
    (annotations / "sampleA___kofam_formatted.csv").write_text(KOFAM_CSV)
    (dbcan / "sampleA_dbCAN_hmm_results.tsv").write_text(DBCAN_HMM_TSV)
    (dbcan / "sampleA_dbCANsub_hmm_results.tsv").write_text(DBCAN_SUB_TSV)

    return annotations, genes, dbcan, tmp_path / "raw-annotations.tsv"


def _run_combine(annotations, genes, dbcan, output):
    env = os.environ.copy()
    env["FASTA_COLUMN"] = "input_fasta"
    env["PYTHONPATH"] = str(ROOT / "bin")
    subprocess.run(
        [
            sys.executable,
            str(COMBINE_SCRIPT),
            "--annotations_dir", str(annotations),
            "--genes_dir", str(genes),
            "--dbcan_dir", str(dbcan),
            "--output", str(output),
        ],
        env=env,
        check=True,
    )


def test_dbcan_columns_present_and_no_bitscore(fixture_dirs):
    annotations, genes, dbcan, output = fixture_dirs
    _run_combine(annotations, genes, dbcan, output)
    df = pd.read_csv(output, sep="\t")

    expected_dbcan = {
        "dbcan_id", "dbcan_i_Evalue",
        "dbcan_sub_id", "dbcan_sub_composition", "dbcan_sub_ec",
        "dbcan_sub_substrate", "dbcan_sub_i_Evalue",
    }
    missing = expected_dbcan - set(df.columns)
    assert not missing, f"missing dbcan columns: {missing}"
    assert "dbcan_bitScore" not in df.columns, \
        "legacy dbcan_bitScore column should be gone with run_dbcan path"


def test_hmm_suffix_stripped(fixture_dirs):
    annotations, genes, dbcan, output = fixture_dirs
    _run_combine(annotations, genes, dbcan, output)
    df = pd.read_csv(output, sep="\t")

    by_query = df.set_index("query_id")["dbcan_id"].to_dict()
    assert by_query["sampleA_scaffold1_1"] == "GH18", \
        "GH18.hmm should be stripped to GH18"
    assert by_query["sampleA_scaffold1_2"] == "CBM50"


def test_sub_columns_join_independently_of_hmm(fixture_dirs):
    """sampleA_scaffold2_1 has a sub_hmm hit but no main hmm hit; both rows
    should still appear and the sub_* columns should be populated for it."""
    annotations, genes, dbcan, output = fixture_dirs
    _run_combine(annotations, genes, dbcan, output)
    df = pd.read_csv(output, sep="\t")

    s2 = df[df["query_id"] == "sampleA_scaffold2_1"].iloc[0]
    assert pd.isna(s2["dbcan_id"]) or s2["dbcan_id"] == ""
    assert s2["dbcan_sub_id"] == "GH5_e3"
    assert s2["dbcan_sub_substrate"] == "cellulose"


def test_input_fasta_tagging_from_dbcan_filename(fixture_dirs):
    """The _dbCAN splitter must derive 'sampleA' from
    'sampleA_dbCAN_hmm_results.tsv', matching the kofam-side 'sampleA'
    derived via the ___ splitter."""
    annotations, genes, dbcan, output = fixture_dirs
    _run_combine(annotations, genes, dbcan, output)
    df = pd.read_csv(output, sep="\t")

    assert (df["input_fasta"] == "sampleA").all(), \
        f"unexpected input_fasta values: {df['input_fasta'].unique()}"


def test_genes_without_any_hit_still_appear(fixture_dirs):
    """Sanity check: every gene from the FAA file appears in the output even
    when no DB matched it. Pre-existing behavior we don't want to regress."""
    annotations, genes, dbcan, output = fixture_dirs
    _run_combine(annotations, genes, dbcan, output)
    df = pd.read_csv(output, sep="\t")

    expected_queries = {
        "sampleA_scaffold1_1", "sampleA_scaffold1_2", "sampleA_scaffold2_1",
    }
    assert set(df["query_id"]) >= expected_queries
