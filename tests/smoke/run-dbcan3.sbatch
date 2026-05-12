#!/usr/bin/env bash
#SBATCH --job-name=dram-dbcan3
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/dbcan3-smoke-%j.out
#SBATCH --error=logs/dbcan3-smoke-%j.err
#
# Driver job for the dbcan3 smoke test. Nextflow runs here as the head
# process and dispatches per-task slurm jobs via `--slurm` (queue size is
# capped by ARRAY_SIZE inside run-dbcan3.sh).
#
# Submit from the repo root on the dramv-dbcan3 branch:
#
#   mkdir -p logs
#   sbatch tests/smoke/run-dbcan3.sbatch
#
# Override defaults via env (sbatch --export=ALL passes these through):
#
#   DBCAN_DB=/path/to/dbcan-v3 sbatch tests/smoke/run-dbcan3.sbatch
#
# All env knobs documented in run-dbcan3.sh apply (DBCAN_DB,
# DBCAN_VERSION_FILE, DBCAN_VERSION, OUTDIR, PARTITION, ARRAY_SIZE).

set -euo pipefail

: "${DBCAN_DB:=/gpfs/helios/home/taavi74/databases/DRAM/databases/dbcan-v3}"
export DBCAN_DB

module purge
module load singularity/3.8.5
# `conda activate` needs the conda shell hook loaded under non-interactive bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate java

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

echo "==> sbatch job ${SLURM_JOB_ID:-?} on $(hostname)"
echo "    repo:    ${REPO_ROOT}"
echo "    branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "    DBCAN_DB=${DBCAN_DB}"

exec bash tests/smoke/run-dbcan3.sh
