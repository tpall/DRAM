#!/usr/bin/env bash
#SBATCH --job-name=dram-dbcan3
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=120:00:00
#SBATCH --output=logs/dbcan3-%x-%j.out
#SBATCH --error=logs/dbcan3-%x-%j.err
#
# eluring metabolic profiling — fresh dbcan3 (run_dbcan v3) run per cohort,
# pinned to the eluring-prod-dev baseline. One sbatch job per cohort.
#
# Submit (from a dir containing this script + a logs/ dir):
#
#   for c in extraves newborn2 samples_65_74_new chulenbayeva curated_metagenomic \
#            sajased-aftekas samples_12_18_new wu_rampelli dengl2025; do
#     sbatch --job-name="dram-$c" --export=ALL,COHORT="$c" eluring-dbcan3.sbatch
#   done
#
# Validate on the smallest cohort (extraves, 70 bins) FIRST, confirm
# ANNOTATE/raw-annotations.tsv has the dbcan_sub_substrate column and
# SUMMARIZE/ is produced, then fire the rest.
#
# Requires ~/dram_dbs.config with the dbcan-v3 paths set (uncommented):
#   dbcan_db           = "${params.db_root}/dbcan-v3/"
#   dbcan_version_file = "${params.db_root}/dbcan-v3/version.txt"
#   dbcan_version      = "5-2_9-13-2025"
set -euo pipefail
: "${COHORT:?set COHORT=<cohort>}"

ROOT=/gpfs/helios/home/taavi74/Projects
INPUT_FASTA="$ROOT/$COHORT/results/binrefine/final_bins"   # native .fa.gz; pipeline decompresses
OUTDIR="$ROOT/$COHORT/results/DRAM/call-annotate-distill"

module purge
module load singularity/3.8.5
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate java

[[ -d "$INPUT_FASTA" ]] || { echo "ERROR: no bins at $INPUT_FASTA" >&2; exit 1; }
mkdir -p "$OUTDIR" logs

echo "==> $COHORT  in=$INPUT_FASTA  out=$OUTDIR  job=${SLURM_JOB_ID:-?}"

nextflow run tpall/DRAM -r eluring-prod-dev \
  -c ~/dram_dbs.config \
  -profile singularity --slurm --partition main \
  --input_fasta "$INPUT_FASTA" --outdir "$OUTDIR" \
  --threads 8 --call --annotate --summarize --qc \
  --use_kofam --use_dbcan --use_merops --use_viral --use_methyl --use_sulfur \
  --array_size 20 --queue_size 20 \
  -with-report   "$OUTDIR/report-${SLURM_JOB_ID}.html" \
  -with-trace    "$OUTDIR/trace-${SLURM_JOB_ID}.txt" \
  -with-timeline "$OUTDIR/timeline-${SLURM_JOB_ID}.html"
