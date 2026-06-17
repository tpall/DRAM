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
: "${RESUME:=}"   # set RESUME=1 to reuse cached tasks (e.g. after a transient node failure)

ROOT=/gpfs/helios/home/taavi74/Projects
INPUT_FASTA="$ROOT/$COHORT/results/binrefine/final_bins"   # native .fa.gz; pipeline decompresses
OUTDIR="$ROOT/$COHORT/results/DRAM/call-annotate-distill"
# Per-cohort Nextflow launch dir. Cohorts run concurrently; if they share one
# launchDir, nf-core dumpParametersToJSON races on a second-resolution temp file
# (.params_<ts>.json) -- one run copies+deletes it, another's copy then fails
# with "No such file or directory". A unique launchDir also isolates each
# cohort's work/ and .nextflow/ history, making RESUME unambiguous.
RUNDIR="$ROOT/$COHORT/results/DRAM/nf-run"

RESUME_ARG=""
case "${RESUME,,}" in 1|true|yes|on) RESUME_ARG="-resume" ;; esac

module purge
module load singularity/3.8.5
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate java

[[ -d "$INPUT_FASTA" ]] || { echo "ERROR: no bins at $INPUT_FASTA" >&2; exit 1; }
mkdir -p "$OUTDIR" "$RUNDIR" logs   # logs/ is relative to the submission dir; create before cd
cd "$RUNDIR"                        # unique launchDir per cohort -- see RUNDIR note above

echo "==> $COHORT  in=$INPUT_FASTA  out=$OUTDIR  run=$RUNDIR  job=${SLURM_JOB_ID:-?}"

# Auto-resume loop: Nextflow's array-jobs + container combo can throw a
# ConcurrentModificationException in BashWrapperBuilder -- a session-level
# framework error with no task exit code, so errorStrategy/maxRetries can't
# catch it and the whole run dies. Each retry forces -resume; since launchDir
# is per-cohort (RUNDIR), the cache from the failed attempt is reused, so only
# the unfinished DB_SEARCH submission phase re-rolls. Bounded so a genuine,
# repeatable failure doesn't loop forever.
resume_arg="$RESUME_ARG"
attempt=1
max_attempts=5
until nextflow run tpall/DRAM -r eluring-prod-dev \
  -c ~/dram_dbs.config \
  ${resume_arg} \
  -profile singularity --slurm --partition main \
  --input_fasta "$INPUT_FASTA" --outdir "$OUTDIR" \
  --threads 8 --call --annotate --summarize --qc \
  --use_kofam --use_dbcan --use_merops --use_viral --use_methyl --use_sulfur \
  --array_size 20 --queue_size 20 \
  -with-report   "$OUTDIR/report-${SLURM_JOB_ID}.html" \
  -with-trace    "$OUTDIR/trace-${SLURM_JOB_ID}.txt" \
  -with-timeline "$OUTDIR/timeline-${SLURM_JOB_ID}.html"; do
  if [ "$attempt" -ge "$max_attempts" ]; then
    echo "ERROR: nextflow failed $max_attempts times; giving up." >&2
    exit 1
  fi
  echo "WARN: nextflow exited non-zero (attempt $attempt/$max_attempts); resuming in 30s..." >&2
  resume_arg="-resume"   # force resume on retries to reuse this run's cache
  attempt=$((attempt + 1))
  sleep 30
done
