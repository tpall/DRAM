#!/usr/bin/env bash
#SBATCH --job-name=dram-dbcan3-real
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/dbcan3-real-%j.out
#SBATCH --error=logs/dbcan3-real-%j.err
#
# Full pipeline (--call --annotate --summarize --qc) against real bins to
# validate the dbcan3 wiring with actual CAZyme hits — the OWC smoke fixture
# is too small to populate dbcan_id rows, so this is the positive-path check.
#
# Mirrors the historical how-to-run.sh layout but pins to the local
# dramv-dbcan3 clone (no `tpall/DRAM -r dev` indirection), adds --call (now
# required by the run_dbcan path), and points --dbcan_db at the v3 install.
#
# Submit from the repo root on the dev branch:
#
#   mkdir -p logs
#   sbatch tests/integration/run-dbcan3-realbins.sbatch
#
# Override via env:
#
#   INPUT_FASTA   dir of decompressed bin FASTAs (default: newborn2 binrefine)
#   OUTDIR        nextflow output dir (default: parallel to existing results)
#   DBCAN_DB      run_dbcan v3 DB dir
#   ARRAY_SIZE    --array_size (default 20; HPC fair-share needs > 0)
#   RESUME        1/true/yes/on to pass -resume (reuse cached tasks; default off)

set -euo pipefail

: "${INPUT_FASTA:=/gpfs/helios/home/taavi74/Projects/newborn2/results/binrefine/input_fasta}"
: "${OUTDIR:=/gpfs/helios/home/taavi74/Projects/newborn2/results/DRAM/call-annotate-distill-dbcan3}"
: "${DBCAN_DB:=/gpfs/helios/home/taavi74/databases/DRAM/databases/dbcan-v3}"
: "${DBCAN_VERSION_FILE:=${DBCAN_DB}/version.txt}"
: "${DBCAN_VERSION:=5-2_9-13-2025}"
: "${ARRAY_SIZE:=20}"
: "${PARTITION:=main}"
: "${RESUME:=}"

# Translate the RESUME env toggle into nextflow's -resume flag. Default off, so
# the integration test runs clean unless a resume is explicitly requested.
RESUME_ARG=""
case "${RESUME,,}" in
    1|true|yes|on) RESUME_ARG="-resume" ;;
esac

module purge
module load singularity/3.8.5
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate java

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

echo "==> sbatch job ${SLURM_JOB_ID:-?} on $(hostname)"
echo "    repo:        ${REPO_ROOT}"
echo "    branch:      $(git rev-parse --abbrev-ref HEAD)"
echo "    head:        $(git rev-parse --short HEAD)"
echo "    INPUT_FASTA: ${INPUT_FASTA}"
echo "    OUTDIR:      ${OUTDIR}"
echo "    DBCAN_DB:    ${DBCAN_DB}"
echo "    resume:      ${RESUME_ARG:-(off)}"

if [[ ! -d "${INPUT_FASTA}" ]]; then
    echo "ERROR: INPUT_FASTA does not exist: ${INPUT_FASTA}" >&2
    exit 1
fi
if [[ ! -d "${DBCAN_DB}" ]]; then
    echo "ERROR: DBCAN_DB does not exist: ${DBCAN_DB}" >&2
    exit 1
fi
if [[ ! -f "${DBCAN_VERSION_FILE}" ]]; then
    echo "ERROR: missing ${DBCAN_VERSION_FILE} (expected: ${DBCAN_VERSION})" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

nextflow run . \
    -profile singularity \
    ${RESUME_ARG} \
    --slurm --partition "${PARTITION}" \
    --input_fasta "${INPUT_FASTA}" \
    --outdir "${OUTDIR}" \
    --threads 8 \
    --call --annotate --summarize --qc \
    --use_kofam --use_dbcan --use_merops --use_viral --use_methyl --use_sulfur \
    --dbcan_db "${DBCAN_DB}" \
    --dbcan_version_file "${DBCAN_VERSION_FILE}" \
    --dbcan_version "${DBCAN_VERSION}" \
    --array_size "${ARRAY_SIZE}" \
    --queue_size "${ARRAY_SIZE}" \
    -with-report "${OUTDIR}/report.html" \
    -with-trace  "${OUTDIR}/trace.txt" \
    -with-timeline "${OUTDIR}/timeline.html"

echo
echo "==> Post-run dbcan check"
RAW="${OUTDIR}/ANNOTATE/raw-annotations.tsv"
if [[ ! -f "${RAW}" ]]; then
    echo "ERROR: raw-annotations.tsv not produced at ${RAW}" >&2
    exit 1
fi

HITS="$(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="dbcan_id") c=i; next} c && $c!="" {n++} END {print n+0}' "${RAW}")"
SUB_HITS="$(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="dbcan_sub_id") c=i; next} c && $c!="" {n++} END {print n+0}' "${RAW}")"
echo "    dbcan_id non-empty rows:     ${HITS}"
echo "    dbcan_sub_id non-empty rows: ${SUB_HITS}"
echo "    top 10 dbcan_id values:"
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="dbcan_id") c=i; next} c && $c!="" {print $c}' "${RAW}" \
    | sort | uniq -c | sort -rn | head -10 | sed 's/^/      /'

if [[ "${HITS}" -eq 0 ]]; then
    echo "ERROR: zero dbcan_id hits across real bins — likely a join/parser regression" >&2
    exit 1
fi
echo "==> dbcan3 real-bins run PASSED"
