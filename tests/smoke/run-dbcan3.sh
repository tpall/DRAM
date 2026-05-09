#!/usr/bin/env bash
# Smoke test for the run_dbcan v3 (RUNDBCAN_EASYSUBSTRATE) path.
#
# Runs --call --annotate --use_dbcan on the OWC fixture (tests/data/owc),
# isolating the dbcan path so failures point clearly at the new wiring
# rather than other DBs. Writes to results/dbcan3-smoke/ by default.
#
# Run from the repo root on HPC, on branch dramv-dbcan3 or later.
#
# Env overrides:
#   DBCAN_DB             dir holding the run_dbcan v3 distribution + version.txt
#   DBCAN_VERSION_FILE   path to one-line version file (default: $DBCAN_DB/version.txt)
#   DBCAN_VERSION        expected version string (default: 5-2_9-13-2025)
#   OUTDIR               nextflow output dir (default: results/dbcan3-smoke)
#   PARTITION            slurm partition (default: main)
#   ARRAY_SIZE           --array_size value (HPC fair-share requires > 0; default: 4)

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

: "${DBCAN_DB:?Set DBCAN_DB to the run_dbcan v3 DB directory before running}"
: "${DBCAN_VERSION_FILE:=${DBCAN_DB}/version.txt}"
: "${DBCAN_VERSION:=5-2_9-13-2025}"
: "${OUTDIR:=results/dbcan3-smoke}"
: "${PARTITION:=main}"
: "${ARRAY_SIZE:=4}"

INPUT_FASTA="${REPO_ROOT}/tests/data/owc/input_fasta"

echo "==> Pre-flight"
echo "    branch:             $(git rev-parse --abbrev-ref HEAD)"
echo "    head:               $(git rev-parse --short HEAD)"
echo "    DBCAN_DB:           ${DBCAN_DB}"
echo "    DBCAN_VERSION_FILE: ${DBCAN_VERSION_FILE}"
echo "    DBCAN_VERSION:      ${DBCAN_VERSION}"
echo "    OUTDIR:             ${OUTDIR}"

if [[ ! -d "${DBCAN_DB}" ]]; then
    echo "ERROR: DBCAN_DB does not exist: ${DBCAN_DB}" >&2
    exit 1
fi
if [[ ! -f "${DBCAN_VERSION_FILE}" ]]; then
    echo "ERROR: version file missing: ${DBCAN_VERSION_FILE}" >&2
    echo "       Create it with:  echo -n '${DBCAN_VERSION}' > '${DBCAN_VERSION_FILE}'" >&2
    exit 1
fi
FOUND_VERSION="$(head -n1 "${DBCAN_VERSION_FILE}")"
if [[ "${FOUND_VERSION}" != "${DBCAN_VERSION}" ]]; then
    echo "WARNING: version file contains '${FOUND_VERSION}', expected '${DBCAN_VERSION}'." >&2
    echo "         The pipeline will refuse to start. Fix before rerunning." >&2
fi
if [[ ! -d "${INPUT_FASTA}" ]]; then
    echo "ERROR: test fixture missing: ${INPUT_FASTA}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

echo
echo "==> Launching nextflow"
nextflow run . \
    -profile singularity \
    --slurm --partition "${PARTITION}" \
    --input_fasta "${INPUT_FASTA}" \
    --fasta_fmt '*.fa' \
    --outdir "${OUTDIR}" \
    --call --annotate \
    --use_dbcan \
    --dbcan_db "${DBCAN_DB}" \
    --dbcan_version_file "${DBCAN_VERSION_FILE}" \
    --dbcan_version "${DBCAN_VERSION}" \
    --threads 4 \
    --array_size "${ARRAY_SIZE}" \
    --queue_size "${ARRAY_SIZE}" \
    -with-report "${OUTDIR}/report.html" \
    -with-trace  "${OUTDIR}/trace.txt"

echo
echo "==> Post-run checks"

RAW="${OUTDIR}/ANNOTATE/raw-annotations.tsv"
RUNDBCAN_DIR="${OUTDIR}/ANNOTATE/RUNDBCAN_EASYSUBSTRATE"

fail=0

if [[ ! -f "${RAW}" ]]; then
    echo "FAIL: raw-annotations.tsv not produced at ${RAW}"
    fail=1
fi

if [[ ! -d "${RUNDBCAN_DIR}" ]]; then
    echo "FAIL: ${RUNDBCAN_DIR} not produced — RUNDBCAN_EASYSUBSTRATE did not publish"
    fail=1
else
    for sample in OWC_0000 OWC_0001; do
        for suffix in dbCAN_hmm_results.tsv dbCANsub_hmm_results.tsv overview.tsv; do
            if ! ls "${RUNDBCAN_DIR}/${sample}_${suffix}" >/dev/null 2>&1; then
                echo "FAIL: missing ${sample}_${suffix} in ${RUNDBCAN_DIR}"
                fail=1
            fi
        done
    done
fi

if [[ -f "${RAW}" ]]; then
    HEADER="$(head -1 "${RAW}")"
    DBCAN_COLS="$(echo "${HEADER}" | tr '\t' '\n' | grep -c '^dbcan_' || true)"
    echo "    dbcan_* column count in header: ${DBCAN_COLS} (expected: 7)"

    for col in dbcan_id dbcan_i_Evalue dbcan_sub_id dbcan_sub_composition \
               dbcan_sub_ec dbcan_sub_substrate dbcan_sub_i_Evalue; do
        if ! echo "${HEADER}" | tr '\t' '\n' | grep -qx "${col}"; then
            echo "FAIL: column '${col}' missing from raw-annotations.tsv header"
            fail=1
        fi
    done

    if echo "${HEADER}" | tr '\t' '\n' | grep -qx "dbcan_bitScore"; then
        echo "FAIL: legacy 'dbcan_bitScore' column still present — should be gone with run_dbcan path"
        fail=1
    fi

    DBCAN_HIT_ROWS="$(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="dbcan_id") c=i; next} c && $c!="" {n++} END {print n+0}' "${RAW}")"
    echo "    rows with non-empty dbcan_id: ${DBCAN_HIT_ROWS}"
    if [[ "${DBCAN_HIT_ROWS}" -eq 0 ]]; then
        echo "FAIL: no dbcan_id hits found — joined 0 rows from run_dbcan output"
        fail=1
    fi
fi

if [[ ${fail} -eq 0 ]]; then
    echo
    echo "==> SMOKE TEST PASSED"
else
    echo
    echo "==> SMOKE TEST FAILED (${fail} issue(s) above)"
    exit 1
fi
