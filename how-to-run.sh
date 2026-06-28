#!/bin/bash
for f in final_bins/*.gz; do
  FA=$(basename ${f%.gz})
  gunzip -c ${f} > input_fasta/${FA}
done

# Auto-resume loop: Nextflow's array-jobs + container combo can hit a
# ConcurrentModificationException in BashWrapperBuilder. That's a session-level
# framework error with no task exit code, so errorStrategy/maxRetries can't
# catch it -- it kills the whole run. Re-running with -resume picks up cached
# work and usually clears the (timing-dependent) race on a later attempt.
# Bounded so a genuine, repeatable failure doesn't loop forever.
attempt=1
max_attempts=5
until nextflow run tpall/DRAM -r eluring-prod-dev \
  -c ~/dram_dbs.config \
  --input_fasta /gpfs/helios/home/taavi74/Projects/newborn/results/binrefine/input_fasta \
  --outdir /gpfs/helios/home/taavi74/Projects/newborn/results/DRAM/call-annotate-distill \
  --threads 8 --call --annotate --summarize --qc \
  --use_kofam --use_dbcan --use_merops --use_viral --use_methyl --use_sulfur \
  -profile singularity --slurm --partition main \
  -with-report -with-trace -with-timeline \
  --array_size 20 \
  --queue_size 20 \
  -resume; do
  if [ "${attempt}" -ge "${max_attempts}" ]; then
    echo "ERROR: nextflow failed ${max_attempts} times; giving up." >&2
    exit 1
  fi
  echo "WARN: nextflow exited non-zero (attempt ${attempt}/${max_attempts}); resuming in 30s..." >&2
  attempt=$((attempt + 1))
  sleep 30
done
