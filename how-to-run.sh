#!/bin/bash
for f in final_bins/*.gz; do
  FA=$(basename ${f%.gz})
  gunzip -c ${f} > input_fasta/${FA}
done

nextflow run tpall/DRAM -r dev \
  --input_fasta /gpfs/helios/home/taavi74/Projects/newborn2/results/binrefine/input_fasta \
  --outdir /gpfs/helios/home/taavi74/Projects/newborn2/results/DRAM/call-annotate-distill \
  --threads 8 --annotate --summarize --qc \
  --use_kofam --use_dbcan --use_merops --use_viral --use_methyl --use_sulfur \
  -profile singularity --slurm --partition main \
  -with-report -with-trace -with-timeline \
  --array_size 20 \
  --queue_size 20

