#!/bin/bash

# replace with own modules
module load singularity 
module load anaconda3
conda activate genome_annotation

THREADS=50

snakemake \
    --snakefile workflow/Snakefile.smk \
    --keep-going \
    --latency-wait 300 \
    --use-envmodules \
    --use-conda \
    --use-singularity \
    -j ${THREADS} \
    --singularity-args "--home $PWD" \
    --singularity-args "--bind $PWD/temp:/tmp" \
    -j ${THREADS} \
    --default-resources "tmpdir='/path/to/tempdir'" \
    --verbose \
    --printshellcmds \
    --reason \
    --nolock \
    --rerun-triggers mtime \
    --stats "./stats.json" 
    
        
        #--report "./report.html"

