#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q large.q,medium.q
#$ -N epi_anot_pri_no_tpsi
#$ -pe smp 41
#$ -e /path/to/logfiles
#$ -o /path/to/logfiles
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m baes

module load anaconda3/2022.05
module load singularity/3.10.5_fix
conda activate genofish # should be snakemake 7.24

mkdir -p $PWD/temp
export TMPDIR="$PWD/temp"
export THREADS=$(expr ${NSLOTS} - 1)

snakemake \
    --snakefile workflow/Snakefile.smk \
    --keep-going \
    --latency-wait 300 \
    --use-envmodules \
    --use-conda \
    --use-singularity \
    -j ${THREADS} \
    --verbose \
    --printshellcmds \
    --nolock \
    --rerun-triggers mtime \
    --singularity-args "--home $PWD" \
    --singularity-args "--bind $TMPDIR:$TMPDIR" \
    --default-resources "tmpdir='$TMPDIR'" 

# Revio:
# bam file with HiFi reads (should be fastq)

# Sequel 2:
# bam file with subreads

# 
