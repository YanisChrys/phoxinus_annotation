#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q medium.q,large.q
#$ -N genofish_HAP1
#$ -pe smp 30
#$ -M I.Chrysostomakis@leibniz-lib.de
#$ -m beas

# T.Oriowo@leibniz-lib.de

module load singularity/3.10.5_fix 
module load anaconda3/2022.05
conda activate genofish

#one core will be used by snakemake to monitore the other processes
THREADS=$(expr ${NSLOTS} - 1)

snakemake -n \
    --snakefile workflow/Snakefile.smk \
    --keep-going \
    --latency-wait 300 \
    --use-envmodules \
    --use-conda \
    --use-singularity \
    --singularity-args "--home $PWD" \
    --singularity-args "--bind $PWD/temp:/tmp" \
    -j ${THREADS} \
    --default-resources "tmpdir='/share/pool/genofish_and_chips/testannotation/hap1_temp'" \
    --verbose \
    --printshellcmds \
    --reason \
    --nolock \
    --rerun-triggers mtime \
    --stats genofish_HAP1-stats_braker3maker2.txt

