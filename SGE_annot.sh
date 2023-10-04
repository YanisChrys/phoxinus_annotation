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

snakemake  \
    --snakefile workflow/Snakefile_HAP1 \
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

# snakemake -s workflow/Snakefile --dag --until braker3busco \
# | dot -Tpdf > graph_of_jobs.pdf

#    --stats genofish_HAP1-stats_abinitio.txt
#    --conda-frontend conda \
#    --resources mem_mb=500000 \
#    --conda-create-envs-only
#    --rerun-incomplete \
#    --keep-incomplete \


# export PATH="/home/ychrysostomakis/.conda/envs/maker/pkgs/evidencemodeler-2.1.0-hdbdd923_1/opt/evidencemodeler-2.1.0/EVidenceModeler:$PATH"
# export PATH="/home/ychrysostomakis/.conda/envs/maker/pkgs/evidencemodeler-2.1.0-hdbdd923_1/opt/evidencemodeler-2.1.0/EvmUtils/evidence_modeler.pl:$PATH"
# export PATH="/home/ychrysostomakis/.conda/envs/maker/pkgs/evidencemodeler-2.1.0-hdbdd923_1/bin/EVidenceModeler:$PATH"

# export PATH="/share/pool/genofish_and_chips/testannotation/.snakemake/conda/f723a857c81a13a829d027679e8a8bfb_/opt/evidencemodeler-2.1.0/EvmUtils/evidence_modeler.pl:$PATH"
# export PATH="/share/pool/genofish_and_chips/testannotation/.snakemake/conda/f723a857c81a13a829d027679e8a8bfb_/bin/EVidenceModeler
