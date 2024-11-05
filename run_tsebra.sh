set -xe

export THREADS=<enter thread # here>

# load singularity in whatever way is best for you
module load singularity

# absolute path to your working directory. This will be bound to the container so all inputs should either be discoverable in here or in your home directory
WORKDIR=/path/to/workdir

# script works with absolute paths

# change these without trailing slash
# location of annotation files
anot_dir="/path/to/results_primary"

# location of galba results
galba_dir="/path/to/galba_pri/GALBA"

# output of this script
tsebra_outdir="/path/to/tsebra_pri"
# id of a specific run with specifc cfg configuration file settings
run_ID="01-close_to_pref"

# input 
# masked genome
genome=${anot_dir}/04_braker3_etp_final/genome.fa

#gtfs
aug_gtf=${anot_dir}/04_braker3_etp_final/augustus.hints.gtf
gm_gtf=${anot_dir}/04_braker3_etp_final/GeneMark-ES/genemark.gtf
keep_all_gtf=${anot_dir}/04_braker3_etp_final/GeneMark-ES/training.gtf
galba_gtf=${galba_dir}/galba.gtf
braker_gtf=${anot_dir}/04_braker3_etp_final/braker.gtf

#gffs
hints_gff1=${anot_dir}/04_braker3_etp_final/hintsfile.gff
hints_gff2=${galba_dir}/hintsfile.gff

# out, log and configs
# download cfg files from the braker github and edit them yourself to see the best output for tsebra
run_outdir=${tsebra_outdir}/${run_ID}
braker_cfg=${tsebra_outdir}/custom_braker3.cfg # edit this
# ${tsebra_outdir}/pref1.cfg
out_gtf=${run_outdir}/braker_galba.gtf
busco_out=${run_outdir}/
out_prefix=braker_galba_busco
log=${run_outdir}/errors/tsebra.stderr


mkdir -p ${run_outdir}/errors/
cp ${braker_cfg} ${run_outdir}/ # for documentation

# ,${aug_gtf},${gm_gtf}
# add or remove gtf files or hintsfiles and see which works best
singularity exec --bind $WORKDIR /path/to/containers/braker3.sif /opt/TSEBRA/bin/tsebra.py \
--gtf ${galba_gtf},${braker_gtf},${aug_gtf},${gm_gtf} \
--hintfiles ${hints_gff2},${hints_gff1} \
--filter_single_exon_genes \
--ignore_tx_phase \
--cfg ${braker_cfg} \
--out ${out_gtf} \
--verbose 1 2> ${log}

# if a gtf file is great , prefare to keep its entries
#--keep_gtf ${keep_all_gtf} \

comb_out_dir=${run_outdir}/braker_galba

singularity exec --bind $WORKDIR /path/to/containers/braker3.sif /opt/conda/bin/python3 \
${anot_dir}/04_braker3_etp_final/scripts/getAnnoFastaFromJoingenes.py \
-g ${genome} \
-f ${out_gtf} \
-o ${comb_out_dir}

# busco
# this takes a long time
# replace with compleasm
conda activate busco

busco_lineage="sauropsida_odb10"
asm=${comb_out_dir}.codingseq

mkdir -p $busco_out

busco -f -m genome \
-c ${THREADS} \
-i ${asm} \
-o ${out_prefix} \
--download_path /share/pool/databases/busco_downloads_v5 \
--out_path ${busco_out} \
-l ${busco_lineage} \
--offline

# code to transform gtf to gff?
# cat ${anot_dir}/04_braker3_etp_final/augustus.hints.gtf | \
# /usr/bin/perl -ne 'if(m/\tAUGUSTUS\t/ or m/\tAnnotationFinalizer\t/ or m/\tGUSHR\t/ or m/\tGeneMark.hmm\t/ or m/\tgmst\t/) {print $_;}' | \
# /usr/bin/perl ${anot_dir}/04_braker3_etp_final/scripts/gtf2gff.pl --gff3 --out=${anot_dir}/04_braker3_etp_final/augustus.hints.gff3 >> ${anot_dir}/04_braker3_etp_final/gtf2gff3.log 2>> ${anot_dir}/04_braker3_etp_final/errors/gtf2gff3.err

conda activate agat
agat_sp_statistics.pl -gff ${GFF_File} -g $Hap1 -p -o Hap1_final_Summary.txt
