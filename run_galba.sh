
set -xe
export THREADS=<enter thread # here>

# load singularity in whatever way is best for you
module load singularity

# absolute path to your working directory. This will be bound to the container so all inputs should either be discoverable in here or in your home directory
WORKDIR=/path/to/workdir

# define where results will be saved. Make specific to your run
PREFIX=<run id here>
outdir=${WORKDIR}/annotation/workflow/galba/${PREFIX}
mkdir -p ${outdir}
cd ${outdir}

asm=${WORKDIR}/my_genome.fasta

# use few (~10 or less) well-annotated proteomes and remove special characters and spaces from the headers of the entries
prot=${WORKDIR}/ncbi_proteomes.fasta
prot2=${WORKDIR}/ncbi_proteomes_header_renamed.fasta

# uncomment and run if you need to
# conda create -y -n seqkit -c conda-forge -c bioconda seqkit 
# conda activate seqkit
# seqkit replace -p "[,\[\]\*]" -r "" "$prot" | seqkit replace -p "[ \|\t]" -r "_" > "$prot2"

# galba/braker will create directories with this name in your home directory
# these need to be unique for each run/genome
species_name="<your unique run id here>"

# for rather distant donors: disable_diamond_filter
# species name should instead be run-specific. Pri and alt need separate species names
singularity exec --bind $WORKDIR //path/to/containers/galba_v1.0.11.2.sif galba.pl \
--species=${species_name} \
--genome=${asm} \
--prot_seq=${prot2} \
--AUGUSTUS_ab_initio \
--threads ${THREADS}
