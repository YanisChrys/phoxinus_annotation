### RNA-seq mapping ###

# Programs and Descriptions in current file:

# - Create an index for the genome for STAR alignment.
#   - Program: STAR.

# - Map RNA-Seq reads to the reference genome.
#   - Program: STAR.

# - Merge BAM files from different samples.
#   - Program: Samtools.

# - Convert BAM files to WIG format for visualization.
#   - Program: STAR.

# - Generate intron hints from RNA-Seq BAM files for gene prediction.
#   - Program: bam2hints (part of Augustus suite).

# - Convert WIG files to exonpart hints for gene prediction.
#   - Program: wig2hints.pl (part of Augustus suite).

# - Assemble transcripts using RNA-Seq BAM files.
#   - Program: Cufflinks.

rule star_index:
    input:
        "results_" + HAP + "/01_masking/genome.masked.fa"
    output:
        "results_" + HAP + "/02_rnaseq/STAR_index/chrLength.txt",
        "results_" + HAP + "/02_rnaseq/STAR_index/chrNameLength.txt",
        "results_" + HAP + "/02_rnaseq/STAR_index/chrName.txt",
        "results_" + HAP + "/02_rnaseq/STAR_index/chrStart.txt",
        "results_" + HAP + "/02_rnaseq/STAR_index/Genome",
        "results_" + HAP + "/02_rnaseq/STAR_index/genomeParameters.txt",
        "results_" + HAP + "/02_rnaseq/STAR_index/Log.out",
        "results_" + HAP + "/02_rnaseq/STAR_index/SA",
        "results_" + HAP + "/02_rnaseq/STAR_index/SAindex"
    threads:
        min(workflow.cores,10)
    envmodules:
        config["modules"]["star"]
    params:
        genomedir="results_" + HAP + "/02_rnaseq/STAR_index/"
    shell: """
        STAR --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {params.genomedir} \
        --genomeSAindexNbases 13 \
        --genomeFastaFiles {input} --outFileNamePrefix starindex
    """

rule star_mapping:
    input:
        R1=config["rnaseq"]["folder"] + "{sampleID}1." + extnsn,
        R2=config["rnaseq"]["folder"] + "{sampleID}2." + extnsn,
        index="results_" + HAP + "/02_rnaseq/STAR_index/chrLength.txt"
    output:
        "results_" + HAP + "/02_rnaseq/aligned/{sampleID}/Aligned.sortedByCoord.out.bam"
    threads:
        min(workflow.cores,10)
    params:
        genomedir="results_" + HAP + "/02_rnaseq/STAR_index/",
        prefix="results_" + HAP + "/02_rnaseq/aligned/{sampleID}/"
    envmodules:
        config["modules"]["star"]
    shell: """
        STAR \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.R1},{input.R2} \
        --outFileNamePrefix {params.prefix} \
        --readFilesCommand zcat \
        --runThreadN {threads} \
        --limitBAMsortRAM 32212254720 \
        --outSAMtype BAM SortedByCoordinate --outWigType wiggle \
        --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical \
        --alignSoftClipAtReferenceEnds No --outWigStrand Unstranded --outSAMstrandField intronMotif
    """

rule merge:
    input:
        expand("results_" + HAP + "/02_rnaseq/aligned/{sampleID}/Aligned.sortedByCoord.out.bam", sampleID=RNA_WLD.sampleID)
    output:
        "results_" + HAP + "/02_rnaseq/aligned/all_samples.bam"
    threads:
        workflow.cores
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools merge -@ {threads} {output} {input}
    """

rule bam2wig:
    input:
        "results_" + HAP + "/02_rnaseq/aligned/all_samples.bam"
    output:
        "results_" + HAP + "/02_rnaseq/aligned/Signal.Unique.str{str}.out.wig"
    threads:
        min(workflow.cores,10)
    params:
        fastqtype=config["rnaseq"]["type"],
        prefix="results_" + HAP + "/02_rnaseq/aligned/"
    envmodules:
        config["modules"]["star"]
    shell: """
        STAR --runMode inputAlignmentsFromBAM \
        --inputBAMfile {input} \
        --outFileNamePrefix {params.prefix} \
        --outWigType wiggle \
        --outWigStrand {params.fastqtype}
    """

rule wig2hints:
    input:
        "results_" + HAP + "/02_rnaseq/aligned/Signal.Unique.str{str}.out.wig"
    output:
        "results_" + HAP + "/02_rnaseq/aligned/hints.rnaseq.exonpart{str}.gff"
    params:
        prefix="results_" + HAP + "/02_rnaseq/aligned/final",
        extra="--width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --radius=4.5 --pri=4",
        strand = lambda wildcards: strandfinder('{strand_type}'.format(strand_type=wildcards.str),'{rnaseq_type}'.format(rnaseq_type=config["rnaseq"]["type"]))
    envmodules:
        config["modules"]["augustus"]
    shell: """
        wig2hints.pl {params.extra} --strand={params.strand} < {input} >  {output}
    """

rule cufflinks:
    input:
        "results_" + HAP + "/02_rnaseq/aligned/{sampleID}/Aligned.sortedByCoord.out.bam"
    output:
        "results_" + HAP + "/02_rnaseq/cufflinks/{sampleID}/transcripts.gtf"
    threads:
        workflow.cores
    params:
        outdir="results_" + HAP + "/02_rnaseq/cufflinks/{sampleID}/",
        extra="--max-frag-multihits 10 --max-bundle-length 5000000 --max-bundle-frags 1000000 --max-intron-length 100000 --overlap-radius 5",
        strand=librarytype(config["rnaseq"]["type"])
    conda:
        "../envs/genofish_cufflinks.yaml"
    shell: """
        cufflinks \
        --no-update-check \
        {params.extra} \
        --num-threads {threads} \
        --library-type {params.strand} \
        -o {params.outdir} {input}
    """
