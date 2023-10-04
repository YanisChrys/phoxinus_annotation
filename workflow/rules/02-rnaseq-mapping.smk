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
