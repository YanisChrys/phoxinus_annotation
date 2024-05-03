

rule hisat2_index:
    input:
        "results_" + HAP + "/01_masking/genome.masked.fa"
    output:
        directory("results_" + HAP + "/02_rnaseq/aligned/index")
    conda:
        "../envs/hisat2.yaml"
    threads: 
        min(workflow.cores,10)
    shell: """
        mkdir -p {output}
        hisat2-build {input} {output}/genome_index
    """

# dta is needed for braker -> transcript assembly
rule hisat2:
    input:
        R1=config["rnaseq"]["folder"] + "{sampleID}1." + extnsn,
        R2=config["rnaseq"]["folder"] + "{sampleID}2." + extnsn,
        # asm="results_" + HAP + "/01_masking/genome.masked.fa",
        h_index_dir="results_" + HAP + "/02_rnaseq/aligned/index"
    output:
        sam="results_" + HAP + "/02_rnaseq/aligned/{sampleID}_accepted_hits.sam",
        summary="results_" + HAP + "/02_rnaseq/aligned/{sampleID}_splicesite.txt"
    threads: 
        min(workflow.cores,10)
    conda:
        "../envs/hisat2.yaml"
    shell: """
        hisat2 \
        --dta \
        --novel-splicesite-outfile {output.summary} \
        -p {threads} \
        -x {input.h_index_dir}/genome_index \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {output.sam}
    """


rule sam_to_bam:
    input:
        "results_" + HAP + "/02_rnaseq/aligned/{sampleID}_accepted_hits.sam"
    output:
        bam=temp("results_" + HAP + "/02_rnaseq/aligned/{sampleID}_accepted_hits.sorted.bam"),
        bai="results_" + HAP + "/02_rnaseq/aligned/{sampleID}_accepted_hits.sorted.bam.bai"
    threads: 
        min(workflow.cores,5)
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools view -bS  {input} | samtools sort -@ {threads} -o {output.bam} - && \
        samtools index {output.bam}
    """


rule merge:
    input:
        expand("results_" + HAP + "/02_rnaseq/aligned/{sampleID}_accepted_hits.sorted.bam", sampleID=RNA_WLD.sampleID)
    output:
        "results_" + HAP + "/02_rnaseq/aligned/all_samples.bam"
    threads: 
        min(workflow.cores,20)
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools merge -@ {threads} {output} {input}

        samtools index {output}
    """

rule samtools_stats:
    input:
        "results_" + HAP + "/02_rnaseq/aligned/all_samples.bam"
    output:
        "results_" + HAP + "/02_rnaseq/aligned/stats/samtools_stats.txt"
    conda:
        "../envs/samtools.yaml"
    threads: 
        min(workflow.cores,10)
    shell: """
        samtools stats --threads {threads} {input} > {output} 
    """

#-hm 3 --collect-overlap-pairs -nr 1000 -nw 500 -outformat PDF:HTML
rule qualimap:
    input:
        "results_" + HAP + "/02_rnaseq/aligned/all_samples.bam"
    output:
        "results_" + HAP + "/02_rnaseq/aligned/stats/qualimap/report.pdf"
    conda:
        "../envs/qualimap.yaml"
    threads:
        min(workflow.cores,10)
    params:
        dir="results_" + HAP + "/02_rnaseq/aligned/stats/qualimap/"
    shell: """
        qualimap bamqc -bam {input} -c -outdir {params.dir} -outformat pdf -nt {threads} --java-mem-size=10G
    """
