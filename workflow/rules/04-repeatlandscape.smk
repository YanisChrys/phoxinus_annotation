
rule repeatlandscape:
# repeat masker .out output is in gff2 format
    input:
        "results_" + HAP + "/01_masking/full_no_TPSI/genome.fas.align"
    output:
        divsum="results_" + HAP + "/01_masking/full_no_TPSI/genome.fas.divsum",
        landscape="results_" + HAP + "/01_masking/full_no_TPSI/genome.fas.html"
    threads: 
        workflow.cores
    container: 
        config["containers"]["repeat_masker"]
    params:
        genome_size=config["genome_size"]
    shell: """
        calcDivergenceFromAlign.pl \
        -s {output.divsum} {input}

        createRepeatLandscape.pl \
        -g {params.genome_size} \
        -div {output.divsum} > \
        {output.landscape}
    """


