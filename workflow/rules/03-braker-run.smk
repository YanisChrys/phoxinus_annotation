### BRAKER3 ###


# input files should be within the directory tree of the working directory 
# (ie directory from which the snakemake command is called)
rule braker3:
    input:
        **conditional_braker_input(config["BRAKER3"]["mode"])
    output:
        "results_" + HAP + "/04_braker3/braker.gff3",
        "results_" + HAP + "/04_braker3/braker.aa",
        "results_" + HAP + "/04_braker3/augustus.hints.gff",
        "results_" + HAP + "/04_braker3/hintsfile.gff",
        "results_" + HAP + "/04_braker3/braker.gtf"
    threads:
        min(workflow.cores,45)
    params:
        species=config["augustus"]["species"],
        outdir="results_" + HAP + "/04_braker3/",
        cfg=config["augustus"]["cfg"],
        augcfg=config["augustus"]["cfg"],
        augconfig="results_" + HAP + "/04_braker3/config",
        augscripts="results_" + HAP + "/04_braker3/scripts",
        input_options=lambda wildcards, input: braker_options(config["BRAKER3"]["mode"],input),
        user_specified_options=config["BRAKER3"]["user_options"]
    container:
        config["containers"]["braker3"]
    shell: """      
        # copy augustus config directory to working dir 
        # in the braker 3 container, the augustus config file is here:
        # /usr/share/augustus/config
        # 
        # sed 's|bin/augustus$|config|'
        # for module:

        #/usr/share/augustus/scripts/
        # bin/augustus -> config/


        cp -rn $(which augustus | sed 's|bin/augustus$|share/augustus/config|') {params.outdir}
        cp -rn $(which augustus | sed 's|bin/augustus$|share/augustus/scripts|') {params.outdir}
        cp {params.cfg} ./
        
        sleep 100

        braker.pl \
        --threads={threads} \
        {params.input_options} \
        --useexisting --species={params.species} \
        --workingdir={params.outdir} \
        {params.user_specified_options} \
        --AUGUSTUS_SCRIPTS_PATH=$PWD/{params.augscripts} \
        --AUGUSTUS_CONFIG_PATH=$PWD/{params.augconfig}
    """

rule braker3busco:
    input:
        "results_" + HAP + "/04_braker3/braker.aa"
    output:
        "results_" + HAP + "/04_braker3/busco/run_" + config["busco"]["lineage"] + "/missing_busco_list.tsv"
    threads:
        workflow.cores
    params:
        dataset_dir=config["busco"]["directory"],
        out_dir="results_" + HAP + "/04_braker3/",
        run_name='busco',
        lineage=config["busco"]["lineage"]
    conda:
        "../envs/busco.yaml"
    shell: """
        busco -f -m prot \
        -c {threads} \
        -i {input} \
        -o {params.run_name} \
        --download_path {params.dataset_dir} \
        --out_path {params.out_dir} \
        -l {params.lineage} \
        --offline
    """
