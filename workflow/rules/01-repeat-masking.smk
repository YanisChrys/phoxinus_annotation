### Repeat Masking ###

# Programs and Descriptions in current file:

# - Split the genome into individual sequences.
#   - Program: Custom script (fasta_split_1.pl).

# - Mask low-complexity regions in the genome.
#   - Program: Dustmasker (from BLAST+).

# - Convert Dustmasker output to GFF3 format.
#   - Program: Custom script (repeat_to_gff.pl).

# - Identify tandem repeats in the genome.
#   - Program: Tandem Repeat Finder (TRF).

# - Convert TRF output to GFF3 format.
#   - Program: Custom script (repeat_to_gff.pl).

# - Build a repeat database and model repeats in the genome.
#   - Program: RepeatModeler.

# - Mask the genome using the self-built repeat database.
#   - Program: RepeatMasker.

# - Mask the genome using a full repeat database.
#   - Program: RepeatMasker.

# - Convert RepeatMasker output to GFF3 format.
#   - Program: rmOutToGFF3.pl (part of RepeatMasker suite).

# - Merge all repeat annotations and mask the genome.
#   - Program: Bedtools.


rule split_genome:
    input:
        config["genome"]["file"]
    output:
        dynamic("results_" + HAP + "/01_masking/split/{SingleSequence}.tfa")
    params:
        run_hap=HAP
    shell: """
        cd results_{params.run_hap}/01_masking/split
        ../../../workflow/scripts/fasta_split_1.pl {input} 1 entries
    """

## Repeat Masking ##

rule dustmasker:
    input:
        "results_" + HAP + "/01_masking/split/{SingleSequence}.tfa"
    output:
        "results_" + HAP + "/01_masking/split/{SingleSequence}.dust"
    conda:
        "../envs/blast.yaml"
    shell: """
       dustmasker -in {input} -out {output} -outfmt interval
    """

rule dustmasker2gff:
    input:
        dynamic("results_" + HAP + "/01_masking/split/{SingleSequence}.dust")
    output:
        dust="results_" + HAP + "/01_masking/dustmasker.dust",
        gff="results_" + HAP + "/01_masking/dustmasker.gff3"
    shell: """
        cat {input} > {output.dust}
        workflow/scripts/repeat_to_gff.pl {output.dust}

        mv {output.dust}.gff {output.gff}
    """

rule tandem_repeat_finder:
    input:
        "results_" + HAP + "/01_masking/split/{SingleSequence}.tfa"
    output:
        "results_" + HAP + "/01_masking/split/{SingleSequence}.tfa.2.7.7.80.10.50.500.dat"
    envmodules:
        config["modules"]["conda"],
        config["modules"]["trf"]
    params:
        run_hap=HAP
    shell: """
        cd results_{params.run_hap}/01_masking/split
        trf ../../../{input} 2 7 7 80 10 50 500 -d -h
    """

rule trf2gff:
    input:
        dynamic("results_" + HAP + "/01_masking/split/{SingleSequence}.tfa.2.7.7.80.10.50.500.dat")
    output:
        dat="results_" + HAP + "/01_masking/trf.dat",
        gff="results_" + HAP + "/01_masking/trf.gff3"
    shell: """
        cat {input} > {output.dat}
        workflow/scripts/repeat_to_gff.pl {output.dat}

        mv {output.dat}.gff {output.gff}
    """

rule copy_genome_locally:
    input:
        config["genome"]["file"]
    output:
        "results_" + HAP + "/01_masking/genome.fas"
    shell: """
        cp {input} {output}
    """

rule repeat_modeler:
    input:
        "results_" + HAP + "/01_masking/genome.fas"
    output:
        DB="results_" + HAP + "/01_masking/MyDatabase-families.fa"
    params:
        database="MyDatabase",
        work="results_" + HAP + "/01_masking",
        engine="ncbi"
    threads: 
        min(workflow.cores,20)
    container: 
        config["containers"]["repeat_modeler"]
    log: 
        "results_" + HAP + "/logs/repeat_modeler.log"
    shell: """
        (cd {params.work}
        
        BuildDatabase -engine {params.engine} -name {params.database} genome.fas

        RepeatModeler -database {params.database} -engine {params.engine} -pa {threads} ) &> {log}
    """


rule combine_dbs:
# combine repModeler and DFam database from inside RepeatMasker container
# Todo: handle absolute paths, and/or custom libraries
    input:
        "results_" + HAP + "/01_masking/MyDatabase-families.fa"
    output:
        species_db="results_" + HAP + "/01_masking/full/dfam_lib_" + config["repeat_library_name"] + ".fa",
        combined_lib="results_" + HAP + "/01_masking/full/custom_lib_" + config["repeat_library_name"] + ".fa"
    params:
        library="-lib " + config["repeat_library"] if os.path.exists(config["repeat_library"]) else "-species " + config["repeat_library"],
        workdir="results_" + HAP + "/01_masking/full/",
        engine="ncbi"
    threads: 
        min(workflow.cores,10)
    container: 
        config["containers"]["repeat_masker"]
    shell: """
        cd {params.workdir}
        queryRepeatDatabase.pl {params.library} > ../../../{output.species_db}
        cat ../../../{input} ../../../{output.species_db} > ../../../{output.combined_lib}
    """

### make sure no proteins overlap the repeats:

prot_filename = os.path.basename(config["prot_db"])
PROT_NAME, PROT_EXT = os.path.splitext(prot_filename)

rule split_uniprot:
    #Split the UniProt/Swissprot protein database into chunks for transposonPSI
    input: 
        config["prot_db"]
    output:
        chunks = temp(expand("results_" + HAP + "/01_masking/split_result/" + PROT_NAME + "_chunk{nr}.fa", nr=range(1, 101)))
    params:
        dir = "results_" + HAP + "/01_masking/split_result/",
        prot = PROT_NAME,
    conda: 
        "../envs/gaas.yaml"
    shell: """
        cd {params.dir}
        cp -p {input} ./

        gaas_fasta_splitter.pl -f ./{params.prot}.fa --nb_chunks 100 -o tmp/
        mv tmp/*fa ../split_result/
        rm -r tmp/
        rm ./{params.prot}.fa
    """

rule transposonPSI:
    #Identify transposons in the UniProt/Swissprot protein dataset#
    input:
        "results_" + HAP + "/01_masking/split_result/" + PROT_NAME + "_chunk{nr}.fa"
    output:
        allHits = temp("results_" + HAP + "/01_masking/split_result/" + PROT_NAME + "_chunk{nr}.fa.TPSI.allHits"),
        topHits = temp("results_" + HAP + "/01_masking/split_result/" + PROT_NAME + "_chunk{nr}.fa.TPSI.topHits")
    params:
        dir = "results_" + HAP + "/01_masking/split_result/"
    conda: 
        "../envs/tepsi.yaml"
    shell: """
        cd {params.dir}
        transposonPSI.pl ../../../{input} prot
    """


rule list_tePSI_hits:
    input:
        expand("results_" + HAP + "/01_masking/split_result/" + PROT_NAME + "_chunk{nr}.fa.TPSI.topHits", nr=range(1, 101))
    output:
        allTopHits = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".TPSI.topHits",
        prot_list = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".TPSI.topHits.accessions.txt"
    shell: """
        cat {input} > {output.allTopHits} &&
        awk '{{if($0 ~ /^[^\/\/.*]/) print $5}}' {output.allTopHits} | sort -u > {output.prot_list}
    """


rule filter_uniprot_fasta:
    #Remove transposons from the UniProt/Swissprot protein dataset
    input:
        prot = config["prot_db"],
        prot_list = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".TPSI.topHits.accessions.txt"
    output:
        prot_filtered = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa"
    params:
        dir = "results_" + HAP + "/01_masking/proteins_tepsi_results/",
        prot_path = "results_" + HAP + "/01_masking/proteins_tepsi_results/"
    conda: 
        "../envs/gaas.yaml"
    threads: 
        min(workflow.cores,20)
    shell: """
        cd {params.dir}
        gaas_fasta_removeSeqFromIDlist.pl -f {params.prot_path} -l {input.prot_list} -o ../../../{output.prot_filtered}
    """


rule filtered_blast_db:
    #Generate BLAST database from filtered UniProt/Swissprot protein dataset
    input: 
        prot_filtered = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa"
    output:
        phr = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa.phr",
        pin = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa.pin",
        psq = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa.psq"
    params:
        dir = "results_" + HAP + "/01_masking/proteins_tepsi_results/"
    conda: 
        "../envs/blast.yaml"
    threads: 
        min(workflow.cores,10)
    shell: """
        makeblastdb -in {input.prot_filtered} -dbtype prot
    """


rule blast_repeat_library:
    #Blastx repeat library to filtered Uniprot/Swissprot database
    input:
        repmo_raw = "results_" + HAP + "/01_masking/full/custom_lib_" + config["repeat_library_name"] + ".fa", 
        blast_db = config["prot_db"]
        # phr = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa.phr",
        # pin = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa.pin",
        # psq = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa.psq",
        # blast_db = "results_" + HAP + "/01_masking/proteins_tepsi_results/" + PROT_NAME + ".noTEs.fa"
    output:
        blast = "results_" + HAP + "/01_masking/full_no_TPSI/blast/custom_lib_" + config["repeat_library_name"] + ".fa" + "blastx.out"
    params:
        dir = "results_" + HAP + "/01_masking/full_no_TPSI/blast/"
    threads: 
        workflow.cores
    conda: 
        "../envs/blast.yaml"
    shell: """
        makeblastdb -in {input.blast_db} -dbtype prot
        blastx -num_threads {threads} -db {input.blast_db} -query {input.repmo_raw} -out {output.blast}
    """

rule protexcluder:
    #Remove blast hits from repeat library
    #script creates file in wd and adds the "noProtFinal" suffix with no flexibility
    input:
        blast = "results_" + HAP + "/01_masking/full_no_TPSI/blast/custom_lib_" + config["repeat_library_name"] + ".fa" + "blastx.out",
        repmo_raw = "results_" + HAP + "/01_masking/full/custom_lib_" + config["repeat_library_name"] + ".fa"
    output:
        "results_" + HAP + "/01_masking/full/custom_lib_" + config["repeat_library_name"] + ".fanoProtFinal"
    params:
        dir = "results_" + HAP + "/01_masking/full_no_TPSI/"
    conda: 
        "../envs/protexcluder.yaml"
    shell: """
        cd {params.dir}
        ProtExcluder.pl ../../../{input.blast} ../../../{input.repmo_raw}
    """

rule repeat_masker_final:
# repeat masker .out output is in gff2 format
    input:
        prot_db="results_" + HAP + "/01_masking/full/custom_lib_" + config["repeat_library_name"] + ".fanoProtFinal",
        genome="results_" + HAP + "/01_masking/genome.fas"
    output:
        "results_" + HAP + "/01_masking/full_no_TPSI/genome.fas.masked",
        "results_" + HAP + "/01_masking/full_no_TPSI/genome.fas.out",
        "results_" + HAP + "/01_masking/full_no_TPSI/genome.fas.align"
    params:
        database="MyDatabase",
        workdir="results_" + HAP + "/01_masking/full_no_TPSI", # no trailing slash
        engine="ncbi",
        run_hap=HAP
    threads: 
        workflow.cores
    container: 
        config["containers"]["repeat_masker"]
    log:
        "results_" + HAP + "/logs/repeat_masker_self.log"
    shell: """
        mkdir -p {params.workdir} # create working directory
        cp -p {input.genome} {params.workdir}/genome.fas
        ( cd {params.workdir}
        RepeatMasker -a -pa {threads} -xsmall -lib ../../../{input.prot_db} genome.fas ) &> {log}
    """


rule repeat_masker2gff:
    input:
        "results_" + HAP + "/01_masking/full_no_TPSI/genome.fas.out"
    output:
        "results_" + HAP + "/01_masking/full_no_TPSI/RMasker_custom_blasted_final_masked.gff3"
    container: 
        config["containers"]["repeat_masker"]
    shell: """
        rmOutToGFF3.pl {input} | \
        perl -ane '$id; if(!/^\#/){{@F=split(/\t/, $_);chomp $F[-1];$id++;$F[-1]="ID=$F[0]-$id;$F[-1]";$F[-1]=~s/ /%20/g;$_=join("\t", @F)."\n"}} print $_'\
        > {output}   
    """

rule final_masking:
    input:
        gffs=("results_" + HAP + "/01_masking/trf.gff3",        
        "results_" + HAP + "/01_masking/dustmasker.gff3",
        "results_" + HAP + "/01_masking/full_no_TPSI/RMasker_custom_blasted_final_masked.gff3"),
        genome=config["genome"]["file"]
    output:
        bed="results_" + HAP + "/01_masking/full_no_TPSI/all_repeats.bed",
        genome="results_" + HAP + "/01_masking/full_no_TPSI/genome.masked.fa"
    conda:
        "../envs/bedtools.yaml"
    threads:
        min(workflow.cores,10)
    shell: """
        cat {input.gffs} | grep -v ">" | sort -k1,1 -k4,4n | bedtools merge -i - > {output.bed}
        cat {output.bed} | bedtools maskfasta -soft -fi {input.genome} -bed - -fo {output.genome}
    """




