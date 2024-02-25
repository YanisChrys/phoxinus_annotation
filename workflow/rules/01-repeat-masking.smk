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
        workflow.cores
    container: 
        config["containers"]["repeat_modeler"]
    log: 
        "results_" + HAP + "/logs/repeat_modeler.log"
    shell: """
        (cd {params.work}
        
        BuildDatabase -engine {params.engine} -name {params.database} genome.fas

        RepeatModeler -database {params.database} -engine {params.engine} -pa {threads} ) &> {log}
    """

rule repeat_masker_self:
    input:
        "results_" + HAP + "/01_masking/MyDatabase-families.fa"
    output:
        "results_" + HAP + "/01_masking/self/genome.fas.masked",
        "results_" + HAP + "/01_masking/self/genome.fas.out"
    params:
        database="MyDatabase",
        workdir="results_" + HAP + "/01_masking/self",
        engine="ncbi",
        run_hap=HAP
    threads: 
        workflow.cores
    container: 
        config["containers"]["repeat_masker"]
    log:
        "results_" + HAP + "/logs/repeat_masker_self.log"
    shell: """
        [ ! -d {params.workdir} ] && mkdir -p {params.workdir} # create working directory
        [ ! -f results_{params.run_hap}/01_masking/self/genome.fas ] && cp results_{params.run_hap}/01_masking/genome.fas results_{params.run_hap}/01_masking/self/genome.fas
        ( cd {params.workdir}
        RepeatMasker -pa {threads} -xsmall -lib ../../../{input} genome.fas ) &> {log}
    """

rule repeat_masker_full:
    input:
        "results_" + HAP + "/01_masking/self/genome.fas.masked"
    output:
        "results_" + HAP + "/01_masking/full/genome.fas.masked.masked",
        "results_" + HAP + "/01_masking/full/genome.fas.masked.out"
    params:
        library="-lib " + config["repeat_library"] if os.path.exists(config["repeat_library"]) else "-species " + config["repeat_library"],
        workdir="results_" + HAP + "/01_masking/full",
        engine="ncbi",
        run_hap=HAP
    threads: 
        workflow.cores
    container: 
        config["containers"]["repeat_masker"]
    log:
        "results_" + HAP + "/logs/repeat_masker_full.log"
    shell: """
        [ ! -d {params.workdir} ] && mkdir -p {params.workdir} # create working directory
        [ ! -f results_{params.run_hap}/01_masking/full/genome.fas ] && cp {input} results_{params.run_hap}/01_masking/full/
        
        ( cd {params.workdir}
        RepeatMasker -a -pa {threads} -xsmall {params.library} genome.fas.masked ) &> {log}
    """

rule repeat_masker2gff:
    input:
        rm_self="results_" + HAP + "/01_masking/self/genome.fas.out",
        rm_full="results_" + HAP + "/01_masking/full/genome.fas.masked.out"
    output:
        rm_self="results_" + HAP + "/01_masking/RMasker_self_masked.gff3",
        rm_full="results_" + HAP + "/01_masking/RMasker_full_masked.gff3"
    container: 
        config["containers"]["repeat_masker"]
    shell: """
        rmOutToGFF3.pl {input.rm_self} | \
        perl -ane '$id; if(!/^\#/){{@F=split(/\t/, $_);chomp $F[-1];$id++;$F[-1]="ID=$F[0]-$id;$F[-1]";$F[-1]=~s/ /%20/g;$_=join("\t", @F)."\n"}} print $_'\
        > {output.rm_self}   

        rmOutToGFF3.pl {input.rm_full} | \
        perl -ane '$id; if(!/^\#/){{@F=split(/\t/, $_);chomp $F[-1];$id++;$F[-1]="ID=$F[0]-$id;$F[-1]";$F[-1]=~s/ /%20/g;$_=join("\t", @F)."\n"}} print $_'\
        > {output.rm_full}   
    """

rule final_masking:
    input:
        gffs=("results_" + HAP + "/01_masking/trf.gff3",        
        "results_" + HAP + "/01_masking/dustmasker.gff3",
        "results_" + HAP + "/01_masking/RMasker_self_masked.gff3",
        "results_" + HAP + "/01_masking/RMasker_full_masked.gff3"),
        genome=config["genome"]["file"]
    output:
        bed="results_" + HAP + "/01_masking/all_repeats.bed",
        genome="results_" + HAP + "/01_masking/genome.masked.fa"
    conda:
        "../envs/bedtools.yaml"
    shell: """
        cat {input.gffs} | grep -v ">" | sort -k1,1 -k4,4n | bedtools merge -i - > {output.bed}
        cat {output.bed} | bedtools maskfasta -soft -fi {input.genome} -bed - -fo {output.genome}
    """
