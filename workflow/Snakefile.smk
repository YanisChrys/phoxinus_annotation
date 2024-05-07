
### Load Configurations ###

import os
import re
import pandas as pd
from glob import glob

configfile: "config/config_copy.yaml"

HAP=config["genome"]["haplotype"]

# Ensure RNA seq data folder path ends with '/'
if config["rnaseq"]["folder"] and not config["rnaseq"]["folder"].endswith('/'):
    raise ValueError("RNA seq data folder path must end with '/', please check your configuration file")

### Wildcard Constraints ###

# 1) handle different extensions for reads
# 2) make sure the strandnumber(1,2) is understood as a number
# 3) ensure that the species wildcard does not contain any slashes("/")
wildcard_constraints:
    ext="(fastq.gz)|(fq.gz)",
    str="[0-9]",
    species="^((?!/).)*$"

### Constants ###

# find the names of our input and save the basename and extension inside wildcards
# if there are two R files for the same sample, this assumes they are symmetrical
rnafolder=config["rnaseq"]["folder"]
RNA_WLD=glob_wildcards(rnafolder + "{sampleID}1.{ext}")

#extract the extension, given that all files have the same one
#appending this to end of our files should take care of the different 
#fastq extensions
extnsn=list(RNA_WLD.ext).pop()

# if input stranded there will be two outputs

strandnmb = ["1"] if config["rnaseq"]["type"] == "Unstranded" else ["1", "2"]

###########################

## Rule All

rule all:
    input:
    # masking
        "results_" + HAP + "/01_masking/all_repeats.bed",
        "results_" + HAP + "/01_masking/genome.masked.fa",
        "results_" + HAP + "/02_rnaseq/aligned/stats/qualimap/report.pdf",
        "results_" + HAP + "/01_masking/full/genome.fas.divsum",
        "results_" + HAP + "/01_masking/full/genome.fas.html",
    # braker
        "results_" + HAP + "/04_braker3/busco/run_" + config["busco"]["lineage"] + "/missing_busco_list.tsv"

### Functions ###

def strandfinder(strand_type, rnaseq_type):
    """
    Determine the appropriate strand option for wig2hints based on strand type and RNA-seq type.

    Parameters:
    - strand_type (str): The type of strand ("1" or "2").
    - rnaseq_type (str): The type of RNA-seq data ("Unstranded" or "Stranded").

    Returns:
    - str: The corresponding strand option (".", "+", or "-") for wig2hints.

    Notes:
    The function is designed to simplify the use of wig2hints by determining the correct strand option 
    based on the provided strand type and RNA-seq type.
    """
    if strand_type == "1" and rnaseq_type == "Unstranded":
        return "."
    elif strand_type == "1" and rnaseq_type == "Stranded":
        return "+"
    elif strand_type == "2" and rnaseq_type == "Stranded":
        return "-"

def librarytype(seqtype):
    """
    Determine the library type based on the provided RNA-seq type.
    """
    if seqtype == "Unstranded":
        return "fr-unstranded"
    elif seqtype == "Stranded":
        return "fr-firststrand"

def conditional_braker_input(mode):
    """
    Generate input paths for BRAKER3 based on the user-specified Genemark mode.
    
    Parameters:
    - mode (str): Mode of operation. Can be one of "ES", "ET", "ET_hints", "EP", or "ETP".
    
    Returns:
    - dict: A dictionary containing the paths for the required inputs based on the mode.
    """
    inputs = {}
    if mode == "ES":
    # genome and rna seq
        inputs["genome"]="results_" + HAP + "/01_masking/genome.masked.fa"
    elif mode == "ET":
    # genome and rna seq
        inputs["genome"]="results_" + HAP + "/01_masking/genome.masked.fa"
        inputs["rnaseq_bam"] = "results_" + HAP + "/02_rnaseq/aligned/all_samples.bam"
    elif mode == "ET_hints":
    # genome and rna hints
        inputs["genome"]="results_" + HAP + "/01_masking/genome.masked.fa"
        inputs["rnaseq_hints"] = "results_" + HAP + "/02_rnaseq/aligned/*.hints"
    elif mode == "EP":
    # genome and proteins
        inputs["genome"]="results_" + HAP + "/01_masking/genome.masked.fa"
        inputs["prot_seq"] = config["prot_db"]
    elif mode == "ETP":
    # genome, protein and RNA sequences
        inputs["genome"]="results_" + HAP + "/01_masking/genome.masked.fa"
        inputs["rnaseq_bam"] = "results_" + HAP + "/02_rnaseq/aligned/all_samples.bam"
        inputs["prot_seq"] = config["prot_db"]
    else:
        raise ValueError("Invalid mode. Please select from ET, ET_hints, EP, or ETP.")
    return inputs


def braker_options(mode,inputs):
    """
    Generate commandline commands for handling the different BRAKER3 input files
    based on the user-specified Genemark mode.
    
    Parameters:
    - mode (str): Mode of operation. Can be one of "ES", "ET", "ET_hints", "EP", or "ETP".
    
    Returns:
    - str: A string containing the command line options for BRAKER based on the mode.
    """
    options = []
    if mode == "ES":
        options.append(f"--genome={inputs['genome']}")
    elif mode == "ET":
        options.append(f"--genome={inputs['genome']}")
        options.append(f"--bam={inputs['rnaseq_bam']}")
    elif mode == "ET_hints":
        options.append(f"--genome={inputs['genome']}")
        options.append(f"--hints={inputs['rnaseq_hints']}")
    elif mode == "EP":
        options.append(f"--genome={inputs['genome']}")
        options.append(f"--prot_seq={inputs['prot_seq']}")
    elif mode == "ETP":
        options.append(f"--genome={inputs['genome']}")
        options.append(f"--bam={inputs['rnaseq_bam']}")
        options.append(f"--prot_seq={inputs['prot_seq']}")
    else:
        raise ValueError("Invalid mode. Please select from ET, EP, or ETP.")
    return options


include: "rules/01-repeat-masking.smk"


if config["rnaseq"]["mapping"] == "HISAT":
    #use hisat for mapping reads
    include: "rules/02-rnaseq-mapping-hisat.smk"
elif config["rnaseq"]["mapping"] == "STAR":
    # use star for mapping reads
    include: "rules/02-rnaseq-mapping.smk"
else:
    # Don't rock the boat...
    raise ValueError(f"Invalid mapping algorithm specified: Please use 'HISAT' or 'STAR'.")

include: "rules/03-braker-run.smk"
include: "rules/04-repeatlandscape.smk"