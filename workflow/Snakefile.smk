
### Load Configurations ###

import os
import re
import pandas as pd
from glob import glob

configfile: "config/config_HAP1.yaml"

HAP=config["genome"]["haplotype"]

# Ensure RNA seq data folder path ends with '/'
if not config["rnaseq"]["folder"].endswith('/'):
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
        rules.final_masking.output,
    # rna seq
        expand("results_" + HAP + "/02_rnaseq/aligned/hints.rnaseq.exonpart{str}.gff", str=strandnmb),
        expand("results_" + HAP + "/02_rnaseq/cufflinks/{sampleID}/transcripts.gtf", sampleID=RNA_WLD.sampleID),
    # braker
        expand("results_" + HAP + "/04_braker3/busco/run_{db}/missing_busco_list.tsv", db=config["busco"]["lineage"]),

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

include: "rules/01-repeat-masking.smk"
include: "rules/02-rnaseq-mapping.smk"
include: "rules/03-braker-run.smk"
