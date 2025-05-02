#!/usr/bin/env python
# Snakemake pipeline for somatic variant calling
# Supports both tumor-only and tumor-normal paired samples
# Author: Jiyang Jiang
# Date: April 1, 2025

import os
import glob
from snakemake.utils import min_version

# Set minimum Snakemake version
min_version("7.18.0")

# Load configuration
configfile: "config/config.yaml"

# Define reference genomes
REFERENCE_GENOMES = {
    "hg38": {
        "fasta": "resources/references/hg38/GRCh38.d1.vd1.fa",
        "dict": "resources/references/hg38/GRCh38.d1.vd1.dict",
        "download_url": "https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834"
    },
    "hg19": {
        "fasta": "resources/references/hg19/ucsc.hg19.fasta",
        "dict": "resources/references/hg19/ucsc.hg19.dict",
        "download_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    }
}

# Get reference genome from config
REFERENCE = REFERENCE_GENOMES[config["reference_genome"]]["fasta"]
REFERENCE_DICT = REFERENCE_GENOMES[config["reference_genome"]]["dict"]

# Import sample information
samples = config["samples"]
units = config["units"]

# Get variant callers to use
VARIANT_CALLERS = config.get("variant_callers", ["mutect2"])
VARIANT_ANNOTATORS = config.get("variant_annotators", ["funcotator"])

# Get output directory
OUTDIR = config.get("output_directory", "results")

# Determine if we need to merge VCFs
MERGE_VCFS = len(VARIANT_CALLERS) > 2

# Wildcard constraints
wildcard_constraints:
    sample="|".join(samples.keys()),
    caller="|".join(VARIANT_CALLERS),
    annotator="|".join(VARIANT_ANNOTATORS)

# Include rules
include: "workflow/rules/common.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/mapping.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/calling.smk"
include: "workflow/rules/filtering.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/report.smk"

# Define final output files
def get_final_output():
    final_output = []
    
    # QC reports
    final_output.append(f"{OUTDIR}/qc/multiqc_report.html")
    
    # BAM QC reports
    final_output.extend(expand(f"{OUTDIR}/qc/bam/{{sample}}/qualimap_report.html", 
                       sample=samples.keys()))
    
    # VCF files per sample and caller
    for sample in samples.keys():
        for caller in VARIANT_CALLERS:
            final_output.append(f"{OUTDIR}/variants/{sample}.{caller}.filtered.vcf.gz")
            final_output.append(f"{OUTDIR}/qc/vcf/{sample}.{caller}.bcftools_stats.txt")
    
    # Merged VCF files if more than 2 callers
    if MERGE_VCFS:
        final_output.extend(expand(f"{OUTDIR}/variants/{{sample}}.merged.vcf.gz", 
                           sample=samples.keys()))
                           
    # Annotated VCF files
    for sample in samples.keys():
        for annotator in VARIANT_ANNOTATORS:
            if MERGE_VCFS:
                final_output.append(f"{OUTDIR}/variants/{sample}.merged.{annotator}.vcf.gz")
            else:
                for caller in VARIANT_CALLERS:
                    final_output.append(f"{OUTDIR}/variants/{sample}.{caller}.{annotator}.vcf.gz")
    
    return final_output

# Target rule
rule all:
    input:
        get_final_output()


