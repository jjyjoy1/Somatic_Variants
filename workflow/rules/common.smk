#!/usr/bin/env python
# Common rules and helper functions for the somatic pipeline
# Author: Jiyang Jiang
# Date: April 1, 2025

import os
from snakemake.utils import validate
import pandas as pd
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# Function to get sample type (tumor or normal)
def get_sample_type(wildcards):
    return samples.loc[wildcards.sample, "sample_type"]

# Function to determine if a sample has a matching normal
def has_normal(wildcards):
    if get_sample_type(wildcards) == "tumor":
        patient = samples.loc[wildcards.sample, "patient"]
        for s, row in samples.iterrows():
            if row["patient"] == patient and row["sample_type"] == "normal":
                return True
    return False

# Function to get matching normal sample for a tumor
def get_normal_sample(wildcards):
    if get_sample_type(wildcards) == "tumor":
        patient = samples.loc[wildcards.sample, "patient"]
        for s, row in samples.iterrows():
            if row["patient"] == patient and row["sample_type"] == "normal":
                return s
    return None

# Function to get all read files for a sample
def get_fastqs(wildcards):
    """Get all FASTQ files for a given sample."""
    sample_units = units.loc[wildcards.sample]
    fastqs = []
    
    if isinstance(sample_units, pd.DataFrame):
        for _, unit in sample_units.iterrows():
            if unit["fq2"]:
                fastqs.extend([unit["fq1"], unit["fq2"]])
            else:
                fastqs.append(unit["fq1"])
    else:
        if sample_units["fq2"]:
            fastqs.extend([sample_units["fq1"], sample_units["fq2"]])
        else:
            fastqs.append(sample_units["fq1"])
            
    return fastqs

# Function to get read 1 files for a sample
def get_fastq_r1(wildcards):
    """Get all R1 FASTQ files for a given sample."""
    sample_units = units.loc[wildcards.sample]
    fastqs = []
    
    if isinstance(sample_units, pd.DataFrame):
        for _, unit in sample_units.iterrows():
            fastqs.append(unit["fq1"])
    else:
        fastqs.append(sample_units["fq1"])
            
    return fastqs

# Function to get read 2 files for a sample
def get_fastq_r2(wildcards):
    """Get all R2 FASTQ files for a given sample."""
    sample_units = units.loc[wildcards.sample]
    fastqs = []
    
    if isinstance(sample_units, pd.DataFrame):
        for _, unit in sample_units.iterrows():
            if unit["fq2"]:
                fastqs.append(unit["fq2"])
    else:
        if sample_units["fq2"]:
            fastqs.append(sample_units["fq2"])
            
    return fastqs

# Function to check if a sample is paired-end
def is_paired_end(wildcards):
    """Check if a sample is paired-end."""
    sample_units = units.loc[wildcards.sample]
    
    if isinstance(sample_units, pd.DataFrame):
        return all(pd.notnull(sample_units["fq2"]))
    else:
        return pd.notnull(sample_units["fq2"])

# Function to get all BAM files for a given caller's input
def get_caller_input_bams(wildcards):
    result = {
        "tumor_bam": f"{OUTDIR}/mapping/{wildcards.sample}.sorted.md.bqsr.bam",
        "tumor_bai": f"{OUTDIR}/mapping/{wildcards.sample}.sorted.md.bqsr.bam.bai",
    }
    
    # Add normal BAM if available
    if has_normal(wildcards):
        normal_sample = get_normal_sample(wildcards)
        result["normal_bam"] = f"{OUTDIR}/mapping/{normal_sample}.sorted.md.bqsr.bam"
        result["normal_bai"] = f"{OUTDIR}/mapping/{normal_sample}.sorted.md.bqsr.bam.bai"
        
    return result

# Function to download reference genomes if needed
def get_reference_genome(wildcards):
    ref_info = REFERENCE_GENOMES[config["reference_genome"]]
    
    if not os.path.exists(ref_info["fasta"]):
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(ref_info["fasta"]), exist_ok=True)
        
        # Download reference genome
        shell(f"wget {ref_info['download_url']} -O {ref_info['fasta']}.gz")
        shell(f"gunzip {ref_info['fasta']}.gz")
        
        # Index reference genome
        shell(f"samtools faidx {ref_info['fasta']}")
        shell(f"picard CreateSequenceDictionary R={ref_info['fasta']} O={ref_info['dict']}")
        shell(f"bwa index {ref_info['fasta']}")
        
    return ref_info["fasta"]


