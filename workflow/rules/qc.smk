#!/usr/bin/env python
# QC rules for the somatic pipeline
# Author: Jiyang Jiang
# Date: April 1, 2025

rule fastqc:
    input:
        get_fastqs
    output:
        html=expand("{outdir}/qc/fastqc/{{sample}}/{{sample}}_{read}_fastqc.html", 
                   outdir=OUTDIR, read=["R1", "R2"]),
        zip=expand("{outdir}/qc/fastqc/{{sample}}/{{sample}}_{read}_fastqc.zip", 
                  outdir=OUTDIR, read=["R1", "R2"])
    params:
        outdir=f"{OUTDIR}/qc/fastqc/{{sample}}"
    log:
        f"{OUTDIR}/logs/fastqc/{{sample}}.log"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input} 2> {log}
        """

rule multiqc:
    input:
        expand("{outdir}/qc/fastqc/{sample}/{sample}_{read}_fastqc.zip", 
              outdir=OUTDIR, sample=samples.keys(), read=["R1", "R2"]),
        expand("{outdir}/qc/bam/{sample}/qualimap_report.html", 
              outdir=OUTDIR, sample=samples.keys()),
        expand("{outdir}/qc/vcf/{sample}.{caller}.bcftools_stats.txt", 
              outdir=OUTDIR, sample=samples.keys(), caller=VARIANT_CALLERS)
    output:
        report=f"{OUTDIR}/qc/multiqc_report.html"
    params:
        outdir=f"{OUTDIR}/qc",
        config="config/multiqc_config.yaml"
    log:
        f"{OUTDIR}/logs/multiqc/multiqc.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc --force --outdir {params.outdir} --config {params.config} {OUTDIR} 2> {log}
        """

# BAM QC with Qualimap
rule bam_qc_qualimap:
    input:
        bam=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bqsr.bam"
    output:
        report=f"{OUTDIR}/qc/bam/{{sample}}/qualimap_report.html"
    params:
        outdir=f"{OUTDIR}/qc/bam/{{sample}}",
        memory="16G"
    log:
        f"{OUTDIR}/logs/qualimap/{{sample}}.log"
    threads: 8
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        qualimap bamqc -bam {input.bam} \
            -outdir {params.outdir} \
            -nt {threads} \
            --java-mem-size={params.memory} \
            -c \
            -gd HUMAN \
            -outformat HTML 2> {log}
        """

# VCF QC with bcftools stats
rule vcf_qc_bcftools:
    input:
        vcf=f"{OUTDIR}/variants/{{sample}}.{{caller}}.filtered.vcf.gz",
        ref=REFERENCE
    output:
        stats=f"{OUTDIR}/qc/vcf/{{sample}}.{{caller}}.bcftools_stats.txt"
    params:
        regions=config.get("target_regions", ""),
        options=lambda wildcards: "--samples " + wildcards.sample
    log:
        f"{OUTDIR}/logs/bcftools/stats/{{sample}}.{{caller}}.log"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        bcftools stats {params.options} \
            {'-t ' + params.regions if params.regions else ''} \
            --fasta-ref {input.ref} \
            {input.vcf} > {output.stats} 2> {log}
        """

# Create custom VCF QC report
rule vcf_qc_report:
    input:
        stats=expand("{outdir}/qc/vcf/{sample}.{caller}.bcftools_stats.txt", 
                    outdir=OUTDIR, sample="{sample}", caller=VARIANT_CALLERS)
    output:
        report=f"{OUTDIR}/qc/vcf/{{sample}}.report.html"
    params:
        sample="{sample}",
        callers=VARIANT_CALLERS,
        outdir=f"{OUTDIR}/qc/vcf"
    log:
        f"{OUTDIR}/logs/vcf_report/{{sample}}.log"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/vcf_report.R"

# Tumor Mutational Burden calculation
rule calculate_tmb:
    input:
        vcf=f"{OUTDIR}/variants/{{sample}}.{{caller}}.filtered.vcf.gz",
        regions=config.get("target_regions", "")
    output:
        tmb=f"{OUTDIR}/qc/tmb/{{sample}}.{{caller}}.tmb.txt"
    params:
        sample="{sample}",
        caller="{caller}",
        genome_size=config.get("effective_genome_size", 2800000000)
    log:
        f"{OUTDIR}/logs/tmb/{{sample}}.{{caller}}.log"
    conda:
        "../envs/variants.yaml"
    script:
        "../scripts/calculate_tmb.py"
