#!/usr/bin/env python
# Report generation rules for the somatic pipeline
# Author: Jiyang Jiang
# Date: April 1, 2025

# Generate report for a single sample
rule sample_report:
    input:
        bam_qc=f"{OUTDIR}/qc/bam/{{sample}}/qualimap_report.html",
        vcf_qc=f"{OUTDIR}/qc/vcf/{{sample}}.report.html",
        tmb=lambda wildcards: expand(f"{OUTDIR}/qc/tmb/{{sample}}.{{caller}}.tmb.txt", 
                                    sample=wildcards.sample, caller=VARIANT_CALLERS),
        vcfs=lambda wildcards: expand(f"{OUTDIR}/variants/{{sample}}.{{caller}}.filtered.vcf.gz", 
                                     sample=wildcards.sample, caller=VARIANT_CALLERS),
        merged_vcf=f"{OUTDIR}/variants/{{sample}}.merged.vcf.gz" if MERGE_VCFS else [],
        annotations=lambda wildcards: expand(
            f"{OUTDIR}/variants/{{sample}}.{'{caller}' if not MERGE_VCFS else 'merged'}.{{annotator}}.vcf.gz",
            sample=wildcards.sample, caller=VARIANT_CALLERS[0], annotator=VARIANT_ANNOTATORS)
    output:
        report=f"{OUTDIR}/reports/{{sample}}/report.html"
    params:
        sample="{sample}",
        callers=VARIANT_CALLERS,
        annotators=VARIANT_ANNOTATORS,
        outdir=f"{OUTDIR}/reports/{{sample}}"
    log:
        f"{OUTDIR}/logs/report/{{sample}}.log"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/generate_report.R"

# Generate project-level report
rule project_report:
    input:
        sample_reports=expand(f"{OUTDIR}/reports/{{sample}}/report.html", sample=samples.keys()),
        multiqc=f"{OUTDIR}/qc/multiqc_report.html"
    output:
        report=f"{OUTDIR}/reports/project_report.html",
        csv=f"{OUTDIR}/reports/variant_summary.csv"
    params:
        samples=samples.keys(),
        callers=VARIANT_CALLERS,
        annotators=VARIANT_ANNOTATORS,
        outdir=f"{OUTDIR}/reports"
    log:
        f"{OUTDIR}/logs/report/project.log"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/project_report.R"

# Generate SVG workflow diagram
rule workflow_diagram:
    output:
        svg=f"{OUTDIR}/reports/workflow.svg"
    params:
        configfile="config/config.yaml"
    log:
        f"{OUTDIR}/logs/report/workflow_diagram.log"
    conda:
        "../envs/snakemake.yaml"
    shell:
        """
        snakemake --configfile {params.configfile} \
            --forceall \
            --dag | \
        dot -Tsvg > {output.svg} 2> {log}
        """

# Generate variant comparison across callers
rule variant_comparison:
    input:
        vcfs=lambda wildcards: expand(f"{OUTDIR}/variants/{{sample}}.{{caller}}.filtered.vcf.gz", 
                                     sample=wildcards.sample, caller=VARIANT_CALLERS)
    output:
        venn=f"{OUTDIR}/reports/{{sample}}/caller_comparison.pdf",
        table=f"{OUTDIR}/reports/{{sample}}/caller_stats.tsv"
    params:
        sample="{sample}",
        callers=VARIANT_CALLERS,
        outdir=f"{OUTDIR}/reports/{{sample}}"
    log:
        f"{OUTDIR}/logs/report/comparison/{{sample}}.log"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/variant_comparison.R"

# Generate mutational signatures report
rule mutational_signatures:
    input:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{VARIANT_CALLERS[0]}.filtered.vcf.gz",
        ref=REFERENCE
    output:
        plot=f"{OUTDIR}/reports/{{sample}}/mutational_signatures.pdf",
        table=f"{OUTDIR}/reports/{{sample}}/signature_contribution.tsv"
    params:
        sample="{sample}",
        outdir=f"{OUTDIR}/reports/{{sample}}"
    log:
        f"{OUTDIR}/logs/report/signatures/{{sample}}.log"
    conda:
        "../envs/signatures.yaml"
    script:
        "../scripts/mutational_signatures.R"


