#!/usr/bin/env python
# Preprocessing rules for the somatic pipeline
# Author: Jiyang JIang
# Date: April 1, 2025

rule create_recalibration_table:
    input:
        bam=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bam",
        bai=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bam.bai",
        ref=REFERENCE,
        known_sites=config["known_sites"]
    output:
        recal=f"{OUTDIR}/preprocessing/{{sample}}/recal_data.table"
    params:
        known_sites=lambda wildcards, input: [f"--known-sites {site}" for site in input.known_sites],
        java_opts="-Xmx16g"
    log:
        f"{OUTDIR}/logs/gatk/bqsr/{{sample}}.recal_table.log"
    envmodules:
        "somatics_variants_bundle/v0.1.0",
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=lambda wildcards, attempt: 60 * 4 * attempt
    shell:
        """
        gatk --java-options "{params.java_opts}" BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            {" ".join(params.known_sites)} \
            -O {output.recal} \
            --use-original-qualities 2> {log}
        """

rule apply_bqsr:
    input:
        bam=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bam",
        bai=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bam.bai",
        ref=REFERENCE,
        recal=f"{OUTDIR}/preprocessing/{{sample}}/recal_data.table"
    output:
        bam=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bqsr.bam"
    params:
        java_opts="-Xmx16g"
    log:
        f"{OUTDIR}/logs/gatk/bqsr/{{sample}}.apply_bqsr.log"
    envmodules:
        "somatics_variants_bundle/v0.1.0",
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=lambda wildcards, attempt: 60 * 4 * attempt
    shell:
        """
        gatk --java-options "{params.java_opts}" ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.recal} \
            -O {output.bam} \
            --create-output-bam-index true \
            --static-quantized-quals 10 \
            --static-quantized-quals 20 \
            --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --use-original-qualities 2> {log}
        """

rule create_interval_lists:
    input:
        ref_dict=REFERENCE_DICT,
        target_bed=config.get("target_regions", "")
    output:
        intervals=f"{OUTDIR}/resources/interval_list.interval_list"
    params:
        java_opts="-Xmx4g"
    log:
        f"{OUTDIR}/logs/picard/interval_list.log"
    envmodules:
        "somatics_variants_bundle/v0.1.0",
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=lambda wildcards, attempt: 60 * 4 * attempt
    shell:
        """
        if [ -f "{input.target_bed}" ]; then
            picard {params.java_opts} BedToIntervalList \
                I={input.target_bed} \
                O={output.intervals} \
                SD={input.ref_dict} 2> {log}
        else
            picard {params.java_opts} ScatterIntervalsByNs \
                R={input.ref_dict} \
                OT=ACGT \
                N=10000 \
                O={output.intervals} 2> {log}
        fi
        """

rule split_intervals:
    input:
        ref=REFERENCE,
        intervals=f"{OUTDIR}/resources/interval_list.interval_list"
    output:
        intervals=directory(f"{OUTDIR}/resources/scattered_intervals")
    params:
        scatter_count=config.get("scatter_count", 50),
        java_opts="-Xmx4g"
    log:
        f"{OUTDIR}/logs/gatk/split_intervals.log"
    envmodules:
        "somatics_variants_bundle/v0.1.0",
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=lambda wildcards, attempt: 60 * 4 * attempt
    shell:
        """
        mkdir -p {output.intervals}
        gatk --java-options "{params.java_opts}" SplitIntervals \
            -R {input.ref} \
            -L {input.intervals} \
            --scatter-count {params.scatter_count} \
            -O {output.intervals} 2> {log}
        """

