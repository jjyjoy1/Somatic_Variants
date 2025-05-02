#!/usr/bin/env python
# Mapping rules for the somatic pipeline
# Author: Jiyang Jiang
# Date: April 1, 2025

rule trim_reads_fastp:
    input:
        r1=get_fastq_r1,
        r2=get_fastq_r2
    output:
        r1=temp(f"{OUTDIR}/preprocessing/{{sample}}/trimmed_R1.fastq.gz"),
        r2=temp(f"{OUTDIR}/preprocessing/{{sample}}/trimmed_R2.fastq.gz"),
        html=f"{OUTDIR}/preprocessing/{{sample}}/fastp_report.html",
        json=f"{OUTDIR}/preprocessing/{{sample}}/fastp_report.json"
    params:
        adapter_r1=config.get("adapter_read1", "auto"),
        adapter_r2=config.get("adapter_read2", "auto"),
        min_qual=config.get("min_quality", 20),
        min_len=config.get("min_length", 50)
    threads: 8
    log:
        f"{OUTDIR}/logs/fastp/{{sample}}.log"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        fastp --in1 {input.r1[0]} \
            --in2 {input.r2[0]} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --adapter_sequence {params.adapter_r1} \
            --adapter_sequence_r2 {params.adapter_r2} \
            --qualified_quality_phred {params.min_qual} \
            --length_required {params.min_len} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --correction \
            --cut_right \
            --cut_right_window_size 4 \
            --cut_right_mean_quality 20 2> {log}
        """

rule bwa_mem_align:
    input:
        r1=f"{OUTDIR}/preprocessing/{{sample}}/trimmed_R1.fastq.gz",
        r2=f"{OUTDIR}/preprocessing/{{sample}}/trimmed_R2.fastq.gz",
        ref=REFERENCE
    output:
        temp(f"{OUTDIR}/mapping/{{sample}}.aligned.bam")
    params:
        read_group=lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA"
    threads: 16
    log:
        f"{OUTDIR}/logs/bwa/{{sample}}.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bwa mem -t {threads} \
            -R '{params.read_group}' \
            {input.ref} \
            {input.r1} {input.r2} 2> {log} | \
        samtools view -bS - > {output}
        """

rule sort_bam:
    input:
        f"{OUTDIR}/mapping/{{sample}}.aligned.bam"
    output:
        temp(f"{OUTDIR}/mapping/{{sample}}.sorted.bam")
    threads: 8
    log:
        f"{OUTDIR}/logs/samtools/sort/{{sample}}.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 2> {log}
        """

rule mark_duplicates:
    input:
        f"{OUTDIR}/mapping/{{sample}}.sorted.bam"
    output:
        bam=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bam",
        metrics=f"{OUTDIR}/mapping/{{sample}}.md_metrics.txt"
    params:
        validation_stringency=config.get("validation_stringency", "LENIENT"),
        remove_duplicates=config.get("remove_duplicates", "false"),
        java_opts="-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=4"
    log:
        f"{OUTDIR}/logs/picard/md/{{sample}}.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        picard {params.java_opts} MarkDuplicates \
            INPUT={input} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
            VALIDATION_STRINGENCY={params.validation_stringency} \
            REMOVE_DUPLICATES={params.remove_duplicates} \
            CREATE_INDEX=true 2> {log}
        """

rule index_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        f"{OUTDIR}/logs/samtools/index/{{prefix}}.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools index {input} {output} 2> {log}
        """

