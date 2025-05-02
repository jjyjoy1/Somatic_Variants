#!/usr/bin/env python
# Variant calling rules for the somatic pipeline
# Author: Jiyang Jiang
# Date: April 1, 2025

# MuTect2 caller
rule mutect2:
    input:
        unpack(get_caller_input_bams),
        ref=REFERENCE,
        intervals=f"{OUTDIR}/resources/interval_list.interval_list",
        germline_resource=config.get("germline_resource", ""),
        panel_of_normals=config.get("panel_of_normals", "")
    output:
        vcf=f"{OUTDIR}/variants/{{sample}}.mutect2.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.mutect2.vcf.gz.tbi",
        stats=f"{OUTDIR}/variants/{{sample}}.mutect2.vcf.gz.stats",
        f1r2=f"{OUTDIR}/variants/{{sample}}.mutect2.f1r2.tar.gz"
    params:
        tumor_name=lambda wildcards: wildcards.sample,
        normal_name=lambda wildcards: get_normal_sample(wildcards) if has_normal(wildcards) else "",
        pon_arg=lambda wildcards, input: f"--panel-of-normals {input.panel_of_normals}" if input.panel_of_normals else "",
        germline_arg=lambda wildcards, input: f"--germline-resource {input.germline_resource}" if input.germline_resource else "",
        normal_arg=lambda wildcards: f"-normal {get_normal_sample(wildcards)}" if has_normal(wildcards) else "",
        java_opts="-Xmx16g"
    log:
        f"{OUTDIR}/logs/gatk/mutect2/{{sample}}.log"
    threads: 8
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java_opts}" Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} \
            -tumor {params.tumor_name} \
            {params.normal_arg} \
            {"-I " + input.normal_bam if has_normal(wildcards) else ""} \
            -L {input.intervals} \
            {params.pon_arg} \
            {params.germline_arg} \
            --f1r2-tar-gz {output.f1r2} \
            --native-pair-hmm-threads {threads} \
            -O {output.vcf} 2> {log}
        """

# VarDict caller
rule vardict:
    input:
        unpack(get_caller_input_bams),
        ref=REFERENCE,
        bed=config.get("target_regions", "")
    output:
        vcf=f"{OUTDIR}/variants/{{sample}}.vardict.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.vardict.vcf.gz.tbi"
    params:
        tumor_name=lambda wildcards: wildcards.sample,
        normal_name=lambda wildcards: get_normal_sample(wildcards) if has_normal(wildcards) else "",
        min_af=config.get("vardict_min_af", "0.01"),
        min_reads=config.get("vardict_min_reads", "4"),
        output_dir=directory(f"{OUTDIR}/variants"),
        bed=lambda wildcards, input: input.bed if input.bed else ""
    log:
        f"{OUTDIR}/logs/vardict/{{sample}}.log"
    threads: 8
    conda:
        "../envs/vardict.yaml"
    shell:
        """
        # Create bed file if not provided
        if [ -z "{params.bed}" ]; then
            TEMP_BED=$(mktemp)
            samtools idxstats {input.tumor_bam} | cut -f 1,2 | grep -v "*" | \
            awk '{{print $1"\\t1\\t"$2}}' > $TEMP_BED
            BED_FILE=$TEMP_BED
        else
            BED_FILE={params.bed}
        fi

        # Run VarDict
        if [ -n "{params.normal_name}" ]; then
            # Paired mode
            vardict-java -G {input.ref} -f {params.min_af} -N {params.tumor_name} \
                -b "{input.tumor_bam}|{input.normal_bam}" -c 1 -S 2 -E 3 -g 4 -th {threads} $BED_FILE | \
            testsomatic.R | \
            var2vcf_paired.pl -N "{params.tumor_name}|{params.normal_name}" -f {params.min_af} > {params.output_dir}/{wildcards.sample}.vardict.vcf
        else
            # Tumor-only mode
            vardict-java -G {input.ref} -f {params.min_af} -N {params.tumor_name} \
                -b {input.tumor_bam} -c 1 -S 2 -E 3 -g 4 -th {threads} $BED_FILE | \
            teststrandbias.R | \
            var2vcf_valid.pl -N {params.tumor_name} -E -f {params.min_af} > {params.output_dir}/{wildcards.sample}.vardict.vcf
        fi

        # Compress and index
        bgzip -f {params.output_dir}/{wildcards.sample}.vardict.vcf
        tabix -p vcf {params.output_dir}/{wildcards.sample}.vardict.vcf.gz

        # Clean up temporary bed if created
        if [ -z "{params.bed}" ]; then
            rm $TEMP_BED
        fi
        """

# FreeBayes caller
rule freebayes:
    input:
        unpack(get_caller_input_bams),
        ref=REFERENCE,
        regions=config.get("target_regions", "")
    output:
        vcf=f"{OUTDIR}/variants/{{sample}}.freebayes.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.freebayes.vcf.gz.tbi"
    params:
        min_af=config.get("freebayes_min_af", "0.01"),
        min_ao=config.get("freebayes_min_ao", "5"),
        pooled_discrete=config.get("freebayes_pooled_discrete", ""),
        region_arg=lambda wildcards, input: f"--target {input.regions}" if input.regions else "",
        normal_arg=lambda wildcards: f"--pooled-continuous" if has_normal(wildcards) else ""
    threads: 8
    log:
        f"{OUTDIR}/logs/freebayes/{{sample}}.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        """
        # Run FreeBayes
        freebayes \
            -f {input.ref} \
            --min-alternate-fraction {params.min_af} \
            --min-alternate-count {params.min_ao} \
            --pooled-continuous \
            --genotype-qualities \
            {params.region_arg} \
            {params.normal_arg} \
            -F 0.01 \
            --theta 0.01 \
            $(if [ -n "{input.normal_bam}" ]; then echo "{input.tumor_bam} {input.normal_bam}"; else echo "{input.tumor_bam}"; fi) | \
        bgzip -c > {output.vcf}

        # Index the VCF
        tabix -p vcf {output.vcf}
        """

# Lofreq caller
rule lofreq:
    input:
        tumor_bam=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bqsr.bam",
        tumor_bai=f"{OUTDIR}/mapping/{{sample}}.sorted.md.bqsr.bam.bai",
        ref=REFERENCE,
        regions=config.get("target_regions", "")
    output:
        vcf=f"{OUTDIR}/variants/{{sample}}.lofreq.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.lofreq.vcf.gz.tbi"
    params:
        min_af=config.get("lofreq_min_af", "0.01"),
        min_bq=config.get("lofreq_min_bq", "30"),
        region_arg=lambda wildcards, input: f"-l {input.regions}" if input.regions else "",
        output_prefix=f"{OUTDIR}/variants/{{sample}}.lofreq"
    threads: 8
    log:
        f"{OUTDIR}/logs/lofreq/{{sample}}.log"
    conda:
        "../envs/lofreq.yaml"
    shell:
        """
        # Run LoFreq
        lofreq call-parallel \
            --pp-threads {threads} \
            -f {input.ref} \
            -o {output.vcf}.tmp \
            -s \
            -S {input.ref}.fai \
            -a {params.min_af} \
            -b {params.min_bq} \
            {params.region_arg} \
            {input.tumor_bam} 2> {log}

        # Compress and index
        bgzip -c {output.vcf}.tmp > {output.vcf}
        tabix -p vcf {output.vcf}
        rm -f {output.vcf}.tmp
        """

# Strelka2 caller
rule strelka2:
    input:
        unpack(get_caller_input_bams),
        ref=REFERENCE,
        regions=config.get("target_regions", "")
    output:
        snvs=f"{OUTDIR}/variants/{{sample}}.strelka2.snvs.vcf.gz",
        indels=f"{OUTDIR}/variants/{{sample}}.strelka2.indels.vcf.gz",
        vcf=f"{OUTDIR}/variants/{{sample}}.strelka2.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.strelka2.vcf.gz.tbi"
    params:
        outdir=directory(f"{OUTDIR}/variants/{{sample}}.strelka2"),
        region_arg=lambda wildcards, input: f"--callRegions={input.regions}" if input.regions else "",
        is_exome=config.get("is_exome", "false"),
        exome_arg=lambda wildcards, params: "--exome" if params.is_exome == "true" else ""
    threads: 8
    log:
        f"{OUTDIR}/logs/strelka2/{{sample}}.log"
    conda:
        "../envs/strelka2.yaml"
    shell:
        """
        # Configure Strelka2
        mkdir -p {params.outdir}
        if [ -n "{input.normal_bam}" ]; then
            # Somatic mode
            configureStrelkaSomaticWorkflow.py \
                --normalBam {input.normal_bam} \
                --tumorBam {input.tumor_bam} \
                --referenceFasta {input.ref} \
                {params.region_arg} \
                {params.exome_arg} \
                --runDir {params.outdir} 2> {log}
        else
            # Germline mode for tumor-only samples
            configureStrelkaGermlineWorkflow.py \
                --bam {input.tumor_bam} \
                --referenceFasta {input.ref} \
                {params.region_arg} \
                {params.exome_arg} \
                --runDir {params.outdir} 2> {log}
        fi

        # Run Strelka2
        {params.outdir}/runWorkflow.py -m local -j {threads} 2>> {log}

        # Copy and rename output files
        if [ -n "{input.normal_bam}" ]; then
            # Somatic mode outputs
            cp {params.outdir}/results/variants/somatic.snvs.vcf.gz {output.snvs}
            cp {params.outdir}/results/variants/somatic.indels.vcf.gz {output.indels}
        else
            # Germline mode outputs
            cp {params.outdir}/results/variants/variants.vcf.gz {output.snvs}
            cp {params.outdir}/results/variants/variants.vcf.gz {output.indels}
        fi

        # Merge SNVs and indels
        bcftools concat -a {output.snvs} {output.indels} | \
        bcftools sort -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

# Samtools/VarScan2 caller
rule varscan2:
    input:
        unpack(get_caller_input_bams),
        ref=REFERENCE,
        regions=config.get("target_regions", "")
    output:
        snvs=f"{OUTDIR}/variants/{{sample}}.varscan2.snp.vcf",
        indels=f"{OUTDIR}/variants/{{sample}}.varscan2.indel.vcf",
        vcf=f"{OUTDIR}/variants/{{sample}}.varscan2.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.varscan2.vcf.gz.tbi"
    params:
        min_af=config.get("varscan_min_af", "0.01"),
        min_cov=config.get("varscan_min_coverage", "8"),
        p_value=config.get("varscan_p_value", "0.05"),
        strand_filter=config.get("varscan_strand_filter", "1"),
        region_arg=lambda wildcards, input: f"-l {input.regions}" if input.regions else "",
        output_prefix=f"{OUTDIR}/variants/{{sample}}.varscan2",
        java_opts="-Xmx8g"
    threads: 4
    log:
        f"{OUTDIR}/logs/varscan2/{{sample}}.log"
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        # Generate mpileup
        if [ -n "{input.normal_bam}" ]; then
            # Paired mode
            samtools mpileup -f {input.ref} {params.region_arg} -B -q 1 {input.normal_bam} {input.tumor_bam} > {params.output_prefix}.mpileup
            
            # Call variants
            java {params.java_opts} -jar $CONDA_PREFIX/share/varscan-*/VarScan.jar somatic \
                {params.output_prefix}.mpileup \
                {params.output_prefix} \
                --mpileup 1 \
                --min-coverage {params.min_cov} \
                --min-var-freq {params.min_af} \
                --p-value {params.p_value} \
                --strand-filter {params.strand_filter} \
                --output-vcf 1 2> {log}
                
            # Process somatic output
            java {params.java_opts} -jar $CONDA_PREFIX/share/varscan-*/VarScan.jar processSomatic \
                {output.snvs} --min-tumor-freq {params.min_af} 2>> {log}
            java {params.java_opts} -jar $CONDA_PREFIX/share/varscan-*/VarScan.jar processSomatic \
                {output.indels} --min-tumor-freq {params.min_af} 2>> {log}
                
            # Merge SNVs and INDELs
            bgzip -c {params.output_prefix}.snp.Somatic.hc.vcf > {params.output_prefix}.snp.Somatic.hc.vcf.gz
            bgzip -c {params.output_prefix}.indel.Somatic.hc.vcf > {params.output_prefix}.indel.Somatic.hc.vcf.gz
            bcftools concat -a {params.output_prefix}.snp.Somatic.hc.vcf.gz {params.output_prefix}.indel.Somatic.hc.vcf.gz | \
            bcftools sort -o {output.vcf}
        else
            # Tumor-only mode
            samtools mpileup -f {input.ref} {params.region_arg} -B -q 1 {input.tumor_bam} > {params.output_prefix}.mpileup
            
            # Call variants
            java {params.java_opts} -jar $CONDA_PREFIX/share/varscan-*/VarScan.jar mpileup2snp \
                {params.output_prefix}.mpileup \
                --min-coverage {params.min_cov} \
                --min-var-freq {params.min_af} \
                --p-value {params.p_value} \
                --strand-filter {params.strand_filter} \
                --output-vcf 1 > {output.snvs} 2> {log}
                
            java {params.java_opts} -jar $CONDA_PREFIX/share/varscan-*/VarScan.jar mpileup2indel \
                {params.output_prefix}.mpileup \
                --min-coverage {params.min_cov} \
                --min-var-freq {params.min_af} \
                --p-value {params.p_value} \
                --strand-filter {params.strand_filter} \
                --output-vcf 1 > {output.indels} 2>> {log}
                
            # Merge SNVs and INDELs
            bgzip -c {output.snvs} > {params.output_prefix}.snp.vcf.gz
            bgzip -c {output.indels} > {params.output_prefix}.indel.vcf.gz
            bcftools concat -a {params.output_prefix}.snp.vcf.gz {params.output_prefix}.indel.vcf.gz | \
            bcftools sort -o {output.vcf}
        fi

        # Index the merged VCF
        tabix -p vcf {output.vcf}
        
        # Cleanup
        rm -f {params.output_prefix}.mpileup
        """


# Merge VCF files from multiple callers if more than 2 are used
rule merge_vcfs:
    input:
        vcfs=lambda wildcards: expand(f"{OUTDIR}/variants/{{sample}}.{{caller}}.filtered.vcf.gz", 
                                      sample=wildcards.sample, caller=VARIANT_CALLERS),
        indexes=lambda wildcards: expand(f"{OUTDIR}/variants/{{sample}}.{{caller}}.filtered.vcf.gz.tbi", 
                                        sample=wildcards.sample, caller=VARIANT_CALLERS)
    output:
        vcf=f"{OUTDIR}/variants/{{sample}}.merged.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.merged.vcf.gz.tbi"
    params:
        min_callers=config.get("min_callers_agreement", 2),
        input_string=lambda wildcards, input: " ".join(input.vcfs)
    log:
        f"{OUTDIR}/logs/bcftools/merge/{{sample}}.log"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        # Merge VCFs with metadata about which caller identified each variant
        bcftools merge -m none --force-samples {params.input_string} | \
        bcftools annotate -x INFO/set | \
        bcftools norm -m+any | \
        bcftools view -i 'count({input_string}) >= {params.min_callers}' -O z -o {output.vcf} 2> {log}
        
        # Index merged VCF
        tabix -p vcf {output.vcf}
        """

