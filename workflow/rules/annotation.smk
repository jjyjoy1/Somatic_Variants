#!/usr/bin/env python
# Variant annotation rules for the somatic pipeline
# Author: Jiyang Jiang
# Date: April 1, 2025

# GATK Funcotator annotation
rule funcotator:
    input:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.filtered.vcf.gz",
        idx=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz.tbi" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.filtered.vcf.gz.tbi",
        ref=REFERENCE,
        data_sources=config["funcotator_data_sources"]
    output:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.funcotator.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.funcotator.vcf.gz",
        idx=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.funcotator.vcf.gz.tbi" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.funcotator.vcf.gz.tbi"
    params:
        ref_version=config["reference_genome"],
        java_opts="-Xmx16g",
        output_format="VCF",
        transcript_selection=config.get("funcotator_transcript_selection", "CANONICAL")
    log:
        lambda wildcards: f"{OUTDIR}/logs/gatk/funcotator/{wildcards.sample}.merged.log" if MERGE_VCFS else f"{OUTDIR}/logs/gatk/funcotator/{wildcards.sample}.{wildcards.caller}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java_opts}" Funcotator \
            --variant {input.vcf} \
            --reference {input.ref} \
            --ref-version {params.ref_version} \
            --data-sources-path {input.data_sources} \
            --output {output.vcf} \
            --output-file-format {params.output_format} \
            --transcript-selection-mode {params.transcript_selection} 2> {log}
        """

# VEP annotation
rule vep:
    input:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.filtered.vcf.gz",
        idx=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz.tbi" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.filtered.vcf.gz.tbi",
        cache=config["vep_cache"],
        fasta=REFERENCE
    output:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.vep.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.vep.vcf.gz",
        idx=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.vep.vcf.gz.tbi" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.vep.vcf.gz.tbi",
        html=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.vep.html" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.vep.html"
    params:
        species=config.get("vep_species", "homo_sapiens"),
        assembly=config.get("vep_assembly", "GRCh38"),
        plugins=config.get("vep_plugins", ""),
        plugin_args=config.get("vep_plugin_args", ""),
        extra=config.get("vep_extra_args", "")
    threads: 8
    log:
        lambda wildcards: f"{OUTDIR}/logs/vep/{wildcards.sample}.merged.log" if MERGE_VCFS else f"{OUTDIR}/logs/vep/{wildcards.sample}.{wildcards.caller}.log"
    conda:
        "../envs/vep.yaml"
    shell:
        """
        vep --input_file {input.vcf} \
            --output_file {output.vcf} \
            --vcf \
            --compress_output bgzip \
            --stats_file {output.html} \
            --fork {threads} \
            --species {params.species} \
            --assembly {params.assembly} \
            --cache \
            --dir_cache {input.cache} \
            --fasta {input.fasta} \
            --everything \
            --allele_number \
            --hgvs \
            --shifted_3prime \
            --af \
            --af_1kg \
            --af_gnomad \
            --check_existing \
            --pick_allele \
            --filter_common \
            --flag_pick \
            --symbol \
            --canonical \
            --biotype \
            --sift b \
            --polyphen b \
            --regulatory \
            --variant_class \
            --mane \
            {params.plugins} {params.plugin_args} \
            {params.extra} 2> {log}
            
        # Index the output VCF
        tabix -p vcf {output.vcf}
        """

# SnpEff annotation
rule snpeff:
    input:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.filtered.vcf.gz",
        idx=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz.tbi" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.filtered.vcf.gz.tbi",
        db=config["snpeff_db"]
    output:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.snpeff.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.snpeff.vcf.gz",
        idx=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.snpeff.vcf.gz.tbi" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.snpeff.vcf.gz.tbi",
        csv=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.snpeff.csv" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.snpeff.csv",
        html=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.snpeff.html" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.snpeff.html"
    params:
        java_opts="-Xmx8g",
        genome_version=config.get("snpeff_genome", "GRCh38.p13"),
        extra=config.get("snpeff_extra_args", "")
    log:
        lambda wildcards: f"{OUTDIR}/logs/snpeff/{wildcards.sample}.merged.log" if MERGE_VCFS else f"{OUTDIR}/logs/snpeff/{wildcards.sample}.{wildcards.caller}.log"
    conda:
        "../envs/snpeff.yaml"
    shell:
        """
        # Decompress VCF for snpEff
        bcftools view {input.vcf} > temp.vcf
        
        # Run SnpEff
        snpEff {params.java_opts} -csvStats {output.csv} \
            -htmlStats {output.html} \
            -dataDir {input.db} \
            {params.genome_version} \
            {params.extra} \
            temp.vcf | \
        bgzip > {output.vcf}
        
        # Index the output VCF
        tabix -p vcf {output.vcf}
        
        # Clean up
        rm temp.vcf
        """

# ANNOVAR annotation
rule annovar:
    input:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.filtered.vcf.gz",
        idx=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.filtered.vcf.gz.tbi" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.filtered.vcf.gz.tbi",
        db=config["annovar_db"]
    output:
        vcf=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.annovar.vcf.gz" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.annovar.vcf.gz",
        idx=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.annovar.vcf.gz.tbi" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.annovar.vcf.gz.tbi",
        txt=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.annovar.txt" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.annovar.txt"
    params:
        prefix=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.merged.annovar" if MERGE_VCFS else f"{OUTDIR}/variants/{wildcards.sample}.{wildcards.caller}.annovar",
        protocol=config.get("annovar_protocol", "refGene,ensGene,knownGene,avsnp150,gnomad211_exome,clinvar_20221231"),
        operation=config.get("annovar_operation", "g,g,g,f,f,f"),
        build=config.get("annovar_build", "hg38")
    log:
        lambda wildcards: f"{OUTDIR}/logs/annovar/{wildcards.sample}.merged.log" if MERGE_VCFS else f"{OUTDIR}/logs/annovar/{wildcards.sample}.{wildcards.caller}.log"
    conda:
        "../envs/annovar.yaml"
    shell:
        """
        # Convert VCF to ANNOVAR format
        bcftools view {input.vcf} > temp.vcf
        convert2annovar.pl -format vcf4 temp.vcf > temp.avinput
        
        # Run ANNOVAR
        table_annovar.pl temp.avinput {input.db} \
            -buildver {params.build} \
            -out {params.prefix} \
            -remove \
            -protocol {params.protocol} \
            -operation {params.operation} \
            -nastring . \
            -vcfinput \
            -polish 2> {log}
            
        # Convert back to VCF and compress
        bgzip -c {params.prefix}.{params.build}_multianno.vcf > {output.vcf}
        tabix -p vcf {output.vcf}
        
        # Clean up
        rm temp.vcf temp.avinput
        """

