#!/usr/bin/env python
# Variant filtering rules for the somatic pipeline
# Author: Jiyang Jiang
# Date: April 1, 2025

# Filter MuTect2 calls
rule filter_mutect2:
    input:
        vcf=f"{OUTDIR}/variants/{{sample}}.mutect2.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.mutect2.vcf.gz.tbi",
        stats=f"{OUTDIR}/variants/{{sample}}.mutect2.vcf.gz.stats",
        ref=REFERENCE,
        contamination=f"{OUTDIR}/variants/{{sample}}.contamination.table",
        segments=f"{OUTDIR}/variants/{{sample}}.segments.table",
        f1r2=f"{OUTDIR}/variants/{{sample}}.mutect2.f1r2.tar.gz"
    output:
        vcf=f"{OUTDIR}/variants/{{sample}}.mutect2.filtered.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.mutect2.filtered.vcf.gz.tbi"
    params:
        java_opts="-Xmx8g"
    log:
        f"{OUTDIR}/logs/gatk/filter_mutect2/{{sample}}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        # Learn read orientation model
        gatk --java-options "{params.java_opts}" LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {OUTDIR}/variants/{wildcards.sample}.read-orientation-model.tar.gz
        
        # Filter calls
        gatk --java-options "{params.java_opts}" FilterMutectCalls \
            -V {input.vcf} \
            -R {input.ref} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segments} \
            --ob-priors {OUTDIR}/variants/{wildcards.sample}.read-orientation-model.tar.gz \
            -O {output.vcf} 2> {log}
        """

# Calculate contamination for MuTect2 filtering
rule calculate_contamination:
    input:
        unpack(get_caller_input_bams),
        ref=REFERENCE,
        variants=config.get("germline_resource", "")
    output:
        pileups_tumor=temp(f"{OUTDIR}/variants/{{sample}}.tumor.pileups.table"),
        pileups_normal=temp(f"{OUTDIR}/variants/{{sample}}.normal.pileups.table") if has_normal else [],
        contamination=f"{OUTDIR}/variants/{{sample}}.contamination.table",
        segments=f"{OUTDIR}/variants/{{sample}}.segments.table"
    params:
        java_opts="-Xmx8g",
        variant_arg=lambda wildcards, input: f"-V {input.variants}" if input.variants else "",
        normal_pileups=lambda wildcards: f"{OUTDIR}/variants/{wildcards.sample}.normal.pileups.table" if has_normal(wildcards) else ""
    log:
        f"{OUTDIR}/logs/gatk/contamination/{{sample}}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        # Get pileup summaries for tumor
        gatk --java-options "{params.java_opts}" GetPileupSummaries \
            -I {input.tumor_bam} \
            -R {input.ref} \
            {params.variant_arg} \
            -O {output.pileups_tumor} 2> {log}
        
        # Get pileup summaries for normal if available
        if [ -n "{input.normal_bam}" ]; then
            gatk --java-options "{params.java_opts}" GetPileupSummaries \
                -I {input.normal_bam} \
                -R {input.ref} \
                {params.variant_arg} \
                -O {output.pileups_normal} 2>> {log}
            
            # Calculate contamination with matched normal
            gatk --java-options "{params.java_opts}" CalculateContamination \
                -I {output.pileups_tumor} \
                -matched {output.pileups_normal} \
                -O {output.contamination} \
                -segments {output.segments} 2>> {log}
        else
            # Calculate contamination without matched normal
            gatk --java-options "{params.java_opts}" CalculateContamination \
                -I {output.pileups_tumor} \
                -O {output.contamination} \
                -segments {output.segments} 2>> {log}
        fi
        """

# Generic VCF filtering rule for other callers
rule filter_variants:
    input:
        vcf=f"{OUTDIR}/variants/{{sample}}.{{caller}}.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.{{caller}}.vcf.gz.tbi",
        ref=REFERENCE
    output:
        vcf=f"{OUTDIR}/variants/{{sample}}.{{caller}}.filtered.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.{{caller}}.filtered.vcf.gz.tbi"
    params:
        min_depth=config.get("min_depth", 10),
        min_alt_reads=config.get("min_alt_reads", 3),
        min_af=config.get("min_af", 0.05),
        max_af=config.get("max_af", 1.0),
        min_qual=config.get("min_qual", 20),
        caller="{caller}"
    log:
        f"{OUTDIR}/logs/bcftools/filter/{{sample}}.{{caller}}.log"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        # Set up filter expressions based on caller
        FILTER_EXPR=""
        case {params.caller} in
            vardict)
                FILTER_EXPR="FORMAT/DP >= {params.min_depth} && FORMAT/AD[0:1] >= {params.min_alt_reads} && FORMAT/AF >= {params.min_af} && QUAL >= {params.min_qual}"
                ;;
            freebayes)
                FILTER_EXPR="FORMAT/DP >= {params.min_depth} && FORMAT/AO >= {params.min_alt_reads} && FORMAT/AO/FORMAT/DP >= {params.min_af} && QUAL >= {params.min_qual}"
                ;;
            lofreq)
                FILTER_EXPR="DP >= {params.min_depth} && AF >= {params.min_af} && QUAL >= {params.min_qual}"
                ;;
            strelka2)
                if [ -n "$(tabix -H {input.vcf} | grep FORMAT/AU)" ]; then
                    # Somatic mode
                    FILTER_EXPR="FORMAT/DP >= {params.min_depth} && QUAL >= {params.min_qual}"
                else
                    # Germline mode
                    FILTER_EXPR="FORMAT/DP >= {params.min_depth} && QUAL >= {params.min_qual}"
                fi
                ;;
            varscan2)
                FILTER_EXPR="FORMAT/DP >= {params.min_depth} && FORMAT/AD[0:1] >= {params.min_alt_reads} && FORMAT/FREQ >= {params.min_af} && QUAL >= {params.min_qual}"
                ;;
            deepvariant)
                FILTER_EXPR="FORMAT/DP >= {params.min_depth} && FORMAT/AD[0:1] >= {params.min_alt_reads} && FORMAT/VAF >= {params.min_af} && QUAL >= {params.min_qual}"
                ;;
            *)
                # Default filter
                FILTER_EXPR="FORMAT/DP >= {params.min_depth} && QUAL >= {params.min_qual}"
                ;;
        esac

        # Apply filters
        bcftools filter -i "$FILTER_EXPR" {input.vcf} | \
        bcftools view -f PASS -O z -o {output.vcf}
        
        # Index filtered VCF
        tabix -p vcf {output.vcf}
        """

# Filter merged VCF
rule filter_merged_vcf:
    input:
        vcf=f"{OUTDIR}/variants/{{sample}}.merged.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.merged.vcf.gz.tbi"
    output:
        vcf=f"{OUTDIR}/variants/{{sample}}.merged.filtered.vcf.gz",
        idx=f"{OUTDIR}/variants/{{sample}}.merged.filtered.vcf.gz.tbi"
    params:
        min_depth=config.get("min_depth", 10),
        min_qual=config.get("min_qual", 20)
    log:
        f"{OUTDIR}/logs/bcftools/filter_merged/{{sample}}.log"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        # Apply filters to merged VCF
        bcftools filter -i "FORMAT/DP >= {params.min_depth} && QUAL >= {params.min_qual}" {input.vcf} | \
        bcftools view -f PASS -O z -o {output.vcf}
        
        # Index filtered VCF
        tabix -p vcf {output.vcf}
        """

