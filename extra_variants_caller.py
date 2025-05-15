def run_variant_callers(params, run_dir, chr21_files):
    """Run multiple variant callers and merge results."""
    # Create a directory for variant caller results
    variant_callers_dir = os.path.join(run_dir, "variant_callers")
    os.makedirs(variant_callers_dir, exist_ok=True)
    
    # Input BAM file from previous step
    input_bam = f"{run_dir}/recalibrated.bam"
    
    # 1. Run VarDict
    vardict_vcf = run_vardict(params, input_bam, chr21_files, variant_callers_dir)
    
    # 2. Run FreeBayes
    freebayes_vcf = run_freebayes(params, input_bam, chr21_files, variant_callers_dir)
    
    # 3. Run LoFreq
    lofreq_vcf = run_lofreq(params, input_bam, chr21_files, variant_callers_dir)
    
    # 4. Run Strelka2
    strelka2_vcf = run_strelka2(params, input_bam, chr21_files, variant_callers_dir)
    
    # Merge the VCF files from all variant callers
    merged_vcf = os.path.join(variant_callers_dir, "merged_variants.vcf.gz")
    merge_vcfs([vardict_vcf, freebayes_vcf, lofreq_vcf, strelka2_vcf], merged_vcf)
    
    return merged_vcf

def run_vardict(params, input_bam, chr21_files, output_dir):
    """Run VarDict variant caller with optimized parameters."""
    output_vcf = os.path.join(output_dir, "vardict_variants.vcf")
    output_vcf_gz = f"{output_vcf}.gz"
    
    # Create BED file for chromosome 21
    bed_file = os.path.join(output_dir, "chr21.bed")
    cmd = f"samtools view -H {input_bam} | grep SQ | grep '{TARGET_CHROM}' | "
    cmd += f"awk '{{print $2}}' | sed 's/SN://' | "
    cmd += f"awk '{{print $1\"\\t0\\t\"}}' | "
    cmd += f"xargs -I {{}} samtools faidx {chr21_files['reference']} {{}} | "
    cmd += f"awk 'NR>1 {{print \"{TARGET_CHROM}\\t0\\t\"$1}}' > {bed_file}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Prepare VarDict command
    vardict_cmd = [
        "vardict-java",
        "-G", chr21_files["reference"],
        "-f", str(params["vardict_min_allele_freq"]),
        "-N", "TUMOR",
        "-b", input_bam,
        "-c", "1",  # Chromosome column
        "-S", "2",  # Start column
        "-E", "3",  # End column
        "-z", "1",  # 0-based coordinates
        "-th", "8",  # Number of threads
        "--nosv" if params["vardict_no_structural_variants"] else "",
        "-q", str(params["vardict_phred_score"]),
        "-Q", str(params["vardict_mapping_quality"]),
        "-F", str(params["vardict_min_vaf_alt_reads"]),
        bed_file,
        "|",
        "teststrandbias.R",
        "|",
        "var2vcf_valid.pl",
        "-N", "TUMOR",
        "-E",
        "-f", str(params["vardict_min_allele_freq"]),
        ">", output_vcf
    ]
    
    subprocess.run(" ".join(filter(None, vardict_cmd)), shell=True, check=True)
    
    # Compress and index VCF
    subprocess.run(f"bgzip -c {output_vcf} > {output_vcf_gz}", shell=True, check=True)
    subprocess.run(f"tabix -p vcf {output_vcf_gz}", shell=True, check=True)
    
    return output_vcf_gz

def run_freebayes(params, input_bam, chr21_files, output_dir):
    """Run FreeBayes variant caller with optimized parameters."""
    output_vcf = os.path.join(output_dir, "freebayes_variants.vcf")
    output_vcf_gz = f"{output_vcf}.gz"
    
    # Prepare FreeBayes command
    freebayes_cmd = [
        "freebayes",
        "-f", chr21_files["reference"],
        "--min-alternate-fraction", str(params["freebayes_min_alt_fraction"]),
        "--min-alternate-count", str(params["freebayes_min_alt_count"]),
        "--min-base-quality", str(params["freebayes_min_base_quality"]),
        "--min-mapping-quality", str(params["freebayes_min_mapping_quality"]),
        "--pooled-continuous" if params["freebayes_pooled_continuous"] else "",
        "--genotype-qualities" if params["freebayes_genotype_qualities"] else "",
        "--report-genotype-likelihood-max" if params["freebayes_report_gl_max"] else "",
        "-r", TARGET_CHROM,  # Restrict to chromosome 21
        input_bam,
        ">", output_vcf
    ]
    
    subprocess.run(" ".join(filter(None, freebayes_cmd)), shell=True, check=True)
    
    # Compress and index VCF
    subprocess.run(f"bgzip -c {output_vcf} > {output_vcf_gz}", shell=True, check=True)
    subprocess.run(f"tabix -p vcf {output_vcf_gz}", shell=True, check=True)
    
    return output_vcf_gz

def run_lofreq(params, input_bam, chr21_files, output_dir):
    """Run LoFreq variant caller with optimized parameters."""
    output_vcf = os.path.join(output_dir, "lofreq_variants.vcf")
    output_vcf_gz = f"{output_vcf}.gz"
    
    # Prepare LoFreq command
    lofreq_cmd = [
        "lofreq call-parallel",
        "--pp-threads", "8",
        "-f", chr21_files["reference"],
        "-r", TARGET_CHROM,  # Restrict to chromosome 21
        "-s", "-d", str(params["lofreq_min_coverage"]),
        "-a", str(params["lofreq_min_af"]),
        "-q", str(params["lofreq_min_bq"]),
        "-Q", str(params["lofreq_min_mq"]),
        "--call-indels" if params["lofreq_call_indels"] else "",
        "-o", output_vcf,
        input_bam
    ]
    
    subprocess.run(" ".join(filter(None, lofreq_cmd)), shell=True, check=True)
    
    # Compress and index VCF
    subprocess.run(f"bgzip -c {output_vcf} > {output_vcf_gz}", shell=True, check=True)
    subprocess.run(f"tabix -p vcf {output_vcf_gz}", shell=True, check=True)
    
    return output_vcf_gz

def run_strelka2(params, input_bam, chr21_files, output_dir):
    """Run Strelka2 variant caller with optimized parameters."""
    strelka_run_dir = os.path.join(output_dir, "strelka_workdir")
    output_vcf = os.path.join(output_dir, "strelka2_variants.vcf.gz")
    
    # First, create a Strelka2 configuration
    strelka_config_cmd = [
        "configureStrelkaGermlineWorkflow.py",
        "--bam", input_bam,
        "--referenceFasta", chr21_files["reference"],
        "--targeted" if params["strelka_targeted"] else "",
        "--exome" if params["strelka_exome"] else "",
        "--callRegions", TARGET_CHROM,  # Restrict to chromosome 21
        "--runDir", strelka_run_dir
    ]
    
    subprocess.run(" ".join(filter(None, strelka_config_cmd)), shell=True, check=True)
    
    # Modify the Strelka2 configuration file to include custom parameters
    strelka_config_file = os.path.join(strelka_run_dir, "workflow.py")
    
    # Run Strelka2
    strelka_run_cmd = [
        f"{strelka_run_dir}/runWorkflow.py",
        "-m", "local",
        "-j", "8"  # Number of parallel jobs
    ]
    
    subprocess.run(" ".join(filter(None, strelka_run_cmd)), shell=True, check=True)
    
    # Find and copy the Strelka2 output VCF
    strelka_output_vcf = os.path.join(strelka_run_dir, "results/variants/variants.vcf.gz")
    subprocess.run(f"cp {strelka_output_vcf} {output_vcf}", shell=True, check=True)
    subprocess.run(f"cp {strelka_output_vcf}.tbi {output_vcf}.tbi", shell=True, check=True)
    
    return output_vcf

def merge_vcfs(vcf_files, output_vcf):
    """Merge multiple VCF files from different variant callers."""
    # Create a merged VCF file
    bcftools_merge_cmd = [
        "bcftools merge",
        "--force-samples",
        "-m", "none",  # Take the first record in case of overlaps
        "-O", "z",
        "-o", output_vcf
    ]
    
    # Add all VCF files to the command
    bcftools_merge_cmd.extend(vcf_files)
    
    subprocess.run(" ".join(bcftools_merge_cmd), shell=True, check=True)
    subprocess.run(f"tabix -p vcf {output_vcf}", shell=True, check=True)
    
    return output_vcf


# Then update extra variants caller parameter space in the objective function to include the new variant callers' parameters:

def objective(trial):
    """Optuna objective function."""
    # Define parameter space
    params = {
        # FastQ pre-processing parameters
        "fastp_cut_mean_quality": trial.suggest_int("fastp_cut_mean_quality", 15, 30),
        # ... other preprocessing parameters ...
        
        # VarDict parameters
        "vardict_min_allele_freq": trial.suggest_float("vardict_min_allele_freq", 0.01, 0.05),
        "vardict_no_structural_variants": trial.suggest_categorical("vardict_no_structural_variants", [True, False]),
        "vardict_phred_score": trial.suggest_int("vardict_phred_score", 10, 30),
        "vardict_mapping_quality": trial.suggest_int("vardict_mapping_quality", 10, 30),
        "vardict_min_vaf_alt_reads": trial.suggest_float("vardict_min_vaf_alt_reads", 0.01, 0.1),
        
        # FreeBayes parameters
        "freebayes_min_alt_fraction": trial.suggest_float("freebayes_min_alt_fraction", 0.01, 0.05),
        "freebayes_min_alt_count": trial.suggest_int("freebayes_min_alt_count", 2, 5),
        "freebayes_min_base_quality": trial.suggest_int("freebayes_min_base_quality", 10, 30),
        "freebayes_min_mapping_quality": trial.suggest_int("freebayes_min_mapping_quality", 10, 30),
        "freebayes_pooled_continuous": trial.suggest_categorical("freebayes_pooled_continuous", [True, False]),
        "freebayes_genotype_qualities": trial.suggest_categorical("freebayes_genotype_qualities", [True, False]),
        "freebayes_report_gl_max": trial.suggest_categorical("freebayes_report_gl_max", [True, False]),
        
        # LoFreq parameters
        "lofreq_min_coverage": trial.suggest_int("lofreq_min_coverage", 5, 20),
        "lofreq_min_af": trial.suggest_float("lofreq_min_af", 0.01, 0.05),
        "lofreq_min_bq": trial.suggest_int("lofreq_min_bq", 10, 30),
        "lofreq_min_mq": trial.suggest_int("lofreq_min_mq", 10, 30),
        "lofreq_call_indels": trial.suggest_categorical("lofreq_call_indels", [True, False]),
        
        # Strelka2 parameters
        "strelka_targeted": trial.suggest_categorical("strelka_targeted", [True, False]),
        "strelka_exome": trial.suggest_categorical("strelka_exome", [True, False]),
        
        # BCFtools parameters
        "bcftools_norm_multiallelic": trial.suggest_categorical("bcftools_norm_multiallelic", ["-", "+"]),
        "bcftools_filter_expression": trial.suggest_categorical("bcftools_filter_expression", 
                                                              ["QUAL>20 && DP>10", 
                                                               "QUAL>30 && DP>20", 
                                                               "QUAL>40 && DP>30 && AF>=0.01"])
    }
    

