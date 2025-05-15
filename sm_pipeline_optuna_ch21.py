import optuna
import subprocess
import os
import time
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support
from scipy.stats import ks_2samp

# Path configurations
REFERENCE_GENOME = "/path/to/reference/genome.fasta"
TRUTH_VCF = "/path/to/truth.vcf.gz"
FASTQ_R1 = "/path/to/sample_R1.fastq.gz"
FASTQ_R2 = "/path/to/sample_R2.fastq.gz"
OUTPUT_DIR = "/path/to/output_directory"
KNOWN_SITES = "/path/to/known_sites.vcf.gz"

# Chromosome for subsampling
TARGET_CHROM = "chr21"  # Use "21" if your reference uses just numbers instead of "chr21"

# Computational resource limits
MAX_RUNTIME_SECONDS = 8 * 3600  # 8 hours (reduced since we're using only chr21)
MAX_MEMORY_GB = 16  # Also reduced for subsampled data

# Expected Ti/Tv ratio (adjust based on whether you're analyzing exome or genome)
EXPECTED_TITV = 3.0  # For exome; use 2.0 for genome

def prepare_chromosome_subset():
    """
    Prepare chromosome 21 subset of reference files and extract chr21 reads from FASTQ.
    This only needs to be done once before optimization starts.
    """
    subset_dir = os.path.join(OUTPUT_DIR, "chr21_subset")
    os.makedirs(subset_dir, exist_ok=True)
    
    # 1. Extract chromosome 21 from reference genome
    chr21_ref = os.path.join(subset_dir, "chr21.fasta")
    if not os.path.exists(chr21_ref):
        cmd = f"samtools faidx {REFERENCE_GENOME} {TARGET_CHROM} > {chr21_ref}"
        subprocess.run(cmd, shell=True, check=True)
        # Index the chromosome 21 reference
        subprocess.run(f"samtools faidx {chr21_ref}", shell=True, check=True)
        subprocess.run(f"bwa-mem2 index {chr21_ref}", shell=True, check=True)
    
    # 2. Extract chromosome 21 from known sites VCF
    chr21_known_sites = os.path.join(subset_dir, "chr21_known_sites.vcf.gz")
    if not os.path.exists(chr21_known_sites):
        cmd = f"bcftools view {KNOWN_SITES} --regions {TARGET_CHROM} -Oz -o {chr21_known_sites}"
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(f"bcftools index {chr21_known_sites}", shell=True, check=True)
    
    # 3. Extract chromosome 21 from truth VCF
    chr21_truth_vcf = os.path.join(subset_dir, "chr21_truth.vcf.gz")
    if not os.path.exists(chr21_truth_vcf):
        cmd = f"bcftools view {TRUTH_VCF} --regions {TARGET_CHROM} -Oz -o {chr21_truth_vcf}"
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(f"bcftools index {chr21_truth_vcf}", shell=True, check=True)
    
    # 4. Extract chromosome 21 reads from FASTQ files
    # This is a two-step process:
    # a) First, align the whole FASTQ to get all alignments
    # b) Then extract only reads that map to chromosome 21
    
    # Check if the chromosome-specific FASTQs already exist
    chr21_fastq_r1 = os.path.join(subset_dir, "chr21_R1.fastq.gz")
    chr21_fastq_r2 = os.path.join(subset_dir, "chr21_R2.fastq.gz")
    
    if not os.path.exists(chr21_fastq_r1) or not os.path.exists(chr21_fastq_r2):
        # First, do a quick alignment to find chr21 reads
        temp_bam = os.path.join(subset_dir, "temp_alignment.bam")
        cmd = (f"bwa-mem2 mem -t 8 {REFERENCE_GENOME} {FASTQ_R1} {FASTQ_R2} | "
               f"samtools view -b -o {temp_bam}")
        subprocess.run(cmd, shell=True, check=True)
        
        # Index the BAM file
        subprocess.run(f"samtools index {temp_bam}", shell=True, check=True)
        
        # Extract chr21 reads and their pairs
        chr21_reads_bam = os.path.join(subset_dir, "chr21_reads.bam")
        cmd = (f"samtools view -b {temp_bam} {TARGET_CHROM} > {chr21_reads_bam}")
        subprocess.run(cmd, shell=True, check=True)
        
        # Extract read names for chr21 alignments
        chr21_read_names = os.path.join(subset_dir, "chr21_read_names.txt")
        cmd = f"samtools view {chr21_reads_bam} | cut -f1 | sort | uniq > {chr21_read_names}"
        subprocess.run(cmd, shell=True, check=True)
        
        # Extract read pairs
        cmd = (f"samtools view -b -N {chr21_read_names} {temp_bam} > {subset_dir}/chr21_read_pairs.bam")
        subprocess.run(cmd, shell=True, check=True)
        
        # Convert BAM to FASTQ
        cmd = (f"samtools fastq -1 {chr21_fastq_r1} -2 {chr21_fastq_r2} {subset_dir}/chr21_read_pairs.bam")
        subprocess.run(cmd, shell=True, check=True)
        
        # Clean up temporary files
        os.remove(temp_bam)
        os.remove(chr21_reads_bam)
        os.remove(chr21_read_names)
        os.remove(f"{subset_dir}/chr21_read_pairs.bam")
    
    # Return paths to chr21-specific files
    return {
        "reference": chr21_ref,
        "known_sites": chr21_known_sites,
        "truth_vcf": chr21_truth_vcf,
        "fastq_r1": chr21_fastq_r1,
        "fastq_r2": chr21_fastq_r2
    }

def run_pipeline(params):
    """Run the complete somatic variant calling pipeline with given parameters, using chr21 subset."""
    # Get chr21 subset files
    chr21_files = prepare_chromosome_subset()
    
    start_time = time.time()
    
    # Create output directory for this run
    run_id = f"run_{int(start_time)}"
    run_dir = os.path.join(OUTPUT_DIR, run_id)
    os.makedirs(run_dir, exist_ok=True)
    
    # STEP 1: FastQC preprocessing with fastp
    fastp_cmd = [
        "fastp",
        "-i", chr21_files["fastq_r1"],
        "-I", chr21_files["fastq_r2"],
        "-o", f"{run_dir}/trimmed_R1.fastq.gz",
        "-O", f"{run_dir}/trimmed_R2.fastq.gz",
        "--cut_mean_quality", str(params["fastp_cut_mean_quality"]),
        "--length_required", str(params["fastp_length_required"]),
        "--low_complexity_filter" if params["fastp_low_complexity_filter"] else "",
        "--adapter_sequence", params["fastp_adapter_sequence"],
        "--detect_adapter_for_pe" if params["fastp_detect_adapter_for_pe"] else "",
        "--poly_x_min_len", str(params["fastp_poly_x_min_len"]),
        "-j", f"{run_dir}/fastp.json",
        "-h", f"{run_dir}/fastp.html"
    ]
    subprocess.run(" ".join(filter(None, fastp_cmd)), shell=True, check=True)
    
    # STEP 2: BWA-mem2 alignment
    bwa_cmd = [
        "bwa-mem2 mem",
        "-t", "8",  # Number of threads, adjust as needed
        "-k", str(params["bwa_seed_length"]),
        "-T", str(params["bwa_min_score"]),
        "-w", str(params["bwa_gap_open_window"]),
        "-M" if params["bwa_picard_compat"] else "",
        chr21_files["reference"],  # Using chr21 reference
        f"{run_dir}/trimmed_R1.fastq.gz",
        f"{run_dir}/trimmed_R2.fastq.gz",
        ">", f"{run_dir}/aligned.sam"
    ]
    subprocess.run(" ".join(filter(None, bwa_cmd)), shell=True, check=True)
    
    # Convert SAM to BAM and sort
    sam_to_bam_cmd = f"samtools view -bS {run_dir}/aligned.sam | samtools sort -o {run_dir}/sorted.bam"
    subprocess.run(sam_to_bam_cmd, shell=True, check=True)
    subprocess.run(f"samtools index {run_dir}/sorted.bam", shell=True, check=True)
    
    # STEP 3: GATK MarkDuplicatesSpark
    markdup_cmd = [
        "gatk MarkDuplicatesSpark",
        "-I", f"{run_dir}/sorted.bam",
        "-O", f"{run_dir}/marked_dups.bam",
        "--optical-duplicate-pixel-distance", str(params["gatk_optical_dup_distance"]),
        "--remove-duplicates" if params["gatk_remove_duplicates"] else "",
        "--spark-master local[8]"  # Adjust based on available cores
    ]
    subprocess.run(" ".join(filter(None, markdup_cmd)), shell=True, check=True)
    
    # STEP 4: GATK BQSR
    bqsr_cmd = [
        "gatk BaseRecalibrator",
        "-I", f"{run_dir}/marked_dups.bam",
        "-R", chr21_files["reference"],  # Using chr21 reference
        "--use-original-qualities" if params["bqsr_use_original_qualities"] else "",
        "--known-sites", chr21_files["known_sites"],  # Using chr21 known sites
        "--interval-padding", str(params["bqsr_interval_padding"]),
        "-O", f"{run_dir}/recal_data.table"
    ]
    subprocess.run(" ".join(filter(None, bqsr_cmd)), shell=True, check=True)
    
    apply_bqsr_cmd = [
        "gatk ApplyBQSR",
        "-I", f"{run_dir}/marked_dups.bam",
        "-R", chr21_files["reference"],  # Using chr21 reference
        "--bqsr-recal-file", f"{run_dir}/recal_data.table",
        "-O", f"{run_dir}/recalibrated.bam"
    ]
    subprocess.run(" ".join(filter(None, apply_bqsr_cmd)), shell=True, check=True)
    
    # STEP 5: Mutect2 variant calling
    mutect2_cmd = [
        "gatk Mutect2",
        "-R", chr21_files["reference"],  # Using chr21 reference
        "-I", f"{run_dir}/recalibrated.bam",
        "-O", f"{run_dir}/somatic_variants.vcf.gz",
        "--min-base-quality-score", str(params["mutect2_min_base_quality"]),
        "--min-pruning", str(params["mutect2_min_pruning"]),
        "--tumor-lod", str(params["mutect2_tumor_lod"]),
        "--initial-tumor-lod", str(params["mutect2_tumor_lod"]),  # Same as tumor-lod for simplicity
        "--initial-tumor-purity", str(params["mutect2_initial_tumor_purity"]),
        "--af-of-alleles-not-in-resource", str(params["mutect2_af_not_in_resource"]),
        "--max-alt-alleles-per-site", str(params["mutect2_max_alt_alleles"]),
        "--max-reads-per-alignment-start", str(params["mutect2_max_reads_per_start"])
    ]
    
    # Add normal sample if using matched mode
    if not params["mutect2_unmatched_mode"]:
        mutect2_cmd.extend(["--normal-sample", "NORMAL"])
    
    subprocess.run(" ".join(filter(None, mutect2_cmd)), shell=True, check=True)
    
    # STEP 6: FilterMutectCalls
    filter_cmd = [
        "gatk FilterMutectCalls",
        "-R", chr21_files["reference"],  # Using chr21 reference
        "-V", f"{run_dir}/somatic_variants.vcf.gz",
        "-O", f"{run_dir}/filtered_variants.vcf.gz"
    ]
    subprocess.run(" ".join(filter(None, filter_cmd)), shell=True, check=True)
    
    # STEP 7: BCFtools normalization and filtering
    bcftools_norm_cmd = [
        "bcftools norm",
        "-m", params["bcftools_norm_multiallelic"],
        "-f", chr21_files["reference"],  # Using chr21 reference
        "-o", f"{run_dir}/normalized.vcf.gz",
        "-O", "z",
        f"{run_dir}/filtered_variants.vcf.gz"
    ]
    subprocess.run(" ".join(filter(None, bcftools_norm_cmd)), shell=True, check=True)
    
    bcftools_filter_cmd = [
        "bcftools filter",
        "-i", params["bcftools_filter_expression"],
        "-o", f"{run_dir}/final_variants.vcf.gz",
        "-O", "z",
        f"{run_dir}/normalized.vcf.gz"
    ]
    subprocess.run(" ".join(filter(None, bcftools_filter_cmd)), shell=True, check=True)
    
    # Index final VCF
    subprocess.run(f"bcftools index {run_dir}/final_variants.vcf.gz", shell=True, check=True)
    
    # Record runtime and estimate memory usage
    runtime = time.time() - start_time
    peak_memory_gb = estimate_peak_memory(run_dir)
    
    # Get metrics on final variant calls
    metrics = evaluate_variants(f"{run_dir}/final_variants.vcf.gz", chr21_files["truth_vcf"])
    metrics["runtime"] = runtime
    metrics["peak_memory_gb"] = peak_memory_gb
    
    return metrics

# The rest of the code (evaluate_variants, compare_variants, etc.) remains the same as in the previous implementation
# Just make sure they work with the chr21-subsampled data

def main():
    """Main function to run the optimization."""
    print("Starting optimization for somatic variant calling pipeline on chromosome 21...")
    
    # Create chr21 subset first (this only needs to be done once)
    chr21_files = prepare_chromosome_subset()
    
    # Create study with TPE sampler
    study = optuna.create_study(
        direction="maximize",  # Maximize composite score
        sampler=optuna.samplers.TPESampler(seed=42),
        pruner=optuna.pruners.MedianPruner(n_startup_trials=5),
        study_name="somatic_pipeline_chr21_optimization",
        storage="sqlite:///chr21_optimization.db",  # Save results to a database
        load_if_exists=True  # Resume if the study already exists
    )
    
    # Run optimization
    try:
        study.optimize(objective, n_trials=50)
    except KeyboardInterrupt:
        print("Optimization interrupted by user. Saving current best results...")
    
    # Print best parameters and results
    print("\nBest trial:")
    trial = study.best_trial
    print(f"  Value: {trial.value}")
    print("  Params: ")
    for key, value in trial.params.items():
        print(f"    {key}: {value}")
    
    # Save best parameters to a file
    with open("best_parameters.txt", "w") as f:
        f.write(f"Best score: {trial.value}\n\n")
        f.write("Parameters:\n")
        for key, value in trial.params.items():
            f.write(f"{key}: {value}\n")
    
    # Create and save visualization plots
    try:
        # Parameter importance
        param_importance_plot = optuna.visualization.plot_param_importances(study)
        param_importance_plot.write_image("param_importances.png")
        
        # Optimization history
        history_plot = optuna.visualization.plot_optimization_history(study)
        history_plot.write_image("optimization_history.png")
        
        # Parallel coordinate plot
        parallel_plot = optuna.visualization.plot_parallel_coordinate(study)
        parallel_plot.write_image("parallel_coordinate.png")
        
        # Slice plot for the most important parameters
        slice_plot = optuna.visualization.plot_slice(study)
        slice_plot.write_image("parameter_slices.png")
    except Exception as e:
        print(f"Warning: Could not create visualization plots: {e}")
    
    print("\nOptimization complete! Results saved to best_parameters.txt and visualization plots.")
    print("To validate these parameters on the full genome, run the final validation script.")

if __name__ == "__main__":
    main()



