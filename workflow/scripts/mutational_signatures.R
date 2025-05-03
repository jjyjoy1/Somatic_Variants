#!/usr/bin/env Rscript

# Load required libraries
library(BSgenome)
library(MutationalPatterns)
library(ggplot2)
library(tidyr)
library(dplyr)

# Capture Snakemake parameters
sample <- snakemake@params[["sample"]]
outdir <- snakemake@params[["outdir"]]
input_vcf <- snakemake@input[["vcf"]]
ref_genome_path <- snakemake@input[["ref"]]
output_plot <- snakemake@output[["plot"]]
output_table <- snakemake@output[["table"]]
log_file <- snakemake@log[[1]]

# Redirect log output
sink(log_file, append = FALSE, split = TRUE)

# Create output directory if it doesn't exist
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Function to load reference genome
load_reference_genome <- function(ref_path) {
  tryCatch({
    # Check if the reference is a BSgenome object
    if (grepl("\\.rds$", ref_path, ignore.case = TRUE)) {
      # If it's an RDS file with a pre-loaded BSgenome object
      ref_genome <- readRDS(ref_path)
    } else {
      # Attempt to load from a standard reference genome
      ref_name <- basename(ref_path)
      
      # Map common reference genome names to BSgenome packages
      genome_map <- list(
        "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
        "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
        "GRCh38" = "BSgenome.Hsapiens.UCSC.hg38",
        "GRCh37" = "BSgenome.Hsapiens.UCSC.hg19"
      )
      
      # Check if the genome is in our predefined map
      if (ref_name %in% names(genome_map)) {
        library(genome_map[[ref_name]], character.only = TRUE)
        ref_genome <- get(genome_map[[ref_name]])
      } else {
        stop("Unsupported reference genome. Please provide a BSgenome object or use a standard genome.")
      }
    }
    
    return(ref_genome)
  }, error = function(e) {
    warning(paste("Error loading reference genome:", e$message))
    return(NULL)
  })
}

# Load reference genome
ref_genome <- load_reference_genome(ref_genome_path)

# Function to analyze mutational signatures
analyze_mutational_signatures <- function(vcf_file, ref_genome, sample_name) {
  tryCatch({
    # Read VCF file
    vcf <- read_vcf(vcf_file, sample_name, ref_genome)
    
    # Get mutations
    mut_type <- mutation_type(vcf)
    mut_context <- mutation_context(vcf, ref_genome)
    
    # Calculate mutation matrix (96 trinucleotide context)
    mut_matrix <- mutation_matrix(vcf, ref_genome)
    
    # Perform NMF to extract signatures
    nmf_res <- extract_signatures(mut_matrix, rank = 3)
    
    # Plot signatures
    plot_signatures <- plot_96_signatures(nmf_res$signatures)
    
    # Save signature plot
    ggsave(output_plot, plot_signatures, width = 10, height = 6)
    
    # Calculate signature contributions
    sig_contributions <- data.frame(
      Signature = colnames(nmf_res$signatures),
      Contribution = colSums(nmf_res$contribution)
    )
    
    # Write signature contributions to TSV
    write.table(sig_contributions, output_table, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Return results for logging
    list(
      signatures = nmf_res$signatures,
      contributions = sig_contributions
    )
  }, error = function(e) {
    warning(paste("Error analyzing mutational signatures:", e$message))
    return(NULL)
  })
}

# Main execution
tryCatch({
  # Validate inputs
  if (is.null(ref_genome)) {
    stop("Could not load reference genome")
  }
  
  # Analyze mutational signatures
  signatures_result <- analyze_mutational_signatures(input_vcf, ref_genome, sample)
  
  # Log results
  if (!is.null(signatures_result)) {
    cat("Mutational Signatures Analysis Complete\n")
    cat("Signatures Extracted:\n")
    print(signatures_result$signatures)
    cat("\nSignature Contributions:\n")
    print(signatures_result$contributions)
  }
}, error = function(e) {
  # Log any unexpected errors
  cat("Error in mutational signatures analysis:", e$message, "\n")
  sink()
  stop(e)
})

# Close log sink
sink()

