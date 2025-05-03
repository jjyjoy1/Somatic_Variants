#!/usr/bin/env Rscript

# Load required libraries
library(rmarkdown)
library(knitr)
library(dplyr)
library(ggplot2)

# Capture Snakemake parameters
sample <- snakemake@params[["sample"]]
callers <- snakemake@params[["callers"]]
annotators <- snakemake@params[["annotators"]]
outdir <- snakemake@params[["outdir"]]
output_report <- snakemake@output[["report"]]
log_file <- snakemake@log[[1]]

# Redirect log output
sink(log_file, append = FALSE, split = TRUE)

# Create output directory if it doesn't exist
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Function to load and process TMB data
load_tmb_data <- function(tmb_files) {
  tmb_data <- lapply(tmb_files, function(file) {
    tryCatch({
      tmb_info <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
      names(tmb_info) <- c("Sample", "Caller", "TMB")
      return(tmb_info)
    }, error = function(e) {
      warning(paste("Error reading TMB file:", file, "-", e$message))
      return(NULL)
    })
  })
  
  # Combine TMB data
  tmb_combined <- do.call(rbind, tmb_data)
  return(tmb_combined)
}

# Function to load VCF information
load_vcf_info <- function(vcf_files) {
  vcf_info <- lapply(vcf_files, function(file) {
    tryCatch({
      # Use bcftools to get basic VCF stats
      stats_cmd <- paste("bcftools stats", file, "| grep ^SN")
      stats_output <- system(stats_cmd, intern = TRUE)
      
      # Parse stats
      parse_stat <- function(pattern) {
        stat_line <- grep(pattern, stats_output, value = TRUE)
        if (length(stat_line) > 0) {
          return(as.numeric(strsplit(stat_line, "\t")[[1]][3]))
        }
        return(0)
      }
      
      list(
        file = file,
        total_variants = parse_stat("number of records:"),
        snps = parse_stat("number of SNPs:"),
        indels = parse_stat("number of indels:"),
        other = parse_stat("number of other type:")
      )
    }, error = function(e) {
      warning(paste("Error processing VCF file:", file, "-", e$message))
      return(NULL)
    })
  })
  
  return(do.call(rbind, lapply(vcf_info, as.data.frame)))
}

# Prepare report data
tryCatch({
  # Load TMB data
  tmb_data <- load_tmb_data(snakemake@input[["tmb"]])
  
  # Load VCF information
  vcf_info <- load_vcf_info(snakemake@input[["vcfs"]])
  
  # Prepare report data
  report_data <- list(
    sample = sample,
    bam_qc = snakemake@input[["bam_qc"]],
    vcf_qc = snakemake@input[["vcf_qc"]],
    tmb_data = tmb_data,
    vcf_info = vcf_info,
    merged_vcf = if(length(snakemake@input[["merged_vcf"]]) > 0) snakemake@input[["merged_vcf"]] else NULL,
    annotations = snakemake@input[["annotations"]]
  )
  
  # Render R Markdown report
  rmarkdown::render(
    input = system.file("templates", "sample_report.Rmd", package = "somatic_pipeline"),
    output_file = output_report,
    params = report_data
  )
  
}, error = function(e) {
  # Log any unexpected errors
  cat("Error in sample report generation:", e$message, "\n")
  sink()
  stop(e)
})

# Close log sink
sink()


