#!/usr/bin/env Rscript

# Load required libraries
library(rmarkdown)
library(knitr)
library(dplyr)
library(tidyr)

# Capture Snakemake parameters
samples <- snakemake@params[["samples"]]
callers <- snakemake@params[["callers"]]
annotators <- snakemake@params[["annotators"]]
outdir <- snakemake@params[["outdir"]]
output_report <- snakemake@output[["report"]]
output_csv <- snakemake@output[["csv"]]
sample_reports <- snakemake@input[["sample_reports"]]
multiqc_report <- snakemake@input[["multiqc"]]
log_file <- snakemake@log[[1]]

# Redirect log output
sink(log_file, append = FALSE, split = TRUE)

# Create output directory if it doesn't exist
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Function to aggregate variant information across samples
aggregate_variant_summary <- function(samples, callers, annotators) {
  # Initialize empty dataframe to store variant summary
  variant_summary <- data.frame()
  
  # Iterate through samples to collect variant information
  for (sample in samples) {
    sample_summary <- tryCatch({
      # Read TMB files for each caller
      tmb_files <- file.path(outdir, "..", "qc", "tmb", 
                             paste0(sample, ".", callers, ".tmb.txt"))
      
      tmb_data <- lapply(tmb_files, function(file) {
        tryCatch({
          tmb_info <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
          names(tmb_info) <- c("Sample", "Caller", "TMB")
          return(tmb_info)
        }, error = function(e) {
          warning(paste("Error reading TMB file:", file, "-", e$message))
          return(NULL)
        })
      }) %>% do.call(rbind, .)
      
      # Read VCF stats for each caller
      vcf_files <- file.path(outdir, "..", "variants", 
                             paste0(sample, ".", callers, ".filtered.vcf.gz"))
      
      vcf_stats <- lapply(vcf_files, function(file) {
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
          
          data.frame(
            Sample = sample,
            Caller = basename(file),
            TotalVariants = parse_stat("number of records:"),
            SNPs = parse_stat("number of SNPs:"),
            Indels = parse_stat("number of indels:"),
            Other = parse_stat("number of other type:")
          )
        }, error = function(e) {
          warning(paste("Error processing VCF file:", file, "-", e$message))
          return(NULL)
        })
      }) %>% do.call(rbind, .)
      
      # Merge TMB and VCF stats
      full_summary <- merge(tmb_data, vcf_stats, by = c("Sample", "Caller"))
      
      return(full_summary)
    }, error = function(e) {
      warning(paste("Error aggregating data for sample", sample, "-", e$message))
      return(NULL)
    })
    
    # Bind sample summaries
    if (!is.null(sample_summary)) {
      variant_summary <- rbind(variant_summary, sample_summary)
    }
  }
  
  return(variant_summary)
}

# Aggregate variant summary
tryCatch({
  # Aggregate variant information
  variant_summary <- aggregate_variant_summary(samples, callers, annotators)
  
  # Write summary to CSV
  write.csv(variant_summary, output_csv, row.names = FALSE)
  
  # Prepare report data
  report_data <- list(
    samples = samples,
    callers = callers,
    annotators = annotators,
    sample_reports = sample_reports,
    multiqc_report = multiqc_report,
    variant_summary = variant_summary
  )
  
  # Render project-level R Markdown report
  rmarkdown::render(
    input = system.file("templates", "project_report.Rmd", package = "somatic_pipeline"),
    output_file = output_report,
    params = report_data
  )
  
}, error = function(e) {
  # Log any unexpected errors
  cat("Error in project report generation:", e$message, "\n")
  sink()
  stop(e)
})

# Close log sink
sink()

