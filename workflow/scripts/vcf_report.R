#!/usr/bin/env Rscript

# Load required libraries
library(rmarkdown)
library(knitr)
library(dplyr)

# Capture parameters from Snakemake
sample <- snakemake@params[["sample"]]
callers <- snakemake@params[["callers"]]
outdir <- snakemake@params[["outdir"]]
stats_files <- snakemake@input[["stats"]]
output_report <- snakemake@output[["report"]]
log_file <- snakemake@log[[1]]

# Redirect log output
sink(log_file, append = FALSE, split = TRUE)

# Function to parse BCFTools stats
parse_bcftools_stats <- function(file) {
  # Read the stats file
  stats_content <- readLines(file)
  
  # Extract key metrics
  summary <- list(
    total_variants = NA,
    insertions = NA,
    deletions = NA,
    snps = NA,
    mnps = NA,
    other = NA
  )
  
  # Parse specific lines for different metrics
  for (line in stats_content) {
    if (grepl("^SN", line)) {
      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) >= 3) {
        key <- parts[2]
        value <- as.numeric(parts[3])
        
        if (key == "number of records:") summary$total_variants <- value
        else if (key == "number of insertions:") summary$insertions <- value
        else if (key == "number of deletions:") summary$deletions <- value
        else if (key == "number of SNPs:") summary$snps <- value
        else if (key == "number of MNPs:") summary$mnps <- value
        else if (key == "number of other type:") summary$other <- value
      }
    }
  }
  
  return(summary)
}

# Parse stats for all callers
caller_stats <- lapply(stats_files, parse_bcftools_stats)
names(caller_stats) <- callers

# Prepare data for rendering
report_data <- list(
  sample = sample,
  caller_stats = caller_stats
)

# Render R Markdown report
rmarkdown::render(
  input = system.file("templates", "vcf_qc_report.Rmd", package = "somatic_pipeline"),
  output_file = output_report,
  params = report_data
)

# Close log sink
sink()

