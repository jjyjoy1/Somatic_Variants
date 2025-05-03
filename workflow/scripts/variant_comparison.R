#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(VennDiagram)
library(dplyr)

# Capture Snakemake parameters
sample <- snakemake@params[["sample"]]
callers <- snakemake@params[["callers"]]
outdir <- snakemake@params[["outdir"]]
vcf_files <- snakemake@input[["vcfs"]]
venn_output <- snakemake@output[["venn"]]
stats_output <- snakemake@output[["table"]]
log_file <- snakemake@log[[1]]

# Redirect log output
sink(log_file, append = FALSE, split = TRUE)

# Create output directory if it doesn't exist
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Function to extract variant positions from VCF files
extract_variant_positions <- function(vcf_file) {
  tryCatch({
    # Use bcftools to extract variant positions
    cmd <- paste("bcftools query -f '%CHROM\\t%POS\\n'", vcf_file)
    variant_positions <- system(cmd, intern = TRUE)
    
    # Convert to data frame
    positions_df <- read.table(text = variant_positions, 
                               col.names = c("Chromosome", "Position"))
    
    # Create unique identifier
    positions_df$Variant_ID <- paste(positions_df$Chromosome, 
                                     positions_df$Position, 
                                     sep = "_")
    
    return(positions_df$Variant_ID)
  }, error = function(e) {
    warning(paste("Error processing VCF file:", vcf_file, "-", e$message))
    return(character(0))
  })
}

# Perform variant comparison
tryCatch({
  # Extract variant positions for each caller
  caller_variants <- lapply(vcf_files, extract_variant_positions)
  names(caller_variants) <- callers
  
  # Create Venn diagram
  venn_plot <- venn.diagram(
    x = caller_variants,
    filename = venn_output,
    fill = rainbow(length(callers)),
    alpha = 0.5,
    cex = 1.5,
    fontfamily = "sans",
    category.names = callers,
    main = paste("Variant Comparison for Sample", sample)
  )
  
  # Calculate variant overlap statistics
  overlap_stats <- data.frame(
    Caller = callers,
    Total_Variants = sapply(caller_variants, length),
    Unique_Variants = sapply(caller_variants, function(x) {
      sum(table(x) == 1)
    })
  )
  
  # Calculate pairwise overlaps
  overlap_matrix <- matrix(0, nrow = length(callers), ncol = length(callers))
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- callers
  
  for (i in 1:(length(callers)-1)) {
    for (j in (i+1):length(callers)) {
      overlap <- length(intersect(caller_variants[[i]], caller_variants[[j]]))
      overlap_matrix[i,j] <- overlap_matrix[j,i] <- overlap
    }
  }
  
  # Prepare detailed overlap statistics
  overlap_stats_detailed <- as.data.frame(as.table(overlap_matrix))
  names(overlap_stats_detailed) <- c("Caller1", "Caller2", "Overlap_Count")
  
  # Write overlap statistics to TSV
  write.table(overlap_stats_detailed, 
              stats_output, 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
  
  # Print summary for logging
  print("Variant Comparison Summary:")
  print(overlap_stats)
  print("Pairwise Overlap Matrix:")
  print(overlap_matrix)
  
}, error = function(e) {
  # Log any unexpected errors
  cat("Error in variant comparison:", e$message, "\n")
  sink()
  stop(e)
})

# Close log sink
sink()


