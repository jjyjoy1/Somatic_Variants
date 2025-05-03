#!/usr/bin/env python3

import os
import vcf
import logging
import sys

# Capture Snakemake parameters
input_vcf = snakemake.input['vcf']
output_tmb = snakemake.output['tmb']
sample = snakemake.params['sample']
caller = snakemake.params['caller']
genome_size = snakemake.params['genome_size']
log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def calculate_tmb(vcf_path, genome_size, regions_path=None):
    """
    Calculate Tumor Mutational Burden (TMB)
    
    Args:
        vcf_path (str): Path to input VCF file
        genome_size (int): Effective genome size
        regions_path (str, optional): Path to target regions file
    
    Returns:
        float: Tumor Mutational Burden (mutations per Mb)
    """
    try:
        # Open VCF reader
        vcf_reader = vcf.Reader(filename=vcf_path)
        
        # Filter variants if regions are specified
        variant_count = 0
        target_bases = genome_size
        
        if regions_path and os.path.exists(regions_path):
            # TODO: Implement region-based filtering if needed
            # This would involve reading the regions file and checking variant locations
            pass
        
        # Count variants
        for record in vcf_reader:
            # Count only somatic variants (you may need to adjust filtering criteria)
            if record.is_snp or record.is_indel:
                variant_count += 1
        
        # Calculate TMB (mutations per Mb)
        tmb = (variant_count / (target_bases / 1e6)) if target_bases else 0
        
        logging.info(f"TMB Calculation for {sample} ({caller}):")
        logging.info(f"Total Variants: {variant_count}")
        logging.info(f"Genome Size: {genome_size}")
        logging.info(f"Tumor Mutational Burden: {tmb:.2f} mutations/Mb")
        
        return tmb
    
    except Exception as e:
        logging.error(f"Error calculating TMB: {e}")
        raise

# Calculate TMB
try:
    tmb = calculate_tmb(input_vcf, genome_size, 
                        snakemake.input.get('regions', None))
    
    # Write TMB to output file
    with open(output_tmb, 'w') as f:
        f.write(f"{sample}\t{caller}\t{tmb:.2f}")
    
except Exception as e:
    logging.error(f"TMB calculation failed: {e}")
    sys.exit(1)

