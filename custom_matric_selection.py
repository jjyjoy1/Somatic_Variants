
#Recommended Metric Adjustments:

#Ultra-Low VAF Detection: Add specific metrics for variants in 0.1-1% range
#Error Suppression: Add metric for false positive rate in control samples
#Sequence Context Awareness: Add metrics for error-prone contexts (homopolymers, etc.)

def calculate_composite_score_liquid_biopsy(metrics):
    # 1. Standard accuracy but with VAF-stratified F1 scores
    accuracy_high_vaf = (0.6 * metrics["snv_f1_high_vaf"]) + (0.4 * metrics["indel_f1_high_vaf"])
    accuracy_low_vaf = (0.6 * metrics["snv_f1_low_vaf"]) + (0.4 * metrics["indel_f1_low_vaf"])
    
    # Combined accuracy with emphasis on low VAF performance
    accuracy = (0.3 * accuracy_high_vaf) + (0.7 * accuracy_low_vaf)
    
    # 2. Error suppression in control samples
    error_suppression = 1 - metrics.get("control_false_positive_rate", 0.01)
    
    # 3. Standard calibration
    calibration = (0.2 * titv_score) + (0.6 * vaf_score) + (0.2 * pass_ratio)
    
    # 4. Resource cost (standard)
    resource_cost = (0.7 * time_score) + (0.3 * memory_score)
    
    # Final score with liquid biopsy priorities
    wcs = (0.45 * accuracy) + (0.25 * error_suppression) + 
          (0.20 * calibration) - (0.10 * resource_cost)
    return wcs



#To make your optimization framework adaptable to different cancer types and applications, implement a strategy selector:

def select_optimization_strategy(cancer_type=None, application=None, **kwargs):
    """Select appropriate scoring function based on cancer type and application."""
    
    # Default to standard scoring if no specific type
    if not cancer_type and not application:
        return calculate_composite_score
        
    # Cancer-type specific scoring functions
    if cancer_type:
        if cancer_type.lower() in ["lung", "nsclc", "sclc"]:
            return calculate_composite_score_lung_cancer
        elif cancer_type.lower() in ["leukemia", "lymphoma", "aml", "cll"]:
            return calculate_composite_score_hematological
        elif cancer_type.lower() in ["colorectal", "crc"]:
            return calculate_composite_score_colorectal
        elif cancer_type.lower() in ["glioma", "gbm", "brain"]:
            return calculate_composite_score_brain
            
    # Application-specific scoring functions
    if application:
        if application.lower() in ["clinical", "diagnostic", "diagnosis"]:
            return calculate_composite_score_clinical
        elif application.lower() in ["research", "discovery"]:
            return calculate_composite_score_research
        elif application.lower() in ["liquid_biopsy", "ctdna", "plasma"]:
            return calculate_composite_score_liquid_biopsy
            
    # Fall back to standard scoring if specific type not found
    return calculate_composite_score


