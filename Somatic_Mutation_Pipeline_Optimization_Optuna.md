# Somatic Mutation Pipeline Optimization with Optuna

This repository contains a framework for optimizing somatic mutation calling pipelines using Optuna. The optimization focuses on a complete workflow from FastQC to variant calling and filtering, with chromosome 21 subsampling for efficient parameter tuning.

## Overview

The pipeline optimization framework uses Optuna to find optimal parameters for each step in a somatic mutation calling pipeline. It implements a weighted composite score to balance variant detection accuracy, biological plausibility, and computational efficiency.

### Pipeline Steps

1. **FastQ Pre-processing** (fastp)
2. **Alignment** (BWA-mem2)
3. **Duplicate Marking** (GATK MarkDuplicatesSpark)
4. **Base Quality Score Recalibration** (GATK BQSR)
5. **Variant Calling** (GATK Mutect2)
6. **Variant Filtering** (GATK FilterMutectCalls)
7. **Variant Normalization** (BCFtools)

### Key Features

- **Chromosome 21 Subsampling**: Accelerates optimization by ~20-30x while maintaining representative genomic contexts
- **Tree-structured Parzen Estimator (TPE)**: Advanced Bayesian optimization algorithm
- **Weighted Composite Score**: Multi-objective evaluation balancing accuracy, calibration, and resource usage
- **Parameter Importance Analysis**: Identifies the most influential parameters
- **Validation Framework**: Tests optimized parameters on different chromosomes

## Tree-structured Parzen Estimator (TPE)

### What is TPE?

The Tree-structured Parzen Estimator is a sophisticated Bayesian optimization algorithm implemented in Optuna for hyperparameter optimization. Unlike random or grid search methods, TPE uses the history of previous trials to make informed decisions about which parameter values to try next.

### How TPE Works

1. **Probability Modeling**: TPE models the conditional probability p(y|x) of the objective value y given hyperparameters x.

2. **Two-stage Approach**:
   - It constructs two probability distributions: l(x) for promising parameter values (those that gave good results) and g(x) for unpromising values.
   - The algorithm selects parameters that maximize the expected improvement l(x)/g(x).

3. **Tree-structured Search**:
   - TPE handles conditional hyperparameters efficiently through a tree-structured approach.
   - It can discover complex parameter interactions that simpler methods might miss.

4. **Adaptive Sampling**:
   - As more trials are conducted, TPE becomes increasingly focused on high-performing regions of the parameter space.
   - This leads to more efficient exploration compared to grid search or random search.

### TPE Advantages for Bioinformatics Pipelines

- **Efficiency**: TPE requires fewer trials to find optimal parameters compared to grid or random search.
- **High-dimensional Space**: Effectively handles the 20+ parameters in our pipeline.
- **Non-uniform Importance**: Recognizes that some parameters have greater impact than others.
- **Parameter Interactions**: Detects interactions between parameters (e.g., how BWA seed length interacts with Mutect2 parameters).

## Weighted Composite Score (WCS)

We developed a custom weighted composite score to balance multiple competing objectives in somatic variant calling. The score is calculated as:

```
WCS = (0.55 × Accuracy) + (0.30 × Calibration) - (0.15 × ResourceCost)
```

### 1. Accuracy Component (55% weight)

```
Accuracy = (0.6 × SNV_F1) + (0.4 × Indel_F1)
```

- **SNV_F1**: F1 score for single nucleotide variants
- **Indel_F1**: F1 score for insertions/deletions

SNVs receive slightly higher weight (0.6) than indels (0.4) because they're typically more numerous, but both are critical for comprehensive variant detection.

### 2. Calibration Component (30% weight)

```
Calibration = (0.3 × TiTv_Score) + (0.5 × VAF_Score) + (0.2 × PASS_Ratio)
```

Where:
- **TiTv_Score** = 1 - min(1, |observed_TiTv - expected_TiTv|/expected_TiTv)
  - Measures how close the transition/transversion ratio is to the expected value (~3.0 for exome, ~2.0 for genome)
  - Deviations suggest systematic errors

- **VAF_Score** = 1 - KS_statistic
  - Compares the variant allele frequency distribution to ground truth using the Kolmogorov-Smirnov statistic
  - Higher values indicate more accurate representation of variant frequencies

- **PASS_Ratio** = Number of PASS variants / Total variants called
  - Measures the proportion of high-quality variants
  - Higher values suggest more effective filtering

The VAF component receives the highest weight (0.5) because accurate allele frequencies are crucial for somatic variant interpretation, especially in cancer genomics.

### 3. Resource Cost Component (15% weight)

```
ResourceCost = (0.7 × Time_Score) + (0.3 × Memory_Score)
```

Where:
- **Time_Score** = runtime / max_acceptable_runtime
- **Memory_Score** = peak_memory / max_acceptable_memory

Runtime is weighted higher (0.7) than memory usage (0.3) because it typically has more practical impact on pipeline usability, especially in high-throughput environments.

### Rationale for Component Weights

- **Accuracy (55%)**: Given highest weight because detecting true variants is the primary goal of the pipeline
- **Calibration (30%)**: Second highest weight because biological plausibility ensures that variants are legitimate
- **ResourceCost (15%)**: Lower weight but still significant to differentiate between similarly performing parameter sets

### VAF Filtering

The framework filters out variants with VAF < 1% before evaluation, as variants below this threshold are generally not reliable with 100X coverage. This focuses optimization on variants that can be confidently called.

## Usage

### Prerequisites

- Python 3.7+
- GATK 4.2+
- BWA-mem2
- SAMtools/BCFtools
- Optuna
- Required Python packages: numpy, pandas, scipy, scikit-learn

### Installation

```bash
git clone https://github.com/yourusername/somatic-pipeline-optimization.git
cd somatic-pipeline-optimization
pip install -r requirements.txt
```

### Configuration

Edit the paths at the top of `optimize_pipeline.py`:

```python
# Path configurations
REFERENCE_GENOME = "/path/to/reference/genome.fasta"
TRUTH_VCF = "/path/to/NA12878.vcf.gz"
FASTQ_R1 = "/path/to/sample_R1.fastq.gz"
FASTQ_R2 = "/path/to/sample_R2.fastq.gz"
OUTPUT_DIR = "/path/to/output_directory"
KNOWN_SITES = "/path/to/known_sites.vcf.gz"
```

### Running Optimization

```bash
python optimize_pipeline.py
```

This will:
1. Create a chromosome 21 subset of your data
2. Run the optimization with 50 trials
3. Save the best parameters to `best_parameters.txt`
4. Generate visualization plots

### Validating Results

After optimization, validate the best parameters on additional chromosomes:

```bash
python validate_parameters.py --params best_parameters.txt --output validation_results
```

## Advanced Configuration

### Modifying the Parameter Space

To modify which parameters are optimized and their ranges, edit the `objective` function in `optimize_pipeline.py`:

```python
def objective(trial):
    params = {
        # FastQ pre-processing parameters
        "fastp_cut_mean_quality": trial.suggest_int("fastp_cut_mean_quality", 15, 30),
        # ...other parameters
    }
    # ...
```

### Adjusting Composite Score Weights

To adjust the weights in the composite score for specific applications, modify the `calculate_composite_score` function:

```python
def calculate_composite_score(metrics, expected_titv=EXPECTED_TITV):
    # 1. Accuracy component
    accuracy = (0.6 * snv_f1) + (0.4 * indel_f1)
    
    # 2. Calibration component
    calibration = (0.3 * titv_score) + (0.5 * vaf_score) + (0.2 * pass_ratio)
    
    # 3. Resource cost component
    resource_cost = (0.7 * time_score) + (0.3 * memory_score)
    
    # Final weighted score
    wcs = (0.55 * accuracy) + (0.30 * calibration) - (0.15 * resource_cost)
    
    return wcs
```

## Visualization and Analysis

The optimization process generates several visualization plots:

1. **Parameter Importance**: Shows which parameters have the greatest impact on performance
2. **Optimization History**: Shows how the score improves over trials
3. **Parallel Coordinate**: Visualizes parameter combinations for top trials
4. **Parameter Slices**: Shows how individual parameters affect performance

## Extending the Framework

### Adding New Pipeline Steps

To add a new step to the pipeline:

1. Add parameters to the `objective` function
2. Add the execution code to the `run_pipeline` function
3. If needed, add new metrics to the `evaluate_variants` function

### Customizing Evaluation Metrics

To add custom evaluation metrics:

1. Implement the metric calculation in the `evaluate_variants` function
2. Add the new metric to the composite score calculation

## Best Practices

1. **Start with a Small Number of Trials**: Run 5-10 trials first to verify everything works
2. **Monitor Early Results**: Check if parameter ranges need adjustment
3. **Resume Interrupted Runs**: Use the `load_if_exists=True` option with database storage
4. **Analyze Parameter Importance**: Focus refinement on the most impactful parameters
5. **Cross-validate Results**: Test optimized parameters on different samples
6. **Scale Parameter Ranges**: Narrow ranges in follow-up optimization for fine-tuning

## References

1. Bergstra, J., Yamins, D., Cox, D. D. (2013). Making a Science of Model Search: Hyperparameter Optimization in Hundreds of Dimensions for Vision Architectures. ICML.
2. Akiba, T., Sano, S., Yanase, T., Ohta, T., & Koyama, M. (2019). Optuna: A Next-generation Hyperparameter Optimization Framework. KDD.
3. McKenna, A., et al. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research.
4. Firtina, C., et al. (2023). Accelerating BWA-MEM2 alignment on GPUs. Bioinformatics.