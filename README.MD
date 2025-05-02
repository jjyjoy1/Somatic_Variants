# Somatic Variant Calling Pipeline

A comprehensive Snakemake pipeline for somatic variant calling that supports both tumor-only and tumor-normal paired samples. This pipeline covers the entire process from raw sequencing data to annotated variants.

## Features

- Supports popular human reference genomes (hg38/GRCh38 and hg19/GRCh37)
- Multiple variant callers with configurable selection:
  - GATK MuTect2
  - VarDict
  - FreeBayes
  - Lofreq
  - Strelka2
  - Samtools/VarScan2
- Multiple annotation tools:
  - GATK Funcotator
  - Ensembl VEP
  - SnpEff
  - ANNOVAR
- Comprehensive QC for FASTQ, BAM, and VCF files
- Automatic VCF merging when multiple callers are used
- Detailed reporting at sample and project levels
- Mutational signature analysis
- Parameter optimization via Optuna (optional)

## Requirements

- Snakemake (≥7.18.0)
- Conda or Mamba
- Python 3.8+
- 16+ GB RAM (recommended)
- 4+ CPU cores (recommended)

## Installation

Clone this repository:

```bash
git clone https://github.com/jjyjoy1/somatic-variant-pipeline.git
cd somatic-variant-pipeline
```

## Usage

### 1. Configure your samples

Edit the `config/config.yaml` file to specify your samples, reference genome, and other parameters.

Example configuration for tumor-normal paired samples:
```yaml
samples:
  tumor1:
    patient: "patient1"
    sample_type: "tumor"
  normal1:
    patient: "patient1"
    sample_type: "normal"
  tumor2:
    patient: "patient2"
    sample_type: "tumor"

units:
  tumor1:
    fq1: "data/tumor1_R1.fastq.gz"
    fq2: "data/tumor1_R2.fastq.gz"
  normal1:
    fq1: "data/normal1_R1.fastq.gz"
    fq2: "data/normal1_R2.fastq.gz"
  tumor2:
    fq1: "data/tumor2_R1.fastq.gz"
    fq2: "data/tumor2_R2.fastq.gz"
```

### 2. Select variant callers and annotators

In the `config/config.yaml` file:

```yaml
# Variant callers to use (comment out to disable)
variant_callers:
  - "mutect2"
  - "vardict"
  - "freebayes"
  - "lofreq"
  - "strelka2"
  - "varscan2"

# Variant annotators to use
variant_annotators:
  - "funcotator"
  - "vep"
  - "snpeff"
  # - "annovar"  # Uncomment to enable
```

### 3. Run the pipeline

```bash
# Dry run (to check the workflow)
snakemake --use-conda -n

# Actual run
snakemake --use-conda --cores 16

# Generate workflow diagram
snakemake --dag | dot -Tsvg > workflow.svg
```

## Directory Structure

```
├── config/
│   └── config.yaml       # Main configuration file
├── workflow/
│   ├── rules/            # Snakemake rules
│   ├── scripts/          # Custom scripts
│   └── envs/             # Conda environments
├── resources/            # Reference genomes, databases
└── results/              # Pipeline outputs
    ├── qc/
    │   ├── fastqc/
    │   ├── bam/
    │   └── vcf/
    ├── mapping/
    ├── variants/
    └── reports/
```

## Pipeline Workflow

1. **Quality Control**: FastQC and MultiQC for raw reads
2. **Preprocessing**: Trimming, alignment, and duplicate marking
3. **BAM Processing**: BQSR, indexing, and quality assessment
4. **Variant Calling**: Multiple callers with sample-specific processing
5. **Variant Filtering**: Caller-specific and general filters
6. **Variant Merging**: When multiple callers are used
7. **Variant Annotation**: Multiple annotators for functional prediction
8. **Reporting**: Comprehensive QC metrics and summary statistics

## Optuna Parameter Optimization

This pipeline supports Optuna for parameter optimization. Below is a list of parameters that can be optimized:

### Alignment Parameters (BWA-MEM2)
- Match score
- Mismatch penalty
- Gap open penalty
- Gap extension penalty
- Seed length

### MuTect2 Parameters
- Tumor LOD threshold
- Normal LOD threshold
- Downsampling threshold
- Minimum pruning threshold

### VarDict Parameters
- Minimum allele frequency
- Minimum variant depth
- Minimum base quality
- Minimum mapping quality

### FreeBayes Parameters
- Minimum alternate fraction
- Minimum mapping quality
- Minimum base quality
- Minimum alternate count

### Lofreq Parameters
- Minimum coverage
- Minimum base quality
- Significance threshold
- Minimum mapping quality

### Strelka2 Parameters
- Maximum indel size
- Minimum trailing base quality
- Indel candidate window size

### VarScan2 Parameters
- Minimum coverage
- Minimum variant allele frequency
- P-value threshold
- Minimum supporting reads

### DeepVariant Parameters
- Model type selection
- Number of shards
- Batch size for inference

### Filtering Parameters
- Min/max depth thresholds
- Minimum quality score
- Strand bias threshold
- Mapping quality bias threshold

To use Optuna with this pipeline, set `optuna.enabled: true` in the config file and specify which parameters to optimize. For detailed Optuna integration, refer to the companion repository [somatic-pipeline-optuna](https://github.com/username/somatic-pipeline-optuna).

## Output Files

### QC Reports
- `results/qc/multiqc_report.html`: Combined QC report
- `results/qc/bam/{sample}/qualimap_report.html`: BAM quality metrics
- `results/qc/vcf/{sample}.report.html`: VCF quality metrics

### Variants
- `results/variants/{sample}.{caller}.filtered.vcf.gz`: Filtered variants per caller
- `results/variants/{sample}.merged.vcf.gz`: Merged variants (if multiple callers)
- `results/variants/{sample}.{caller}.{annotator}.vcf.gz`: Annotated variants

### Reports
- `results/reports/{sample}/report.html`: Sample-level report
- `results/reports/project_report.html`: Project-level summary
- `results/reports/{sample}/mutational_signatures.pdf`: Mutational signature analysis
- `results/reports/{sample}/caller_comparison.pdf`: Caller comparison visualization

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

jiyang Jiang email: jiyang.jiang@gmail.com
For questions or support, please open an issue on the GitHub repository.
