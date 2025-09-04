# TESFaC

**TesFac** (**T**herapeutic **E**fficacy **S**tudies **Fa**ilure **C**lassification) is a Nextflow pipeline for detecting treatment failure in malaria longitudinal studies using genomic data analysis and machine learning.

## Overview

This pipeline analyzes malaria genomic data to distinguish between treatment failure (recrudescent infections) and new infections. It processes genomic amplicon data, generates synthetic training datasets, calculates genomic features, and applies machine learning models to predict treatment outcomes.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Input Data](#input-data)
- [Parameters](#parameters)
- [Output](#output)
- [Pipeline Steps](#pipeline-steps)
- [Requirements](#requirements)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (≥21.04)
- [Singularity](https://sylabs.io/singularity/) (≥3.6)
- Minimum 32GB RAM recommended
- 20+ CPU cores recommended for optimal performance

### Building the Singularity Container

```bash
# Clone the repository
git clone https://github.com/your-username/TESFaC.git
cd TESFaC

# Build the Singularity container
singularity build TESFaC.sif Singularity.def
```

## Quick Start

```bash
# Run the complete pipeline
nextflow run main.nf \
    --site Tete \
    --raw_data_dir data/raw/ \
    -profile singularity
```

## Pipeline Overview

The TESFaC pipeline consists of 9 sequential processes:

1. **Data Cleaning** - Quality control and preprocessing
2. **COI Calculation** - Complexity of infection estimation
3. **FNR Calculation** - False negative rate calculation
4. **Clone Generation** - Monoclonal strain extraction
5. **Mix Generation** - Synthetic mixture creation
6. **Pairs Generation** - Training pair construction
7. **Feature Calculation** - Genomic feature extraction
8. **Target Augmentation** - Training data enhancement
9. **ML Modeling** - Machine learning classification

## Input Data

### Required Files

Place the following files in your `--raw_data_dir`:

```
data/raw/
├── metadata_tes_<SITE>.csv          # Sample metadata of the failure pairs
├── genomic_site_<SITE>.csv          # Site genomic data (aka contextual samples from the site the failure pairs belong to)
└── *_RESULTS*/                      # TES failure pairs sequencing run directory
    └── allele_data_global_max_0_filtered.csv
```

### Required Resources

```
resources/
└── pfPHAST_20amps.csv               # Amplicon selection file
```

### Metadata Format

The metadata file must contain these columns:
- `NIDA`: Sample identifier
- `PairsID`: Paired sample identifier
- `time_point`: D0 or Dx
- `SampleID`: Additional sample ID

### Genomic Data Format

Genomic files must contain:
- `sampleID`: Sample identifier
- `locus`: Amplicon/locus name
- `pseudo_cigar`: Sequence variant
- `reads`: Read counts
- `norm.reads.locus`: Normalized read frequency

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--site` | Study site name | `Tete` |
| `--raw_data_dir` | Input data directory | `data/raw/` |

### Optional Parameters

#### Data Processing
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--maf_filter` | 0.01 | Minor allele frequency threshold |
| `--min_allele_read_count` | 10 | Minimum reads per allele |
| `--cum_curve_threshold` | 0.95 | Cumulative curve threshold |
| `--optim_ampset` | false | Use optimized amplicon set |

#### Computational Resources
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | 20 | Maximum CPU cores |
| `--max_memory` | 60.GB | Maximum memory |
| `--chunk_size` | 100 | Processing chunk size |

#### Training Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--n_clones_target` | 1000 | Target number of clones |
| `--initial_sample_size` | 200 | Initial mixture sample size |
| `--train_test_split` | 0.7 | Training/test split ratio |
| `--cv_folds` | 10 | Cross-validation folds |

## Output

Results are saved to `results_<SITE>/`:

```
results_<SITE>/
├── processed/                       # Intermediate files
│   ├── genomic_updated_<SITE>.csv
│   ├── metadata_updated_<SITE>.csv
│   ├── clones_genomic_data_<SITE>.csv
│   ├── MIXES_METADATA_<SITE>.RDS
│   └── FNRs_<SITE>.RDS
└── results/                         # Final outputs
    ├── <SITE>_REAL_DATA_PREDICTIONS.csv    # Main predictions
    ├── GLMNET_model_<SITE>.RDS             # Trained model
    ├── training_Results_GLMNET_<SITE>.csv  # Training results
    ├── <SITE>_model_results_GLMNET.png     # Performance plots
    └── feat_importance_<SITE>.png          # Feature importance
```

### Key Output Files

- **Predictions**: `<SITE>_REAL_DATA_PREDICTIONS.csv` contains treatment failure predictions
- **Model**: `GLMNET_model_<SITE>.RDS` is the trained GLMNET model
- **Results**: Training performance metrics and decision thresholds
- **Plots**: Model performance and feature importance visualizations

## Pipeline Steps

### 1. Data Cleaning
- Quality control filtering
- Sample ID formatting
- Amplicon selection
- Data merging

### 2. COI Calculation
- Complexity of infection estimation using offset-naive method (second-highest allele count per locus)

### 3. FNR Calculation
- False negative rate calculation based on strain proportions
- Experimentally derived FNR values

### 4. Clone Generation
- Natural clone identification (1 allele per locus)
- Synthetic clone generation from polyclonal samples
- Redundant clone removal

### 5. Mix Generation
- Rarefaction curve analysis for clone sufficiency
- Clone culling (natural clones are favored over synthetic clones)
- Strain mixture creation across COI levels

### 6. Pairs Generation
- Training pair construction (D0-Dx)
- Recrudescent vs. new infection labeling

### 7. Feature Calculation
- Genomic feature extraction:
  - Jaccard similarity
  - Allele retention rate
  - Allele gain
  - Allele loss
  - Locus discordance rate
  - Transition asymmetry

### 8. Target Augmentation
- FNR-based data augmentation
- Synthetic sample generation for pair types that with low frequency clones in the real data

### 9. ML Modeling
- GLMNET model training with cross-validation
- Optimal threshold determination
- Real data prediction
- Performance evaluation

## Requirements

### System Requirements
- Linux/Unix system
- Minimum 32GB RAM
- 20+ CPU cores (recommended)
- 100GB+ free disk space

## Troubleshooting

### Common Issues

#### Memory Errors
```bash
# Increase memory allocation
nextflow run main.nf --max_memory 128.GB
```

#### Missing Input Files
Check that all required files exist with the correct names:
- `metadata_tes_<SITE>.csv`
- `genomic_site_<SITE>.csv` 
- Sequencing run with the TES failure pairs with the correct naming (see [Input Data](#input-data)).

#### Singularity Issues
```bash
# Build container if needed
singularity build --force TESFaC.sif Singularity.def
```

### Performance Optimization

#### For Large Datasets
```bash
# Increase chunk size and cores
nextflow run main.nf \
    --chunk_size 200 \
    --max_cpus 40 \
    --max_memory 128.GB
```

#### For Limited Resources
```bash
# Reduce targets and sample sizes
nextflow run main.nf \
    --n_clones_target 500 \
    --initial_sample_size 100 \
    --max_cpus 8
```

## Configuration Profiles

### Singularity Profile
```bash
nextflow run main.nf -profile singularity
```

### Custom Configuration
Create a custom config file:

```nextflow
// custom.config
params {
    max_cpus = 40
    max_memory = '128.GB'
    chunk_size = 200
}

process {
    withName: 'FEATURE_CALCULATION' {
        cpus = 40
        memory = '128.GB'
        time = '24.h'
    }
}
```

Run with:
```bash
nextflow run main.nf -c custom.config
```

## Advanced Usage

### Running Specific Processes
```bash
# Run only data processing steps
nextflow run main.nf --site Tete -entry DATA_CLEANING
```

### Resume Failed Runs
```bash
# Resume from last checkpoint
nextflow run main.nf -resume
```

### Multiple Sites
```bash
# Process multiple sites
for site in Tete Maputo Zambezia; do
    nextflow run main.nf --site $site --raw_data_dir data/raw/
done
```

## Citations

If you use this pipeline, please cite the TESFaC repository.

## Support

For issues and questions:
- Open an issue on GitHub
- Check the troubleshooting section
- Review the Nextflow documentation

## License

[Add your license information here]
