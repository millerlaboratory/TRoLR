# TRoLR - Tandem Repeat outliers (identified with) Long Reads

TRoLR is a comprehensive pipeline for identifying and analyzing tandem repeat outliers in long-read sequencing data, with a focus on detecting potential pathogenic expansions.


## Overview

TRoLR processes haplotype-aware BAM files to:
- Identify tandem repeat variations using VAMOS
- Detect outlier expansions compared to control populations
- Annotate variants with gene/exon/intron/UTR information
- Identify potentially pathogenic STR expansions
- Generate comprehensive HTML reports with visualization

## Installation

### Prerequisites

-Linux/Unix environment
-Conda or Mamba package manager
-R (>= 4.0)
-Python (>= 3.8)
### Required Software

The following tools need to be installed:
- `vamos` for tandem repeat calling
- `bedtools`
 (>= 2.31.1)
- `samtools`
 (>= 1.22)
- `bgzip`
 (part of htslib)

### Quick Install

```
bash

# Clone the repository

git clone https://github.com/yourusername/TRoLR.git
cd TRoLR


# Create conda environment

conda env create -f environment.yml
conda activate trolr-env

# Install R packages

Rscript scripts/install_R_packages.R

```

## Usage

### Basic Usage

```
bash



```

###
 Parameters

- `BAM_HP1`
: Path to first haplotype BAM file (hg38 aligned)
- `BAM_HP2`
: Path to second haplotype BAM file (hg38 aligned)
- `KARYOTYPE`
: Sample karyotype (XX or XY)
- `OUTPUT_DIR`
: Output directory (optional, defaults to current directory)
###
 Example

```
bash

# Run analysis for an XX sample

./TRoLR.sh sample_hp1.hg38.bam sample_hp2.hg38.bam XX ./results/

# Run analysis for an XY sample

./TRoLR.sh sample_hap1.hg38.bam sample_hap2.hg38.bam XY ./results/

```

## Output

The pipeline generates a sample-specific directory containing:
```

<sample_name>/
├── <sample>_hp1.vcf.gz              # VAMOS calls for haplotype 1
├── <sample>_hp2.vcf.gz              # VAMOS calls for haplotype 2
├── <sample>_hp1.strchive.vcf        # STRchive-specific calls
├── <sample>_hp2.strchive.vcf        # STRchive-specific calls
├── <sample>_lps.bed                 # Longest per site (LPS) repeats
├── <sample>_lps_annotated.bed       # Gene-annotated LPS
├── <sample>_lps_annotated_outliers.bed  # Statistical outliers
├── <sample>_pathogenic_results.tsv  # Pathogenic STR calls
├── <sample>_outlier_report.html     # Interactive HTML report
└── plots/                            # Distribution plots
    └── *.png

```

## Reference Data

