![Logo](banner/TRoLR_logo.png)
# TRoLR - Tandem Repeat outliers (identified with) Long Reads

TRoLR is a comprehensive pipeline for identifying and analyzing tandem repeat outliers in long-read sequencing data, with a focus on detecting potential pathogenic expansions.

TRoLR uses VAMOS developed by the Chaisson Lab (https://github.com/ChaissonLab/vamos, Ren et al., 2023)

*TRoLR v.1 has been depricated in favor of TRComp DB v2.0 and vamos_v3.0.6. See releases if you would prefer to install the old version*

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

### Quick Install

```
bash

# Clone the repository

git clone git@github.com:millerlaboratory/TRoLR.git
cd TRoLR


# Create conda environment

conda env create -f environment.yaml
conda activate trolr-env

chmod +x TRoLR.sh

```

### Install Reference Data

The reference data is hosted on AWS as it is too large to host here.

In the TRoLR directory run the following:

```
bash

wget wget https://s3.amazonaws.com/1000g-ont/Gibson_etal_TRoLR_preprint_supplemental/TRoLR_reference_data.tar.gz

tar -xzf TRoLR_reference_data.tar.gz
```

## Usage

### Basic Usage

**TRoLR takes haplotype-resolved assemblies aligned to hg38 as its input, it is recommened to generate these from fastq files using hifiasm:https://github.com/chhylp123/hifiasm**

```
bash
./TRoLR.sh SAMPLE_NAME <BAM_HP1> <BAM_HP2> <KARYOTYPE: XX|XY> [OUTPUT_DIR]

```

###
 Parameters
- `SAMPLE_NAME`
: Sample name for output naming scheme
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

./TRoLR.sh SAMPLE_NAME sample_hp1.hg38.bam sample_hp2.hg38.bam XX ./results/

# Run analysis for an XY sample

./TRoLR.sh SAMPLE_NAME sample_hap1.hg38.bam sample_hap2.hg38.bam XY ./results/

```
## Test Sample

Using GIAB sample data from *chr9:82230983-92919503*, you can test the program by running

```
./TRoLR.sh HG002 HG002_sample_data/HG002_asm_1.hg38.bam HG002_sample_data/HG002_asm_2.hg38.bam XY ./results/
```

**The `results` directory must be an already existing filepath**

## Output

The pipeline generates a sample-specific directory in the output containing:
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
├── <sample>_performance.log         # Performance log
├── plots/                            # Distribution plots
    └── *.png

```
An example html output from the HG002 sample data can be found in `example_html`

## Reference Data

```
DATA_FILE_1_Curated_TRCompV2_eff0.1_noHp_noSTRchive.tsv   #TR catalog
DATA_FILE_2_Curated_STRchive_loci_v1.0_20260330.bed       #STRChive Catalog
TR_GENCODE_v.45_ANNOTATION_PHENOTYPE_SORTED.bed           #Gencode v45 annotation
STRCHIVE_locus_information_annotation.tsv                 #STRChive information including pathogenic thresholds
Supp_LPS_MOTIF_ALLELE_COUNTS.tsv.gz                       #The count of each LPS length at each locus. For generating plots
DATA_FILE_5_ALL_LOCI_SUMMARY_STATS_PER_MOTIF.tsv.gz       #Summary statistics for identifing outliers

```
## Citing TRoLR

*Preprint information coming soon!*

