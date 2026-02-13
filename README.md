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

### Quick Install

```
bash

# Clone the repository

git clone github.com/millerlaboratory/TRoLR/TRoLR.git
cd TRoLR


# Create conda environment

conda env create -f environment.yml
conda activate trolr-env

chmod +x TRoLR.sh

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
├── <sample>_performance.log         # Performance log
├── plots/                            # Distribution plots
    └── *.png

```

## Reference Data

```
vamos.motif.hg38.v2.1.e0.1.noSTRCHIVE.nohp.bed.gz                   #TR catalog
STRchive-disease-vamos-20260121.bed                                 #STRChive Catalog
TR_GENCODE_v.45_ANNOTATION.bed                                      #Gencode v45 annotation
STRchive-disease-loci.hg38.general.20260121.tsv                     #STRChive information including pathogenic thresholds
vamos_asm_lps_e0.1_247_catalog_control_length_counts.tsv.gz         #The count of each LPS length at each locus. For generating plots
amos_asm_lps_control_summary.tsv.gz                                 #Summary statistics for identifing outliers

```


