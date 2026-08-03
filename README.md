#Human Whole Exome Sequencing Variant Calling using GATK Best Practices

##Project

Germline Variant Discovery from Human Exome Sequencing Data using GATK

##Overview

This project demonstrates a complete germline variant calling workflow using Human Whole Exome Sequencing (WES) data and GATK Best Practices.

This includes quality control, alignment to the human reference genome (GRCh38), duplicate removal, base quality score recalibration (BQSR), variant calling, variant filtering, and generation of high-confidence SNPs and INDELs.

##Objective

To identify high-confidence germline variants from human exome sequencing data using a reproducible bioinformatics pipeline.

##Tools

Linux
FastQC
fastp
BWA-MEM
SAMtools
bcftools
gatk4

##Dataset

reference:
GRCh38

Known Variant Sites:
Homo_sapiens_assembly38.dbsnp138.vcf
Homo_sapiens_assembly38.known_indels.vcf.gz

##sample

BioProject:
PRJNA59853

Run Used:
SRR098356

Platform:
Illumina HiSeq 2000

Layout:
Paired-End

##Workflow

FASTQ
↓
Quality Control (FastQC)
↓
Trimming (Fastp)
↓
Quality Control (FastQC)
↓
Alignment (BWA-MEM)
↓
SAM to BAM
↓
BAM sorting
↓
Mark Duplicates (gatk MarkDuplicates)
↓ 
BQSR (gatk ApplyBQSR)
↓
Variant Calling (gatk HaplotypeCaller)
↓
SNP and INDEL Separation (gatk selectVariants)
↓
Varient Filtering (gatk VariantFiltration)
↓
PASS Variant Selection

## Variant Filtering Criteria

### SNP Filters

- QD < 2.0
- MQ < 40.0
- FS > 60.0
- SOR > 3.0

### INDEL Filters

- QD < 2.0
- FS > 200.0
- SOR > 10.0

##Results

Raw Variants identified: 18201

SNPs: 17112
INDELs: 1089

Afterapplying GATK hard-filtering criteria:

PASS SNPs: 16176
PASS INDELS: 1088

Total High-Confidence Variants: 17264

Total result 94.85% varients retention rate.

##Author

Sheetal






























