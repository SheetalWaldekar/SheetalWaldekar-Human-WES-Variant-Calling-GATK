#!/bin/bash

# ==========================================
# Human WES Variant Calling Pipeline (GATK)
# Sample: SRR098356
# Reference: GRCh38
# ==========================================

# -----------------------------
# Step 1: Quality Control
# -----------------------------
fastqc SRR098356_1.fastq.gz SRR098356_2.fastq.gz

# -----------------------------
# Step 2: Trimming
# -----------------------------
fastp -i SRR098356_1.fastq.gz -I SRR098356_2.fastq.gz -o SRR098356_1_trimmed.fastq.gz -O SRR098356_2_trimmed.fastq.gz

# -----------------------------
# Step 3: Quality Control again
# -----------------------------
fastqc SRR098356_1_trimmed.fastq.gz SRR098356_2_trimmed.fastq.gz

# -----------------------------
# Step 4: BWA index
# -----------------------------
bwa index Homo_sapiens_assembly38.fasta

# -----------------------------
# Step 5: BWA Alignment
# -----------------------------
bwa mem \
-R "@RG\tID:SRR098356\tSM:SRR098356\tPL:ILLUMINA" \
Homo_sapiens_assembly38.fasta \
SRR098356_1.fastq.gz \
SRR098356_2.fastq.gz \
> SRR098356.sam

# -----------------------------
# Step 6: SAM -> BAM
# -----------------------------
samtools view -Sb SRR098356.sam > SRR098356.bam

# -----------------------------
# Step 7: Sort BAM
# -----------------------------
samtools sort \
SRR098356.bam \
-o SRR098356.sorted.bam

# -----------------------------
# Step 8: BAM Index
# -----------------------------
samtools index SRR098356.sorted.bam

# -----------------------------
# Step 9: Mark Duplicates
# -----------------------------
gatk MarkDuplicates \
-I SRR098356.sorted.bam \
-O SRR098356.dedup.bam \
-M duplicate_metrics.txt

samtools index SRR098356.dedup.bam

# -----------------------------
# Step 10: Base Recalibration
# -----------------------------
gatk BaseRecalibrator \
-R Homo_sapiens_assembly38.fasta \
-I SRR098356.dedup.bam \
--known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
--known-sites Homo_sapiens_assembly38.known_indels.vcf.gz \
-O recal.table

# -----------------------------
# Step 11: Apply BQSR
# -----------------------------
gatk ApplyBQSR \
-R Homo_sapiens_assembly38.fasta \
-I SRR098356.dedup.bam \
--bqsr-recal-file recal.table \
-O SRR098356.recal.bam

samtools index SRR098356.recal.bam

# -----------------------------
# Step 12: Variant Calling
# -----------------------------
gatk HaplotypeCaller \
-R Homo_sapiens_assembly38.fasta \
-I SRR098356.recal.bam \
-O SRR098356.raw.vcf.gz

# -----------------------------
# Step 13: Separate SNPs
# -----------------------------
gatk SelectVariants \
-R Homo_sapiens_assembly38.fasta \
-V SRR098356.raw.vcf.gz \
--select-type SNP \
-O raw_snps.vcf.gz

# -----------------------------
# Step 14: Separate INDELs
# -----------------------------
gatk SelectVariants \
-R Homo_sapiens_assembly38.fasta \
-V SRR098356.raw.vcf.gz \
--select-type INDEL \
-O raw_indels.vcf.gz

# -----------------------------
# Step 15: Filter SNPs
# -----------------------------
gatk VariantFiltration \
-R Homo_sapiens_assembly38.fasta \
-V raw_snps.vcf.gz \
-O snps.filtered.vcf.gz \
--filter-name "LowQD" --filter-expression "QD < 2.0" \
--filter-name "LowMQ" --filter-expression "MQ < 40.0" \
--filter-name "HighFS" --filter-expression "FS > 60.0" \
--filter-name "HighSOR" --filter-expression "SOR > 3.0"

# -----------------------------
# Step 16: Filter INDELs
# -----------------------------
gatk VariantFiltration \
-R Homo_sapiens_assembly38.fasta \
-V raw_indels.vcf.gz \
-O indels.filtered.vcf.gz \
--filter-name "LowQD" --filter-expression "QD < 2.0" \
--filter-name "HighFS" --filter-expression "FS > 200.0" \
--filter-name "HighSOR" --filter-expression "SOR > 10.0"

# -----------------------------
# Step 17: Extract PASS SNPs
# -----------------------------
bcftools view \
-f PASS \
snps.filtered.vcf.gz \
-Oz \
-o pass_snps.vcf.gz

# -----------------------------
# Step 18: Extract PASS INDELs
# -----------------------------
bcftools view \
-f PASS \
indels.filtered.vcf.gz \
-Oz \
-o pass_indels.vcf.gz

# -----------------------------
# Step 19: Index PASS VCFs
# -----------------------------
tabix -p vcf pass_snps.vcf.gz
tabix -p vcf pass_indels.vcf.gz

# -----------------------------
# Step 20: Merge PASS Variants
# -----------------------------
bcftools concat \
-a \
-Oz \
-o final_variants.vcf.gz \
pass_snps.vcf.gz \
pass_indels.vcf.gz

# -----------------------------
# Step 21: Variant Summary
# -----------------------------
echo "Raw variants:"
bcftools view -H SRR098356.raw.vcf.gz | wc -l

echo "Raw SNPs:"
bcftools view -H raw_snps.vcf.gz | wc -l

echo "Raw INDELs:"
bcftools view -H raw_indels.vcf.gz | wc -l

echo "PASS SNPs:"
bcftools view -H pass_snps.vcf.gz | wc -l

echo "PASS INDELs:"
bcftools view -H pass_indels.vcf.gz | wc -l

echo "Final variants:"
bcftools view -H final_variants.vcf.gz | wc -l
