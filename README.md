# BAM to Variant Calling Nextflow Pipeline

A Nextflow pipeline for processing BAM files through variant calling with VQSR filtering, based on GATK best practices.

## Pipeline Structure

```
main.nf                           # Main workflow
modules/
├── bam_prep.nf                   # CleanSam → SortSam → MarkDuplicates → region filter
├── bam_stats.nf                  # DepthOfCoverage, InsertSizeMetrics, samtools stats
├── haplotype_caller.nf           # HaplotypeCaller → gVCF (per sample, per chromosome)
├── genomicsdb_import.nf          # GenomicsDBImport (per chromosome)
├── genotype_gvcf.nf              # GenotypeGVCFs (per chromosome, per region)
├── gather_vcfs.nf                # GatherVcfs (per chromosome)
├── variant_recalibrator.nf       # VariantRecalibrator (INDEL + SNP)
└── apply_vqsr.nf                 # ApplyVQSR (INDEL → SNP chain)
nextflow.config                   # Resource configuration
```

## Workflow

```
BAM → BAM_PREP → CleanSam → SortSam → MarkDuplicates → [BED filter]
                    ↓
              BAM_STATS (parallel)
                    ↓
         HAPLOTYPE_CALLER (per sample × chr)
                    ↓
         GENOMICSDB_IMPORT (per chr)
                    ↓
         GENOTYPE_GVCF (per chr × region)
                    ↓
              GATHER_VCFS
                    ↓
         VARIANT_RECALIBRATOR (INDEL)
                    ↓
              APPLY_VQSR (INDEL)
                    ↓
         VARIANT_RECALIBRATOR (SNP)
                    ↓
              APPLY_VQSR (SNP) → Final VCFs
```

## Usage

### Basic run:
```bash
nextflow run main.nf \
    --bam "*.bam" \
    --ref "Pf3D7.fasta" \
    --resource_vcf "Strains.vcf.gz"
```

### With BED region filter:
```bash
nextflow run main.nf \
    --bam "*.bam" \
    --ref "Pf3D7.fasta" \
    --resource_vcf "Strains.vcf.gz" \
    --bed "Pf3D7_core.bed"
```

### With custom parameters:
```bash
nextflow run main.nf \
    --bam "*.bam" \
    --ref "Pf3D7.fasta" \
    --resource_vcf "Strains.vcf.gz" \
    --ploidy 6 \
    --chromosomes "1,2,3,4,5,6,7,8,9,10,11,12,13,14" \
    --hc_threads 16 \
    --db_import_threads 24
```

### Resume:
```bash
nextflow run main.nf -resume --bam "*.bam" --ref "Pf3D7.fasta" --resource_vcf "Strains.vcf.gz"
```

## Required Input Files

| File | Description |
|------|-------------|
| `*.bam` | Input BAM files |
| `reference.fasta` | Reference genome FASTA (with .fai index) |
| `Strains.vcf.gz` | Resource/training VCF for VQSR |
| `core_chr*.list` | Chromosome interval lists for HaplotypeCaller |
| `Genomic_region_list.tsv` | Genomic regions for parallel genotyping |
| `Pf3D7_core.bed` | (Optional) BED file for region filtering |

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bam` | Input BAM file pattern | `*.bam` |
| `--ref` | Reference FASTA file | Required |
| `--resource_vcf` | Training VCF for VQSR | Required |
| `--bed` | BED file for region filtering | Optional |
| `--ploidy` | Sample ploidy | `6` |
| `--chromosomes` | Comma-separated chromosome list | `1-14` |
| `--hc_threads` | HaplotypeCaller threads | `16` |
| `--db_import_threads` | GenomicsDBImport threads | `24` |

## Output Files

- `ReadCoverage_final.tsv` - Coverage statistics
- `InsertSize_Final.txt` - Insert size metrics
- `Bam_stats_Final.tsv` - BAM quality metrics
- `Chr*.raw.vcf.gz` - Raw combined VCFs per chromosome
- `Chr*.raw.recal.vcf.gz` - Final VQSR-filtered VCFs

## Requirements

- GATK 4.2.2.0+
- samtools
- bgzip + tabix
- datamash
