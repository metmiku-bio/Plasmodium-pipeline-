# BAM to Variant Calling Nextflow Pipeline

A Nextflow DSL2 pipeline for processing BAM files through joint variant calling and VQSR filtering,
based on GATK best practices. Designed for polyploid organisms but works with any ploidy.

---

## Pipeline Structure

```
main.nf                           # Main workflow
modules/
├── make_windows.nf               # Generate per-chromosome interval files (Mode A or B)
├── haplotype_caller.nf           # HaplotypeCaller → gVCF (per sample × chromosome)
├── genomicsdb_import.nf          # GenomicsDBImport (per chromosome)
├── genotype_gvcf.nf              # GenotypeGVCFs (per chromosome × region, parallelised)
├── gather_vcfs.nf                # GatherVcfs (merge region VCFs → one per chromosome)
├── variant_recalibrator.nf       # VariantRecalibrator (INDEL + SNP)
└── apply_vqsr.nf                 # ApplyVQSR (INDEL → SNP chain)
nextflow.config                   # Resource and executor configuration
intervals/                        # (Optional) Pre-made GATK interval list files
```

---

## Workflow

```
BAM files
    │
    ▼
MAKE_WINDOWS  ◄── core_chr<N>.list (Mode A)  OR  ref.fasta.fai (Mode B)
    │
    ├── whole.list  ──────────────────────────────────────────────────────┐
    │                                                                     │
    └── windows.list  ────────────────────────────────────────────────┐  │
                                                                      │  │
HAPLOTYPE_CALLER  (per sample × chromosome) ◄────────── whole.list ──┘  │
    │                                                                     │
    ▼                                                                     │
GENOMICSDB_IMPORT  (per chromosome) ◄──────────────── windows.list ──────┘
    │
    ▼
GENOTYPE_GVCF  (per chromosome × region, one job per line of windows.list)
    │
    ▼
GATHER_VCFS  (merge region VCFs → one raw VCF per chromosome)
    │
    ▼
VARIANT_RECALIBRATOR  (INDEL)
    │
    ▼
APPLY_VQSR  (INDEL)
    │
    ▼
VARIANT_RECALIBRATOR  (SNP)
    │
    ▼
APPLY_VQSR  (SNP)  →  Final VCFs
```

---

## Interval Modes

The pipeline supports two ways to define the genomic intervals used for variant calling.
Both modes are handled by the `MAKE_WINDOWS` process and controlled by `--intervals_dir`.

### Mode A — Curated interval lists (recommended when available)

Set `--intervals_dir` to a directory containing one GATK interval list file per chromosome,
named `core_chr<N>.list`.

Each file contains the pre-defined callable regions for that chromosome — for example,
the euchromatic core excluding telomeres, centromeres, and hypervariable repeats.
These regions are used directly as-is, with no re-tiling.

```
core_chr1.list   (e.g. ~700 curated regions for chr1)
core_chr2.list
...
```

**What each file feeds:**
- `HaplotypeCaller -L`    → the whole `core_chr<N>.list` (call only within curated regions)
- `GenomicsDBImport`      → the whole `core_chr<N>.list` (import only those regions)
- `GenotypeGVCFs`         → one parallel job **per line** of `core_chr<N>.list`

> `--num_genome_chunks` is **ignored** in this mode. The curated regions in the list
> file are already the unit of parallelism for `GenotypeGVCFs`.

### Mode B — Automatic tiling from the reference `.fai` (portable, no pre-made files needed)

Omit `--intervals_dir`. The pipeline reads the reference `.fai` and generates intervals
on the fly using `bedtools makewindows`.

**What is generated per chromosome:**
- `whole.list`   → a single `chrN:1-<length>` interval for `HaplotypeCaller`
- `windows.list` → `--num_genome_chunks` equal-width tiles for `GenomicsDBImport`
                   and one parallel `GenotypeGVCFs` job per tile

> `--num_genome_chunks` controls the number of parallel `GenotypeGVCFs` jobs per
> chromosome in this mode. Default: `100`.

---

## Usage

### Mode A — with curated interval lists

```bash
nextflow run main.nf \
    --bam "path/to/*.bam" \
    --ref reference.fasta \
    --resource_vcf training.vcf.gz \
    --intervals_dir ./intervals \
    --chromosomes "1,2,3,4,5,6,7,8,9,10,11,12,13,14" \
    --ploidy 6
```

### Mode B — automatic tiling from .fai

```bash
nextflow run main.nf \
    --bam "path/to/*.bam" \
    --ref reference.fasta \
    --resource_vcf training.vcf.gz \
    --chromosomes "1,2,3,4,5,6,7,8,9,10,11,12,13,14" \
    --num_genome_chunks 200 \
    --ploidy 6
```

### Resume a previous run

```bash
nextflow run main.nf -resume \
    --bam "path/to/*.bam" \
    --ref reference.fasta \
    --resource_vcf training.vcf.gz \
    --intervals_dir ./intervals
```

### Run with a specific profile

```bash
# Local (default)
nextflow run main.nf -profile server ...

# HPC (SLURM)
nextflow run main.nf -profile hpc ...

# Docker container
nextflow run main.nf -profile docker ...

# Singularity container
nextflow run main.nf -profile singularity ...
```

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bam` | Input BAM file glob pattern | `*.bam` |
| `--ref` | Reference FASTA (must have `.fai` and `.dict`) | **Required** |
| `--resource_vcf` | Training/resource VCF for VQSR | **Required** |
| `--chromosomes` | Comma-separated list of chromosomes to process | `1,2,...,14` |
| `--intervals_dir` | Directory with `core_chr<N>.list` files **(Mode A)** | `""` (unset) |
| `--num_genome_chunks` | Windows per chromosome via `bedtools` **(Mode B only)** | `100` |
| `--ploidy` | Sample ploidy passed to HaplotypeCaller | `6` |
| `--max_gaussians` | Max Gaussians for VQSR | `4` |
| `--db_import_threads` | Reader threads for GenomicsDBImport | `4` |
| `--output` | Output directory | `./results` |

---

## Required Input Files

| File | Description |
|------|-------------|
| `*.bam` + `*.bai` | Input BAM files and their indexes |
| `reference.fasta` | Reference genome FASTA |
| `reference.fasta.fai` | FASTA index (samtools faidx) |
| `reference.dict` | Sequence dictionary (GATK CreateSequenceDictionary) |
| `training.vcf.gz` + `.tbi` | Resource VCF for VQSR training |
| `intervals/core_chr<N>.list` | *(Mode A only)* Per-chromosome GATK interval lists |

---

## Output Files

All outputs are written under `--output` (default: `./results`):

```
results/
├── gvcfs/           # Per-sample per-chromosome gVCFs from HaplotypeCaller
├── genomicsdb/      # GenomicsDB workspaces (intermediate)
├── genotyped/       # Per-region raw VCFs from GenotypeGVCFs
├── raw_vcfs/        # Per-chromosome gathered raw VCFs
├── vqsr/            # Recalibration tables, tranches, R plots
└── final_vcfs/      # Final VQSR-filtered VCFs (SNP + INDEL applied)
```

---

## Requirements

| Tool | Version |
|------|---------|
| Nextflow | 22.10+ |
| GATK | 4.2.2.0+ |
| bedtools | 2.29+ *(Mode B only)* |
| samtools | 1.10+ |
| bgzip / tabix | htslib 1.10+ |

### The main code used to run in the server
nextflow run main.nf --bam "../bam_2/*.bqsr.bam" --ref ~/plasmodium_falciparum/Pfalciparum.genome.fasta --resource_vcf ~/plasmodium_falciparum/3d7_hb3.combined.final.vcf.gz --ploidy 6 --chromosomes "13" --output ../result_pipeline/ -profile server  -resume --intervals_dir ./intervals