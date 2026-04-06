# Workflow Diagram

```mermaid
graph TD
    A[BAM Files] --> B[BAM_PREP]
    B -->|CleanSam| C[SortSam]
    C --> D[MarkDuplicatesSpark]
    D -->|Optional BED filter| E[Final BAM]
    
    E --> F[BAM_STATS]
    E --> G[HAPLOTYPE_CALLER]
    
    F -->|ReadCoverage_final.tsv| Z[Results]
    F -->|InsertSize_Final.txt| Z
    F -->|Bam_stats_Final.tsv| Z
    
    G -->|Per sample × chr| H[gVCF files]
    H --> I[GENOMICSDB_IMPORT]
    I -->|Per chr| J[GenomicsDB]
    J --> K[GENOTYPE_GVCF]
    K -->|Per chr × region| L[Part VCFs]
    L --> M[GATHER_VCFS]
    M -->|Chr*.raw.vcf.gz| Z
    
    M --> N[VARIANT_RECALIBRATOR INDEL]
    N --> O[APPLY_VQSR INDEL]
    O --> P[VARIANT_RECALIBRATOR SNP]
    P --> Q[APPLY_VQSR SNP]
    Q -->|Chr*.raw.recal.vcf.gz| Z
    
    style A fill:#e1f5fe
    style B fill:#fff3e0
    style G fill:#fff3e0
    style K fill:#f3e5f5
    style Q fill:#e8f5e9
    style Z fill:#c8e6c9
```

## Pipeline Steps

1. **BAM_PREP**: CleanSam → SortSam → MarkDuplicatesSpark → [BED region filter]
2. **BAM_STATS**: DepthOfCoverage, InsertSizeMetrics, samtools stats
3. **HAPLOTYPE_CALLER**: Per-sample per-chromosome gVCF generation
4. **GENOMICSDB_IMPORT**: Combine gVCFs into GenomicsDB per chromosome
5. **GENOTYPE_GVCF**: Joint genotyping per chromosome per region
6. **GATHER_VCFS**: Merge part VCFs into per-chromosome raw VCFs
7. **VARIANT_RECALIBRATOR + APPLY_VQSR**: VQSR filtering (INDEL → SNP chain)
