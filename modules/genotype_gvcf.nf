process GENOTYPE_GVCF {
    // One process instance per chromosome × sub-region — replaces the SLURM fan-out
    tag { "chr${chr}_${region}" }

    input:
    // database: the GenomicsDB workspace directory for this chromosome
    // region:   a sub-interval string, e.g. "Pf3D7_13_v3:1-100000"
    tuple val(chr), path(database), val(region), path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    // region label is preserved so GatherVcfs can sort parts in the correct order
    tuple val(chr), val(region), path("chr${chr}_${region.replaceAll(':','_').replaceAll('-','_')}.raw.vcf.gz"), emit: vcf
    tuple val(chr), val(region), path("chr${chr}_${region.replaceAll(':','_').replaceAll('-','_')}.raw.vcf.gz.tbi"), emit: vcf_idx

    script:
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    def safe   = region.replaceAll(':','_').replaceAll('-','_')
    """
    gatk --java-options "-Xmx${mem_gb}G" GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -V gendb://${database} \\
        --max-genotype-count 1024 \\
        -stand-call-conf 30 \\
        -L ${region} \\
        -O chr${chr}_${safe}.raw.vcf.gz

    tabix -p vcf -f chr${chr}_${safe}.raw.vcf.gz
    """
}
