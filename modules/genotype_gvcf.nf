process GENOTYPE_GVCFS {
    tag { "chr${chr}" }

    publishDir "${params.output ?: '.'}/raw_vcfs", mode: 'copy'

    input:
    tuple val(chr), path(gvcf), path(gvcf_idx), path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    tuple val(chr), path("Chr${chr}.raw.vcf.gz"), emit: vcf
    tuple val(chr), path("Chr${chr}.raw.vcf.gz.tbi"), emit: vcf_idx

    script:
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    """
    # Index gVCF if not already indexed
    if [ ! -f "${gvcf}.tbi" ] && [ ! -f "${gvcf}.idx" ]; then
        gatk IndexFeatureFile -I ${gvcf}
    fi

    gatk --java-options "-Xmx${mem_gb}G" GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -V ${gvcf} \\
        --max-genotype-count 1024 \\
        -O Chr${chr}.raw.vcf.gz \\
        -stand-call-conf 30

    tabix -p vcf -f Chr${chr}.raw.vcf.gz
    """
}
