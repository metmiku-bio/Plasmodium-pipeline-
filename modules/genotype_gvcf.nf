process GENOTYPE_GVCF {
    tag { "chr${chr}" }

    publishDir "${params.output ?: '.'}/genotyped", mode: 'copy'

    input:
    tuple path(database), path(ref_fasta), path(ref_fai), path(ref_dict), val(chr)

    output:
    tuple val(chr), val("chr${chr}"), path("chr${chr}.raw.vcf.gz"), emit: vcf
    tuple val(chr), val("chr${chr}"), path("chr${chr}.raw.vcf.gz.tbi"), emit: vcf_idx

    script:
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    """
    # Index combined gVCF if needed
    if [ ! -f "${database}.tbi" ] && [ ! -f "${database}.idx" ]; then
        gatk IndexFeatureFile -I ${database}
    fi

    gatk --java-options "-Xmx${mem_gb}G" GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -V ${database} \\
        --max-genotype-count 1024 \\
        -O chr${chr}.raw.vcf.gz \\
        -stand-call-conf 30

    tabix -p vcf -f chr${chr}.raw.vcf.gz
    """
}
