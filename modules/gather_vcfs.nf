process GATHER_VCFS {
    tag { "chr${chr}" }

    publishDir "${params.output ?: '.'}/raw_vcfs", mode: 'copy'

    input:
    tuple val(chr), path(vcf_list), path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    tuple val(chr), path("Chr${chr}.raw.vcf.gz"), emit: vcf
    tuple val(chr), path("Chr${chr}.raw.vcf.gz.tbi"), emit: vcf_idx

    script:
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    def vcf_args = vcf_list.collect { "-V $it" }.join(" ")
    """
    gatk --java-options "-Xmx${mem_gb}G" GatherVcfs \\
        -R ${ref_fasta} \\
        $vcf_args \\
        -O Chr${chr}.raw.vcf.gz \\
        --CREATE_INDEX true

    tabix -p vcf -f Chr${chr}.raw.vcf.gz
    """
}
