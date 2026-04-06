process GATHER_VCFS {
    tag { "chr${chr}" }

    input:
    path vcf_list
    path ref_fasta
    val chr

    output:
    path "Chr${chr}.raw.vcf.gz", emit: vcf
    path "Chr${chr}.raw.vcf.gz.tbi", emit: vcf_idx

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga() - 4}G" GatherVcfs \\
        -R ${ref_fasta} \\
        -I ${vcf_list} \\
        -O Chr${chr}.raw.vcf.gz \\
        --CREATE_INDEX true
    
    tabix -p vcf Chr${chr}.raw.vcf.gz
    """
}
