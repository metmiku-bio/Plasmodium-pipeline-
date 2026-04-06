process COMBINE_GVCFS {
    tag { "chr${chr}" }

    input:
    tuple val(chr), path(gvcf_files), path(ref_fasta), path(ref_dict)

    output:
    tuple val(chr), path("*.combined.g.vcf.gz"), emit: combined_gvcf
    tuple val(chr), path("*.combined.g.vcf.gz.tbi"), emit: combined_gvcf_idx

    script:
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    def gvcf_list = gvcf_files instanceof List ? gvcf_files : [gvcf_files]
    def gvcf_args = gvcf_list.collect { "-V $it" }.join(" ")
    """
    gatk --java-options "-Xmx${mem_gb}G" CombineGVCFs \\
        -R ${ref_fasta} \\
        ${gvcf_args} \\
        -O chr${chr}.combined.g.vcf
    """
}
