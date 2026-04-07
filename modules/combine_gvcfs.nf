process COMBINE_GVCFS {
    tag { "chr${chr}" }

    publishDir "${params.output ?: '.'}/combined_gvcfs", mode: 'copy'

    input:
    tuple val(chr), path(gvcf_list), path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    tuple val(chr), path("Chr${chr}.combined.g.vcf.gz"), emit: combined_gvcf
    tuple val(chr), path("Chr${chr}.combined.g.vcf.gz.tbi"), emit: combined_gvcf_idx

    script:
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    def vcf_args = gvcf_list.collect { "-V $it" }.join(" ")
    """
    # Index gVCFs if needed
    for f in ${gvcf_list.join(" ")}; do
        if [ ! -f "\${f}.tbi" ] && [ ! -f "\${f}.idx" ]; then
            gatk IndexFeatureFile -I \$f
        fi
    done

    gatk --java-options "-Xmx${mem_gb}G" CombineGVCFs \\
        -R ${ref_fasta} \\
        ${vcf_args} \\
        -O Chr${chr}.combined.g.vcf

    bgzip -f Chr${chr}.combined.g.vcf
    tabix -p vcf Chr${chr}.combined.g.vcf.gz
    """
}
