process APPLY_VQSR {
    tag { "chr${chr}.${mode}" }

    publishDir "${params.output ?: '.'}/final_vcfs", mode: 'copy'

    input:
    tuple path(vcf), path(recal), path(tranches), val(chr), val(mode)

    output:
    path "*.recal.vcf.gz", emit: vcf
    path "*.recal.vcf.gz.tbi", emit: vcf_idx

    script:
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    """
    gatk --java-options "-Xmx${mem_gb}g" ApplyVQSR \\
        -V ${vcf} \\
        --recal-file ${recal} \\
        --tranches-file ${tranches} \\
        --create-output-variant-index true \\
        --lod-score-cutoff 0.0 \\
        --exclude-filtered true \\
        -mode ${mode} \\
        -O Chr${chr}.raw.${mode.toLowerCase()}.recal.vcf.gz

    tabix -p vcf -f Chr${chr}.raw.${mode.toLowerCase()}.recal.vcf.gz
    """
}
