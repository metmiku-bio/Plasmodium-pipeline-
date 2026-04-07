process VARIANT_RECALIBRATOR {
    tag { "chr${chr}.${mode}" }

    publishDir "${params.output ?: '.'}/vqsr", mode: 'copy'

    input:
    tuple path(vcf), path(ref_fasta), path(ref_fai), path(ref_dict), val(chr), val(mode), path(resource_vcf), val(annotations)

    output:
    tuple val(chr), path("*.recal"), path("*.tranches"), path("*.plots.R"), emit: recal_data

    script:
    // max_gaussians: use the param value (default 4).
    // For small/regional datasets that fail to converge, pass --max_gaussians 1 on the CLI.
    def max_gaussians = params.max_gaussians ?: 4
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    """
    # Index input VCF if needed
    if [ ! -f "${vcf}.tbi" ] && [ ! -f "${vcf}.idx" ]; then
        gatk IndexFeatureFile -I ${vcf}
    fi
    # Index resource VCF if needed
    if [ ! -f "${resource_vcf}.tbi" ] && [ ! -f "${resource_vcf}.idx" ]; then
        gatk IndexFeatureFile -I ${resource_vcf}
    fi

    gatk --java-options "-Xmx${mem_gb}g" VariantRecalibrator \\
        -R ${ref_fasta} \\
        -V ${vcf} \\
        --trust-all-polymorphic \\
        -an QD -an DP -an FS -an SOR -an MQ \\
        -mode ${mode} \\
        --max-gaussians ${max_gaussians} \\
        --resource:Brown,known=true,training=true,truth=true,prior=15.0 ${resource_vcf} \\
        -O Chr${chr}.raw.${mode.toLowerCase()}.recal \\
        --tranches-file Chr${chr}.raw.${mode.toLowerCase()}.tranches \\
        --rscript-file Chr${chr}.raw.${mode.toLowerCase()}.plots.R
    """
}
