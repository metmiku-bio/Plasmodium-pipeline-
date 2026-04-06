process VARIANT_RECALIBRATOR {
    tag { "chr${chr}.${mode}" }

    publishDir "${params.output ?: '.'}/vqsr", mode: 'copy'

    input:
    tuple path(vcf), path(ref_fasta), path(ref_fai), path(ref_dict), val(chr), val(mode), path(resource_vcf), val(annotations)

    output:
    tuple val(chr), path("*.recal"), path("*.tranches"), path("*.plots.R"), emit: recal_data

    script:
    def max_gaussians = params.max_gaussians ?: 4
    def annot_args = annotations.join(" -an ")
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
        -an ${annot_args} \\
        -mode ${mode} \\
        --max-gaussians ${max_gaussians} \\
        --resource:Brown,known=true,training=true,truth=true,prior=15.0 ${resource_vcf} \\
        -O Chr${chr}.raw.${mode.toLowerCase()}.recal \\
        --tranches-file Chr${chr}.raw.${mode.toLowerCase()}.tranches \\
        --rscript-file Chr${chr}.raw.${mode.toLowerCase()}.plots.R \\
        --dont-run-rscript
    """
}
