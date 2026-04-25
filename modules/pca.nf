process PCA {
    tag "pca_analysis"

    publishDir "${params.output ?: '.'}/pca", mode: 'copy'

    input:
    tuple path(combined_vcf), path(metadata)

    output:
    path("plink_out_*"), emit: pca_outputs

    script:
    """
    # Filter to SNPs only
    bcftools view -v snps ${combined_vcf} -o snps_only.vcf.gz -O z
    bcftools index snps_only.vcf.gz

    # Run PLINK distance matrix
    plink --vcf snps_only.vcf.gz --distance square --double-id --allow-extra-chr --out plink_out_

    # Run R script for PCA
    Rscript ${projectDir}/scripts/plink_pca.R --workdir . --prefix plink_out_ --metadata ${metadata} --color_by country
    """
}