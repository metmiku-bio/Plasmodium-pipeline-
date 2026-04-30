process PCA {
    tag "pca_analysis"

    publishDir "${params.output ?: '.'}/pca", mode: 'copy'

    input:
    tuple path(combined_vcf), path(metadata), val(color_by), val(pc_count), val(plink_memory)

    output:
    path("plink_out_*"), emit: pca_outputs
    path("pca_plot_*.pdf"), emit: pca_plots_pdf
    path("pca_plot_*.png"), emit: pca_plots_png
    path("*.newick"), emit: tree_file
    path("pca_summary.txt"), emit: pca_summary
    path("*_pca_results.rds"), emit: pca_rds
    path("variance_explained.pdf"), emit: variance_plot
    path("scree_plot.pdf"), emit: scree_plot

    script:
    def color_param = color_by ?: "country"
    def pc_param = pc_count ?: 10
    def mem_param = plink_memory ?: 8000
    
    """
    echo "=========================================="
    echo "RUNNING PCA ANALYSIS"
    echo "Input VCF:      ${combined_vcf}"
    echo "Metadata:       ${metadata}"
    echo "Color by:       ${color_param}"
    echo "PCs to compute: ${pc_param}"
    echo "PLINK memory:   ${mem_param} MB"
    echo "=========================================="
    
    # Filter to SNPs only
    bcftools view -v snps ${combined_vcf} -o snps_only.vcf.gz -O z 2>&1 | tee bcftools_filter.log
    bcftools index snps_only.vcf.gz 2>&1 | tee -a bcftools_filter.log
    
    # Run PLINK PCA
    plink --vcf snps_only.vcf.gz \\
          --pca ${pc_param} header \\
          --double-id \\
          --allow-extra-chr \\
          --memory ${mem_param} \\
          --out plink_out_ 2>&1 | tee plink_pca.log
    
    # Run R script for PCA visualization and tree building
    Rscript ${projectDir}/scripts/plink_pca.R \\
        --workdir . \\
        --prefix plink_out_ \\
        --metadata ${metadata} \\
        --color_by ${color_param} \\
        --pc_count ${pc_param} \\
        --output_dir .
    
    echo "PCA analysis completed successfully"
    echo "Results: plink_out_.eigenvec, plink_out_.eigenval, pca_plot_*.pdf, *.newick"
    """
}