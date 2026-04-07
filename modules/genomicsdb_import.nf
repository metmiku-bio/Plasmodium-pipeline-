process GENOMICSDB_IMPORT {
    tag { "chr${chr}" }

    // GenomicsDB workspaces are intermediate directories — not published
    // (they are consumed by GENOTYPE_GVCF and not needed as final outputs)

    input:
    // sample_map: path to a two-column TSV  <sample_name>\t<gvcf_path>
    tuple val(chr), path(sample_map), path(gvcfs), path(gvcf_tbis), path(interval_list)

    output:
    tuple val(chr), path("chr${chr}_database"), path(interval_list), emit: database

    script:
    def threads = params.db_import_threads ?: 24
    def mem_gb  = Math.max(task.memory.toGiga() - 4, 2)
    """
    gatk --java-options "-Xmx${mem_gb}G" GenomicsDBImport \\
        --sample-name-map ${sample_map} \\
        --genomicsdb-workspace-path chr${chr}_database \\
        --intervals ${interval_list} \\
        --batch-size 100 \\
        --reader-threads ${threads} \\
        --genomicsdb-segment-size 8048576 \\
        --genomicsdb-vcf-buffer-size 160384
    """
}
