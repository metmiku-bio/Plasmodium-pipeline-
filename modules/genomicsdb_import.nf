process GENOMICSDB_IMPORT {
    tag { "chr${chr}" }

    input:
    val chr
    path gvcf_files
    path interval_list

    output:
    tuple val(chr), path("chr${chr}_database"), emit: database

    script:
    def threads = params.db_import_threads ?: 24
    def mem_gb = Math.min(task.memory.toGiga() - 4, 36)
    def gvcf_args = gvcf_files.collect { "--variant $it" }.join(" ")
    """
    gatk --java-options "-Xmx${mem_gb}G" GenomicsDBImport \\
        $gvcf_args \\
        --genomicsdb-workspace-path chr${chr}_database \\
        --intervals ${interval_list} \\
        --batch-size 100 \\
        --reader-threads ${threads} \\
        --genomicsdb-segment-size 8048576 \\
        --genomicsdb-vcf-buffer-size 160384
    """
}
