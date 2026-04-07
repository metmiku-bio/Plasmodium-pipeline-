process GATHER_VCFS {
    tag { "chr${chr}" }

    publishDir "${params.output ?: '.'}/raw_vcfs", mode: 'copy'

    input:
    // vcf_parts: list of per-region VCF files, already sorted by start coordinate
    // (sorting is done in main.nf before groupTuple so the list arrives in order)
    tuple val(chr), path(vcf_parts), path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    tuple val(chr), path("Chr${chr}.raw.vcf.gz"),     emit: vcf
    tuple val(chr), path("Chr${chr}.raw.vcf.gz.tbi"), emit: vcf_idx

    script:
    def mem_gb  = Math.max(task.memory.toGiga() - 4, 2)
    // Build -I args in the order the files arrive (caller guarantees genomic order)
    def vcf_args = vcf_parts.collect { vcf -> "-I ${vcf}" }.join(" \\\n        ")
    """
    gatk --java-options "-Xmx${mem_gb}G" GatherVcfs \\
        ${vcf_args} \\
        -O Chr${chr}.raw.vcf.gz

    tabix -p vcf -f Chr${chr}.raw.vcf.gz
    """
}
