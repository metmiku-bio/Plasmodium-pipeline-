process GATHER_ALL_VCFS {
    tag "gather_all"

    publishDir "${params.output ?: '.'}/combined_vcfs", mode: 'copy'

    input:
    tuple val(dummy), path(vcf_list), path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    path("all_chromosomes.vcf.gz"), emit: combined_vcf

    script:
    def mem_gb  = Math.max(task.memory.toGiga() - 4, 2)
    // Sort VCFs by chromosome number for proper order
    def sorted_vcfs = vcf_list.sort { a, b ->
        def chr_a = (a =~ /Chr(\d+)/)[0][1] as int
        def chr_b = (b =~ /Chr(\d+)/)[0][1] as int
        chr_a <=> chr_b
    }
    def vcf_args = sorted_vcfs.collect { vcf -> "-I ${vcf}" }.join(" \\\n        ")
    """
    gatk --java-options "-Xmx${mem_gb}G" GatherVcfs \\
        ${vcf_args} \\
        -O all_chromosomes.vcf.gz

    tabix -p vcf -f all_chromosomes.vcf.gz
    """
}