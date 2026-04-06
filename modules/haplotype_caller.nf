process HAPLOTYPE_CALLER {
    tag { "${prefix}.chr${chr}" }

    publishDir "${params.output ?: '.'}/gvcfs", mode: 'copy'

    input:
    tuple val(prefix), val(chr), path(bam), path(bam_idx), path(ref_fasta), path(ref_fai), path(ref_dict), path(interval_list)

    output:
    tuple val(prefix), val(chr), path("*.g.vcf.gz"), emit: gvcf
    tuple val(prefix), val(chr), path("*.g.vcf.gz.tbi"), emit: gvcf_idx

    script:
    def threads = task.cpus ?: 16
    def mem_gb = Math.max(task.memory.toGiga() - 4, 2)
    def ploidy = params.ploidy ?: 6
    """
    gatk --java-options "-Xmx${mem_gb}G" HaplotypeCaller \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -ERC GVCF \\
        -ploidy ${ploidy} \\
        --native-pair-hmm-threads ${threads} \\
        -O ${prefix}.chr${chr}.g.vcf \\
        --assembly-region-padding 100 \\
        --max-num-haplotypes-in-population 128 \\
        --kmer-size 10 --kmer-size 25 \\
        --min-dangling-branch-length 4 \\
        --heterozygosity 0.0029 --indel-heterozygosity 0.0017 \\
        --min-assembly-region-size 100 \\
        -L ${interval_list} \\
        -mbq 5 --base-quality-score-threshold 12

    bgzip ${prefix}.chr${chr}.g.vcf
    tabix -p vcf ${prefix}.chr${chr}.g.vcf.gz
    """
}
