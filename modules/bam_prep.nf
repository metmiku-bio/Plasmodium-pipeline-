process BAM_PREP {
    tag { prefix }

    input:
    tuple val(prefix), path(bam), path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    tuple val(prefix), path("*.pfbam"), emit: bam
    tuple val(prefix), path("*.pfbam.bai"), emit: bai
    tuple val(prefix), path("*.bamstats"), emit: stats

    script:
    def threads = task.cpus ?: 8
    def mem_gb = Math.max(task.memory.toGiga() - 2, 2)
    """
    gatk --java-options "-Xmx${mem_gb}g" \\
        CleanSam -R ${ref_fasta} -I ${bam} -O ${prefix}.clean.bam

    gatk --java-options "-Xmx${mem_gb}g" \\
        SortSam -R ${ref_fasta} -I ${prefix}.clean.bam -O ${prefix}.sorted.bam -SO coordinate --CREATE_INDEX true

    gatk --java-options "-Xmx${mem_gb}g" \\
        MarkDuplicatesSpark -R ${ref_fasta} -I ${prefix}.sorted.bam -O ${prefix}.sorted.dup.bam -M ${prefix}.metrics.txt

    samtools index ${prefix}.sorted.dup.bam
    samtools flagstat ${prefix}.sorted.dup.bam > ${prefix}.bamstats

    rm ${prefix}.clean.bam ${prefix}.sorted.bam 2>/dev/null || true
    """
}
