process BAM_STATS {
    tag { prefix }

    input:
    tuple val(prefix), path(bam), path(ref_fasta)

    output:
    path "ReadCoverage_final.tsv", emit: coverage
    path "InsertSize_Final.txt", emit: insert_size
    path "Bam_stats_Final.tsv", emit: bam_stats

    script:
    def chroms = params.chromosomes ?: "1,2,3"
    def stat_dir = "stats_${prefix}"
    """
    mkdir -p ${stat_dir}

    # Depth of Coverage per chromosome
    for chr in ${chroms.replace(",", " ")}; do
        gatk --java-options "-Xmx${task.memory.toGiga() - 4}g" DepthOfCoverage \\
            -R ${ref_fasta} \\
            -O ${stat_dir}/chr\${chr} \\
            --omit-locus-table true \\
            -I ${bam}

        awk -F "\\t" -v OFS="\\t" -v chr="chr\${chr}" '{print \$0, chr}' ${stat_dir}/chr\${chr}.sample_summary > ${stat_dir}/chr\${chr}.sample2_summary
    done

    cat ${stat_dir}/*.sample2_summary | awk '!/sample_id/ {print \$0}' | \\
        sed '1isample_id,total,mean,third_quartile,median,first_quartile,bases_perc_above_15' > ReadCoverage_final.tsv

    # Insert Size Metrics
    gatk CollectInsertSizeMetrics -I ${bam} -O ${stat_dir}/${prefix}_insert.txt -H ${stat_dir}/${prefix}_histo.pdf -M 0.05
    awk 'FNR>=8 && FNR<=8 {print \$0}' ${stat_dir}/${prefix}_insert.txt > InsertSize_Final.txt

    # BAM Stats
    samtools stats ${bam} | grep ^SN | cut -f 2- | awk -F"\\t" '{print \$2}' > ${stat_dir}/${prefix}_bamstat.tsv
    datamash transpose < ${stat_dir}/${prefix}_bamstat.tsv | awk -F '\\t' -v OFS='\\t' -v id="${prefix}" '{ \$(NF+1) = id; print }' > Bam_stats_Final.tsv
    """
}
