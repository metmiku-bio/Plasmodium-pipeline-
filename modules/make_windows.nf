process MAKE_WINDOWS {
    // Generates per-chromosome interval files in two modes:
    //
    // MODE A  --intervals_dir provided  (curated static regions)
    //   interval_list = path to core_chr<N>.list  (GATK interval format, one region per line)
    //   - whole.list  : the static file copied as-is → HaplotypeCaller -L
    //   - windows.list: the static file copied as-is → GenomicsDBImport intervals
    //                   and fanned-out line-by-line to GenotypeGVCFs
    //   bedtools is NOT used in this mode — the curated regions are the windows.
    //
    // MODE B  no --intervals_dir  (derive everything from the reference .fai)
    //   interval_list = NO_FILE sentinel
    //   - whole.list  : single "chrN:1-<len>" line extracted from the .fai
    //   - windows.list: bedtools makewindows -n <num_genome_chunks> tiles the full chromosome
    //
    // In both modes the region format in windows.list is:
    //   <contig>:<start>-<end>   (1-based, GATK-compatible)

    tag { "chr${chr}" }

    input:
    tuple val(chr), path(ref_fai), path(interval_list)   // interval_list is NO_FILE in Mode B

    output:
    tuple val(chr), path("chr${chr}_whole.list"),   emit: whole
    tuple val(chr), path("chr${chr}_windows.list"), emit: windows

    script:
    def n_chunks   = params.num_genome_chunks ?: 100
    def use_static = (interval_list.name != 'NO_FILE')
    """
    if ${use_static}; then
        # ── MODE A: curated interval list — use it directly, no re-tiling ──
        cp ${interval_list} chr${chr}_whole.list
        cp ${interval_list} chr${chr}_windows.list

    else
        # ── MODE B: derive everything from the .fai ──────────────────────────

        # Restrict .fai to this chromosome only
        awk -v c="${chr}" '\$1==c || \$1=="Chr"c || \$1=="chr"c' \
            ${ref_fai} > chr${chr}.fai

        # whole.list → single whole-chromosome interval for HaplotypeCaller
        awk '{printf "%s:1-%s\\n", \$1, \$2}' chr${chr}.fai > chr${chr}_whole.list

        # windows.list → equal tiles across the full chromosome
        bedtools makewindows -n ${n_chunks} -g chr${chr}.fai \\
            | awk '{printf "%s:%s-%s\\n", \$1, \$2+1, \$3}' \\
            > chr${chr}_windows.list
    fi
    """
}
