nextflow.enable.dsl=2

include { MAKE_WINDOWS }                                       from './modules/make_windows'
include { HAPLOTYPE_CALLER }                                   from './modules/haplotype_caller'
include { GENOMICSDB_IMPORT }                                  from './modules/genomicsdb_import'
include { GENOTYPE_GVCF }                                      from './modules/genotype_gvcf'
include { GATHER_VCFS }                                        from './modules/gather_vcfs'
include { VARIANT_RECALIBRATOR }                               from './modules/variant_recalibrator'
include { APPLY_VQSR }                                         from './modules/apply_vqsr'
include { VARIANT_RECALIBRATOR as VARIANT_RECALIBRATOR_SNP }   from './modules/variant_recalibrator'
include { APPLY_VQSR as APPLY_VQSR_SNP }                       from './modules/apply_vqsr'
include { GATHER_ALL_VCFS }                                    from './modules/gather_all_vcfs'
include { PCA }                                                from './modules/pca'

// ===================== Parameters =====================
params.bam               = "*.bam"
params.ref               = ""
params.resource_vcf      = ""
params.chromosomes       = "1,2,3,4,5,6,7,8,9,10,11,12,13,14"
params.intervals_dir     = ""           // optional: path to dir with core_chr<N>.list files
                                        //   set → MODE A: list files used directly as intervals (no re-tiling)
                                        //   unset → MODE B: intervals generated from .fai via bedtools
params.num_genome_chunks = 100          // MODE B only: number of equal windows per chromosome
                                        //   ignored when --intervals_dir is set (the list file IS the windows)
params.ploidy            = 6
params.run_pca          = false
params.output            = "./results"
params.max_gaussians     = 4
params.metadata          = ""

// ===================== Validation =====================
if (!params.ref)          exit 1, "Error: --ref is required"
if (!params.resource_vcf) exit 1, "Error: --resource_vcf is required for VQSR"

// ===================== Main Workflow =====================
workflow {

    // ---- Reference files ----
    def ref_fasta = file(params.ref)
    def ref_fai   = file(params.ref + ".fai")
    def ref_dict  = file(params.ref.replaceAll(/\.fasta$/, ".dict"))
    def res_vcf   = file(params.resource_vcf)

    // ---- Chromosome channel ----
    def chr_ch = channel.of( params.chromosomes.toString().split(",") ).flatten()

    // ---- BAM input channel ----
    def bam_files = channel.fromPath(params.bam)
        | map { bam ->
            def prefix = bam.baseName.replaceAll(/\.(sorted|bqsr|dup|pfbam|mkdup)*\.bam$/, "")
            tuple(prefix, bam)
          }

    // ========== 0. MAKE_WINDOWS: generate intervals per chromosome ==========
    //
    // MODE A  --intervals_dir supplied:
    //   whole.list   → the static core_chr<N>.list passed straight to HaplotypeCaller
    //   windows.list → bedtools tiles WITHIN those curated regions
    //
    // MODE B  no --intervals_dir:
    //   whole.list   → single "chrN:1-<len>" line from the .fai
    //   windows.list → bedtools tiles the full chromosome from the .fai
    //
    // We resolve the static file (or a NO_FILE sentinel) here in Groovy so
    // the process input is always a path — never a conditional branch in the channel.
    def use_intervals = params.intervals_dir ? true : false

    chr_ch
        | map { chr ->
            def ivl = use_intervals
                ? file("${params.intervals_dir}/core_chr${chr}.list")
                : file("NO_FILE")            // sentinel — detected by name in the process
            tuple(chr, ref_fai, ivl)
          }
        | MAKE_WINDOWS
        | set { windows_out }
    // windows_out.whole:   tuple val(chr), path(chr<N>_whole.list)
    // windows_out.windows: tuple val(chr), path(chr<N>_windows.list)

    // ========== 1. HaplotypeCaller: per sample x per chromosome ==========
    // Combine BAMs with the whole-chromosome interval from MAKE_WINDOWS
    bam_files
        | combine(windows_out.whole)          // [prefix, bam] × [chr, whole.list]
        | map { prefix, bam, chr, interval ->
            def bai = file("${bam}".replaceAll(/\.bam$/, ".bai"))
            tuple(prefix, chr.toString(), bam, bai, ref_fasta, ref_fai, ref_dict, interval)
          }
        | HAPLOTYPE_CALLER
        | set { hc_out }

    // ========== 2. GenomicsDBImport: per chromosome ==========
    // Build the sample-name-map TSV on the fly and collect all gVCFs + indexes.
    // Join with the tiled windows list from MAKE_WINDOWS as the interval file.
    hc_out.gvcf
        | map { prefix, chr, gvcf ->
            def tbi = file("${gvcf}.tbi")
            tuple(chr.toString(), prefix, gvcf, tbi)
          }
        | groupTuple()
        | join( windows_out.windows )          // attach chr<N>_windows.list
        | map { chr, prefixes, gvcfs, tbis, window_list ->
            def map_file = file("gvcf_chr${chr}_list.tsv")
            def lines = [prefixes, gvcfs].transpose().collect { p, g -> "${p}\t${g}" }
            map_file.text = lines.join("\n") + "\n"
            tuple(chr, map_file, gvcfs, tbis, window_list)
          }
        | GENOMICSDB_IMPORT
        | set { db_out }
    // db_out.database: tuple val(chr), path(database_dir), path(window_list)

    // ========== 3. GenotypeGVCFs: per chromosome x sub-region ==========
    // Each line of the interval list file becomes one parallel GenotypeGVCFs job.
    db_out.database
        | flatMap { chr, db, interval_list ->
            interval_list.readLines()
                .collect { line -> line.trim() }
                .findAll { line -> line && !line.startsWith('@') }
                .collect { region ->
                    tuple(chr, db, region, ref_fasta, ref_fai, ref_dict)
                }
          }
        | GENOTYPE_GVCF
        | set { gt_out }
    // gt_out.vcf: tuple val(chr), val(region), path(vcf)

    // ========== 4. GatherVcfs: merge region VCFs into one per chromosome ==========
    // Sort tiles by start coordinate before handing to GatherVcfs.
    gt_out.vcf
        | map { chr, region, vcf ->
            def start = (region =~ /:(\d+)-/)[0][1] as long
            tuple(chr, start, vcf)
          }
        | groupTuple(by: 0)
        | map { chr, starts, vcfs ->
            def sorted_vcfs = [starts, vcfs].transpose()
                                             .sort { a, b -> a[0] <=> b[0] }
                                             .collect { s, v -> v }
            tuple(chr, sorted_vcfs, ref_fasta, ref_fai, ref_dict)
          }
        | GATHER_VCFS
        | set { gathered_out }
    // gathered_out.vcf: tuple val(chr), path(Chr<N>.raw.vcf.gz)

    // ========== 5. VQSR — Indels ==========
    def indel_annotations = ["QD", "DP", "FS", "SOR", "MQ"]

    gathered_out.vcf
        | map { chr, raw_vcf ->
            tuple(raw_vcf, ref_fasta, ref_fai, ref_dict, chr.toString(), "INDEL", res_vcf, indel_annotations)
          }
        | VARIANT_RECALIBRATOR
        | set { indel_recal_out }

    def indel_recal_keyed = indel_recal_out.recal_data
        | map { chr, recal, tranches, plots -> tuple(chr.toString(), recal, tranches) }

    def raw_vcf_keyed = gathered_out.vcf
        | map { chr, raw_vcf -> tuple(chr.toString(), raw_vcf) }

    indel_recal_keyed
        | join(raw_vcf_keyed)
        | map { chr, recal, tranches, raw_vcf ->
            tuple(raw_vcf, recal, tranches, chr, "INDEL")
          }
        | APPLY_VQSR
        | set { indel_out }

    // ========== 6. VQSR — SNPs ==========
    def snp_annotations = ["QD", "DP", "FS", "SOR", "MQ"]

    def indel_vcf_keyed = indel_out.vcf
        | map { vcf ->
            def chr = (vcf.name =~ /Chr(\d+)/)[0][1]
            tuple(chr, vcf)
          }

    indel_vcf_keyed
        | map { chr, indel_vcf ->
            tuple(indel_vcf, ref_fasta, ref_fai, ref_dict, chr, "SNP", res_vcf, snp_annotations)
          }
        | VARIANT_RECALIBRATOR_SNP
        | set { snp_recal_out }

    def snp_recal_keyed = snp_recal_out.recal_data
        | map { chr, recal, tranches, plots -> tuple(chr.toString(), recal, tranches) }

    snp_recal_keyed
        | join(indel_vcf_keyed)
        | map { chr, recal, tranches, indel_vcf ->
            tuple(indel_vcf, recal, tranches, chr, "SNP")
          }
        | APPLY_VQSR_SNP
        | set { final_vcfs }

    // ========== 7. Gather all chromosomes into one VCF ==========
    final_vcfs.vcf
        | map { chr, vcf -> vcf }
        | collect
        | map { vcfs -> tuple("all", vcfs, ref_fasta, ref_fai, ref_dict) }
        | GATHER_ALL_VCFS
        | set { combined_vcf_out }

    if (params.run_pca) {
        def metadata = file(params.metadata)

        // ========== 8. PCA Analysis on SNPs ==========
        combined_vcf_out
            | combine(channel.of(metadata))
            | map { vcf, meta -> [vcf, meta] }
            | PCA
    }

    // ========== Progress logging ==========
    hc_out.gvcf.subscribe { t ->
        if (t && t[2]) log.info "gVCF produced:          ${t[2].name}"
    }
    db_out.database.subscribe { chr, db, ivl ->
        log.info "GenomicsDB built:       chr${chr}"
    }
    gathered_out.vcf.subscribe { chr, vcf ->
        log.info "Gathered VCF:           ${vcf.name}"
    }
    indel_out.vcf.subscribe { vcf ->
        if (vcf) log.info "INDEL recal VCF:        ${vcf.name}"
    }
    final_vcfs.vcf.subscribe { vcf ->
        if (vcf) log.info "Final SNP recal VCF:    ${vcf.name}"
    }
}
