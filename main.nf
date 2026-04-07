nextflow.enable.dsl=2

include { HAPLOTYPE_CALLER }                               from './modules/haplotype_caller'
include { COMBINE_GVCFS }                                  from './modules/combine_gvcfs'
include { GENOTYPE_GVCF }                                  from './modules/genotype_gvcf'
include { VARIANT_RECALIBRATOR }                           from './modules/variant_recalibrator'
include { APPLY_VQSR }                                     from './modules/apply_vqsr'
include { VARIANT_RECALIBRATOR as VARIANT_RECALIBRATOR_SNP } from './modules/variant_recalibrator'
include { APPLY_VQSR as APPLY_VQSR_SNP }                   from './modules/apply_vqsr'

// ===================== Parameters =====================
params.bam           = "*.bam"
params.ref           = ""
params.resource_vcf  = ""
params.chromosomes   = "1,2,3,4,5,6,7,8,9,10,11,12,13,14"
params.ploidy        = 6
params.output        = "./results"
params.max_gaussians = 4

// ===================== Validation =====================
if (!params.ref)          exit 1, "Error: --ref is required"
if (!params.resource_vcf) exit 1, "Error: --resource_vcf is required for VQSR"

// ===================== Main Workflow =====================
workflow {

    // ---- Reference files (declared inside workflow with def) ----
    def ref_fasta = file(params.ref)
    def ref_fai   = file(params.ref + ".fai")
    def ref_dict  = file(params.ref.replaceAll(/\.fasta$/, ".dict"))
    def res_vcf   = file(params.resource_vcf)

    // ---- Chromosome channel ----
    // FIX: channel.of() replaces deprecated Channel.from()
    def chr_ch = channel.of( params.chromosomes.toString().split(",") ).flatten()

    // ---- BAM input channel ----
    // FIX: channel.fromPath() replaces deprecated Channel.fromPath()
    // FIX: def used for all variables inside closures
    def bam_files = channel.fromPath(params.bam)
        | map { bam ->
            def prefix = bam.baseName.replaceAll(/\.(sorted|bqsr|dup|pfbam|mkdup)*\.bam$/, "")
            tuple(prefix, bam)
          }

    // ========== 1. HaplotypeCaller: per sample x per chromosome ==========
    bam_files
        | combine(chr_ch)
        | map { prefix, bam, chr ->
            def chr_int  = chr as int
            def bai      = file("${bam}".replaceAll(/\.bam$/, ".bai"))
            def interval = file("core_chr${chr}.list")
            tuple(prefix, chr_int, bam, bai, ref_fasta, ref_fai, ref_dict, interval)
          }
        | HAPLOTYPE_CALLER
        | set { hc_out }

    // ========== 2. Combine gVCFs: per chromosome ==========
    hc_out.gvcf
        | map { prefix, chr, gvcf -> tuple(chr, gvcf) }
        | groupTuple()
        | map { chr, gvcf_list ->
            tuple(chr, gvcf_list, ref_fasta, ref_fai, ref_dict)
          }
        | COMBINE_GVCFS
        | set { combined_out }

    // ========== 3. GenotypeGVCFs: per chromosome ==========
    // COMBINE_GVCFS emits: tuple(chr, combined_gvcf)
    // GENOTYPE_GVCF expects: tuple(database, ref_fasta, ref_fai, ref_dict, chr)
    combined_out.combined_gvcf
        | map { chr, combined_gvcf ->
            tuple(combined_gvcf, ref_fasta, ref_fai, ref_dict, chr)
          }
        | GENOTYPE_GVCF
        | set { gt_out }
    // GENOTYPE_GVCF emits: tuple val(chr), val(region), path(vcf)

    // ========== 4. VQSR: Indels ==========
    // FIX: replaced fragile map-into-Groovy-Map pattern with proper join()
    def indel_annotations = ["QD", "DP", "FS", "SOR", "MQ"]

    gt_out.vcf
        | map { chr, region, raw_vcf ->
            // VARIANT_RECALIBRATOR input: tuple path(vcf), path(ref), path(fai), path(dict), val(chr), val(mode), path(resource_vcf), val(annotations)
            tuple(raw_vcf, ref_fasta, ref_fai, ref_dict, chr.toString(), "INDEL", res_vcf, indel_annotations)
          }
        | VARIANT_RECALIBRATOR
        | set { indel_recal_out }
    // VARIANT_RECALIBRATOR emits: tuple val(chr), path(recal), path(tranches), path(plots)

    // Key indel recal outputs by chr (String), then join with raw VCFs (also String key)
    def indel_recal_keyed = indel_recal_out.recal_data
        | map { chr, recal, tranches, plots -> tuple(chr.toString(), recal, tranches) }

    def raw_vcf_keyed = gt_out.vcf
        | map { chr, region, raw_vcf -> tuple(chr.toString(), raw_vcf) }

    // join() on chr key — safe, no fragile Groovy Map collection
    indel_recal_keyed
        | join(raw_vcf_keyed)
        // now: tuple(chr, recal, tranches, raw_vcf)
        | map { chr, recal, tranches, raw_vcf ->
            tuple(raw_vcf, recal, tranches, chr, "INDEL")
          }
        | APPLY_VQSR
        | set { indel_out }
    // APPLY_VQSR emits: path(vcf), path(vcf_idx)  [no chr — keyed by filename]

    // ========== 5. VQSR: SNPs ==========
    // FIX: re-attach chr by parsing the deterministic output filename,
    //      then use join() the same way as the indel step.
    def snp_annotations = ["QD", "DP", "FS", "SOR", "MQ"]

    // Re-key indel VCFs by chr (filename: Chr<N>.raw.indel.recal.vcf.gz)
    // Use String key so join() keys match the String chr flowing from VARIANT_RECALIBRATOR output
    def indel_vcf_keyed = indel_out.vcf
        | map { vcf ->
            def chr = (vcf.name =~ /Chr(\d+)/)[0][1]   // String, e.g. "13"
            tuple(chr, vcf)
          }

    indel_vcf_keyed
        | map { chr, indel_vcf ->
            // VARIANT_RECALIBRATOR_SNP input: tuple path(vcf), path(ref), path(fai), path(dict), val(chr), val(mode), path(resource_vcf), val(annotations)
            tuple(indel_vcf, ref_fasta, ref_fai, ref_dict, chr, "SNP", res_vcf, snp_annotations)
          }
        | VARIANT_RECALIBRATOR_SNP
        | set { snp_recal_out }
    // VARIANT_RECALIBRATOR_SNP emits: tuple val(chr), path(recal), path(tranches), path(plots)

    def snp_recal_keyed = snp_recal_out.recal_data
        | map { chr, recal, tranches, plots -> tuple(chr.toString(), recal, tranches) }

    // Join SNP recal with indel VCFs by chr (both now String keys)
    snp_recal_keyed
        | join(indel_vcf_keyed)
        // now: tuple(chr, recal, tranches, indel_vcf)
        | map { chr, recal, tranches, indel_vcf ->
            tuple(indel_vcf, recal, tranches, chr, "SNP")
          }
        | APPLY_VQSR_SNP
        | set { final_vcfs }

    // ========== Output Summary ==========
    hc_out.gvcf.subscribe { t ->
        if (t && t[2]) log.info "gVCF produced:        ${t[2].name}"
    }
    combined_out.combined_gvcf.subscribe { t ->
        if (t && t[1]) log.info "Combined gVCF:        ${t[1].name}"
    }
    gt_out.vcf.subscribe { t ->
        if (t && t[2]) log.info "Genotyped VCF:        ${t[2].name}"
    }
    indel_out.vcf.subscribe { vcf ->
        if (vcf) log.info "INDEL recal VCF:      ${vcf.name}"
    }
    final_vcfs.vcf.subscribe { vcf ->
        if (vcf) log.info "Final SNP recal VCF:  ${vcf.name}"
    }
}
