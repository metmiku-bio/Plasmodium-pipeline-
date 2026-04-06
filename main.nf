nextflow.enable.dsl=2

include { HAPLOTYPE_CALLER } from './modules/haplotype_caller'
include { GENOTYPE_GVCFS } from './modules/genotype_gvcf'
include { VARIANT_RECALIBRATOR } from './modules/variant_recalibrator'
include { APPLY_VQSR } from './modules/apply_vqsr'
include { VARIANT_RECALIBRATOR as VARIANT_RECALIBRATOR_SNP } from './modules/variant_recalibrator'
include { APPLY_VQSR as APPLY_VQSR_SNP } from './modules/apply_vqsr'

// ===================== Parameters =====================
params.bam = "*.bam"
params.ref = ""
params.resource_vcf = ""
params.chromosomes = "13"
params.ploidy = 6
params.output = "./result"

// ===================== Validation =====================
if (!params.ref) exit 1, "Error: --ref is required"
if (!params.resource_vcf) exit 1, "Error: --resource_vcf is required for VQSR"

// ===================== Reference Files =====================
ref_fasta = file(params.ref)
ref_fai  = file(params.ref + ".fai")
ref_dict = file(params.ref.replaceAll(/\.fasta$/, ".dict"))

// ===================== Chromosome Channel =====================
chr_list = params.chromosomes.toString().split(",")
chr_ch = Channel.from(chr_list)

// ===================== Main Workflow =====================
workflow {

    // --- BAM input channel ---
    bam_files = Channel.fromPath(params.bam)
        | map { bam ->
            prefix = bam.baseName.replaceAll(/\.(sorted|bqsr)?\.(dup|pfbam|mkdup)?\.bam$/, "")
            tuple(prefix, bam)
          }

    // ========== 1. HaplotypeCaller: per sample x per chromosome ==========
    bam_files
        | combine(chr_ch)
        | map { prefix, bam, chr ->
            chr_int = chr as int
            bai = file("${bam}".replaceAll(/\.bam$/, ".bai"))
            interval = file("core_chr${chr}.list")
            tuple(prefix, chr_int, bam, bai, ref_fasta, ref_fai, ref_dict, interval)
          }
        | HAPLOTYPE_CALLER
        | set { hc_out }

    // ========== 2. GenotypeGVCFs: per chromosome ==========
    hc_out.gvcf
        | map { prefix, chr, gvcf ->
            gvcf_idx = file("${gvcf}".replaceAll(/\.gz$/, ".tbi"))
            tuple(chr, gvcf, gvcf_idx, ref_fasta, ref_fai, ref_dict)
          }
        | GENOTYPE_GVCFS
        | set { gt_out }

    // ========== 3. VQSR: Indels ==========
    gt_out.vcf
        | map { t ->
            def chr = t[0]
            def vcf = t[1]
            tuple(vcf, ref_fasta, ref_fai, ref_dict, chr as String, "INDEL", params.resource_vcf, ["QD", "DP", "FS", "SOR", "MQ"])
          }
        | VARIANT_RECALIBRATOR
        | set { indel_recal_out }

    // ApplyVQSR for Indels
    indel_recal_out.recal_data
        | map { t ->
            def chr = t[0]
            def recal = t[1]
            def tranches = t[2]
            def vcf = file("Chr${chr}.raw.vcf.gz")
            tuple(vcf, recal, tranches, chr, "INDEL")
          }
        | APPLY_VQSR
        | set { indel_out }

    // ========== 4. VQSR: SNPs ==========
    indel_out.vcf
        | map { vcf ->
            def chr = vcf.name.replaceAll(".*Chr([0-9]+).*", "\$1")
            tuple(vcf, ref_fasta, ref_fai, ref_dict, chr, "SNP", params.resource_vcf, ["QD", "DP", "FS", "SOR", "MQ"])
          }
        | VARIANT_RECALIBRATOR_SNP
        | set { snp_recal_out }

    // ApplyVQSR for SNPs
    snp_recal_out.recal_data
        | map { t ->
            def chr = t[0]
            def recal = t[1]
            def tranches = t[2]
            def vcf = file("Chr${chr}.raw.indel.recal.vcf.gz")
            tuple(vcf, recal, tranches, chr, "SNP")
          }
        | APPLY_VQSR_SNP
        | set { final_vcfs }

    // ========== Output Summary ==========
    hc_out.gvcf.subscribe { t -> if(t && t[2]) println "✓ gVCF: ${t[2].name}" }
    gt_out.vcf.subscribe { t -> if(t && t[1]) println "✓ Raw VCF: ${t[1].name}" }
    final_vcfs.vcf.subscribe { vcf -> if(vcf) println "✓ Final VCF: ${vcf.name}" }
}
