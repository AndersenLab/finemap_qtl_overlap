#!/usr/bin/env nextflow

//nextflow.preview.dsl=2
nextflow.enable.dsl=2


/*
~ ~ ~ > * Optional Parameters setup - GENERAL
*/
params.rename_chroms = "${workflow.projectDir}/bin/rename_chromosomes"
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.maf = "0.05"
params.common_strains_script = "${workflow.projectDir}/bin/common_strains.py"


/*
~ ~ ~ > * Required Parameters setup - GENERAL
*/

// import the subworkflows
include { common_strains} from './scripts/overlaps.nf'
include { make_finemap_vcf} from './scripts/overlaps.nf'
include { interval_ld} from './scripts/overlaps.nf'



// Define the workflow
workflow {
    def vcf = params.vcf
    def strains = params.strains
    def peak_a = params.peak_a
    def peak_b = params.peak_b
    
    // Check required parameters
    if (!vcf || !strains || !peak_a || !peak_b) {
        error "Missing required parameters. Please provide: vcf, strains, peak_a, and peak_b"
    }
    
    // Run the filter_vcf process
    filter_vcf(vcf, strains)
    
    // Run the calculate_ld process using the output from filter_vcf
    calculate_ld(filter_vcf.out.vcf, peak_a, peak_b)
}

// Process to filter VCF file
process filter_vcf {
    label 'bcftools_filter_vcf'
    
    input:
    path vcf
    path strains
    
    output:
    path "finemap.vcf.gz", emit: vcf
    
    script:
    """
    bcftools view -S ${strains} -Ou ${vcf} | \\
    bcftools filter -i N_MISSING=0 -o finemap.vcf.gz
    """
}

// Process to calculate LD between peak markers
process calculate_ld {
    label 'plink_recode_vcf'
    
    input:
    path vcf
    val peak_a
    val peak_b
    
    output:
    path "${peak_a}_${peak_b}.log", emit: log
    
    script:
    """
    plink --vcf ${vcf} \\
        --threads 5 \\
        --snps-only \\
        --maf 0.05 \\
        --biallelic-only \\
        --allow-extra-chr \\
        --set-missing-var-ids @:# \\
        --ld ${peak_a} ${peak_b} \\
        --out ${peak_a}_${peak_b}
    """
}



