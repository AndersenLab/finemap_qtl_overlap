#!/usr/bin/env nextflow

//nextflow.preview.dsl=2
nextflow.enable.dsl=2




// import the subworkflows
include { common_strains} from './scripts/overlaps.nf'
include { make_finemap_vcf} from './scripts/overlaps.nf'
include { interval_ld} from './scripts/overlaps.nf'



// Define the workflow
workflow {
    def vcf = params.vcf
    def strains = params.strains ?: file("${workflow.projectDir}/test_data/strains.txt")
    def peak_a = params.peak_a
    def peak_b = params.peak_b
    
    // Check required parameters
    if (!vcf || !peak_a || !peak_b) {
        error "Missing required parameters. Please provide: vcf, peak_a, and peak_b"
    }
    
    // Run the filter_vcf process
    filter_vcf(vcf, strains)
    
    // Run the calculate_ld process using the output from filter_vcf
    calculate_ld(filter_vcf.out.vcf, peak_a, peak_b)
}

// Process to filter VCF file
process filter_vcf {
    label 'filter_vcf'
    
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
    label 'calculate_ld'
    
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



