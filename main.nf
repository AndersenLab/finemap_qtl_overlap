#!/usr/bin/env nextflow

//nextflow.preview.dsl=2
nextflow.enable.dsl = 2

// Needed to publish results
nextflow.preview.output = true

// import the subworkflows
include { common_strains } from './scripts/overlaps.nf'
include { make_finemap_vcf } from './scripts/overlaps.nf'
include { interval_ld } from './scripts/overlaps.nf'



// Define the workflow
workflow {
    main:
    def vcf = params.vcf
    def strains = params.strains ?: file("${workflow.projectDir}/test_data/strains.txt")
    def marker_pairs = params.marker_pairs ?: file("${workflow.projectDir}/test_data/test_peaks.csv")

    // Check required parameters
    if (!vcf || !marker_pairs) {
        error("Missing required parameters. Please provide: vcf and marker_pairs")
    }

    // Parse the marker pairs CSV file
    Channel.fromPath(marker_pairs)
        .splitCsv(header: true)
        .map { row -> [row.peak_idA, row.peakidB] }
        .set { peak_pairs }

    // Run the filter_vcf process
    filter_vcf(vcf, strains)

    // Run the calculate_ld process for each peak pair
    filtered_vcf = filter_vcf.out.vcf
    calculate_ld(filtered_vcf.combine(peak_pairs))

    publish:
    calculate_ld.out.ld_files.flatten().collect() >> "."
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
    tuple path(vcf), val(peak_a), val(peak_b)

    output:
    tuple path("${peak_a}.${peak_b}.log"), path("${peak_a}.${peak_b}.nosex"), emit: ld_files

    script:
    def peak_a_transformed = peak_a.replaceAll("_", ":")
    def peak_b_transformed = peak_b.replaceAll("_", ":")
    """
    plink --vcf ${vcf} \\
        --threads 5 \\
        --snps-only \\
        --maf 0.05 \\
        --biallelic-only \\
        --allow-extra-chr \\
        --set-missing-var-ids @:# \\
        --ld ${peak_a_transformed} ${peak_b_transformed} \\
        --out ${peak_a}.${peak_b}
    """
}
