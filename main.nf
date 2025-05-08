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
    def outdir = params.outdir ?: "${workflow.projectDir}/results"
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
        .map { row -> tuple(row.peak_idA, row.peakidB) } // Emit tuples
        .set { peak_pairs_ch }

    // Run the filter_vcf process
    filter_vcf(vcf, strains)

    // Prepare input for calculate_ld
    filtered_vcf_ch = filter_vcf.out.vcf
    input_for_calculate_ld = filtered_vcf_ch.combine(peak_pairs_ch)

    // Run the calculate_ld process for each peak pair
    calculate_ld(input_for_calculate_ld)
    
    // Extract R-squared values
    ch_for_extraction = calculate_ld.out.ld_outputs
        .map { peak_a, peak_b, log_file -> tuple(peak_a, peak_b, log_file) }
    extract_r_squared(ch_for_extraction)

    // Aggregate R-squared values and prepare CSV content
    out_ch = extract_r_squared.out.rsq_values // Channel: [peak_a, peak_b, rsq_string_from_stdout]
    
    out_ch.map { it.join(',') }
    .collectFile(name:"${outdir}/peak_ld_rsq.csv", sort:false).view()
}
output {
    "." {
        mode "copy"
    }
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
    tuple val(peak_a), val(peak_b), path("${peak_a}.${peak_b}.log"), emit: ld_outputs

    script:
    def peak_a_transformed = peak_a.replaceAll("_", ":")
    def peak_b_transformed = peak_b.replaceAll("_", ":")
    """
    plink --vcf ${vcf} \\
        --threads 6 \\
        --snps-only \\
        --maf 0.05 \\
        --biallelic-only \\
        --allow-extra-chr \\
        --set-missing-var-ids @:# \\
        --ld ${peak_a_transformed} ${peak_b_transformed} \\
        --out ${peak_a}.${peak_b}
    """
}

// New process to extract R-squared value
process extract_r_squared {
    label 'extract_r_squared'

    input:
    tuple val(peak_a), val(peak_b), path(log_file)

    output:
    tuple val(peak_a), val(peak_b), stdout, emit: rsq_values

    script:
    """
    #!/bin/bash
    rsq_val=\$(grep 'R-sq = ' "${log_file}" | awk '{print \$3}')
    echo "\$rsq_val"
    """
}

// New process to write the final CSV file
process write_final_csv {
    label 'write_final_csv'

    input:
    val csv_content

    output:
    path "peak_ld_rsq.csv", emit: file

    script:
    """
    echo "${csv_content}" > peak_ld_rsq.csv
    """
}
