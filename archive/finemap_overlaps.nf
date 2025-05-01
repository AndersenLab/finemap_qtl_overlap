#! usr/bin/env nextflow
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

//nextflow.preview.dsl=2
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

/*
~ ~ ~ > * Parameters setup - GENERAL
*/

params.rename_chroms = "${workflow.projectDir}/bin/rename_chromosomes"
params.out = "Analysis_Results-${date}"
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.vcf = "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.impute.isotype.vcf.gz" // must have an index in the same directory
params.vcf_index = "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.impute.isotype.vcf.gz.tbi"
//params.pheno_1 = "${workflow.projectDir}/data/test_data/c_elegans/phenotypes/test_trait_1.tsv"
//params.pheno_2 = "${workflow.projectDir}/data/test_data/c_elegans/phenotypes/test_trait_2.tsv" 
//params.interval_chrom = "I"
//params.interval_start = "1000"
//params.interval_end = "100000"
//params.trait_id = "test_trait_1-test_trait_2"
params.common_strains_script = "${workflow.projectDir}/bin/common_strains.py"
params.maf = "0.05"
//params.peak_pos = "10000"

// *** Hard coded for now - copied input data from NemaScan directory
//params.genes = "${params.data_dir}/${params.species}/annotations/${params.species}.gff"
params.genes = "/projects/b1059/projects/Ryan/Toxin_GWAS/qtl_overlaps/ld_between_peak_markers/input_data/c_elegans/annotations/c_elegans.gff"

// Define the annotation file channel
//ann_file = Channel.fromPath("/projects/b1059/data/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.strain-annotation.tsv")
ann_file = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.strain-annotation.tsv")

params.sparse_cut = 0.05

include {make_finemap_vcf; interval_ld} from './scripts/overlaps.nf'

qtl_overlap = "/projects/b1059/projects/Ryan/Toxin_GWAS/qtl_overlaps/ld_between_peak_markers/data/test_data/toxin_length_finemapping_input_with_pheno_files.tsv"

// UNIQUE TO FINEMAPS_OVERLAPS
// On quest need to specify the path to the conda environment
process common_strains {
    executor 'local'
    conda '/home/rjm6024/.conda/envs/vcf_stats_1.0'
    publishDir "${params.out}/${trait_pair}/${peak}/Data", pattern: "common_strains.txt", overwrite: true

    input:
        tuple val(trait_pair), file(pheno_1), file(pheno_2), file(common_strains_script), val(peak)

    output:
        tuple val(trait_pair), val(peak), file('common_strains.txt')

    """
    python ${common_strains_script} ${pheno_1} ${pheno_2} common_strains.txt
    """
}

process prep_ld_files {
    container 'mckeowr1/prep_sims:1.1'
    executor 'slurm'

    cpus 5
    
    tag {trait_pair}

    publishDir "${params.out}/${trait_pair}/${peak}/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix_inbred.tsv"
    publishDir "${params.out}/${trait_pair}/${peak}/Data", mode: 'copy', pattern: "*LD_inbred.tsv"

    publishDir "${params.out}/${trait_pair}/${peak}/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix_loco.tsv"
    publishDir "${params.out}/${trait_pair}/${peak}/Data", mode: 'copy', pattern: "*LD_loco.tsv"

    input:
        tuple val(trait_pair), val(peak), val(chromosome), val(start_pos), val(end_pos), val(peak_pos), file(common_strains), file(num_chroms), file(imputed_vcf), file(vcf_index), val(maf), val(algorithm)

    output:
        tuple val(trait_pair), val(peak), val(peak_pos), val(chromosome), val(start_pos), val(end_pos), file("*ROI_Genotype_Matrix*.tsv"), file("*LD*.tsv"), file("*.bim"), file("*.bed"), file("*.fam"), emit: finemap_preps
        //tuple val(trait_pair), file(pheno), file("*ROI_Genotype_Matrix*.tsv"), file("*LD*.tsv"), file("*.bim"), file("*.bed"), file("*.fam"), val(algorithm), emit: finemap_preps
        tuple val(trait_pair), file("*ROI_Genotype_Matrix_inbred.tsv"), file("*LD_inbred.tsv"), emit: finemap_LD_inbred, optional: true
        tuple val(trait_pair), file("*ROI_Genotype_Matrix_loco.tsv"), file("*LD_loco.tsv"), emit: finemap_LD_loco, optional: true

    """
        bcftools view --regions ${chromosome}:${start_pos}-${end_pos} ${imputed_vcf} \
        -S ${common_strains} |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools filter -e 'GT="het"' |\\
        bcftools annotate --rename-chrs ${num_chroms} |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.txt
        bcftools view --regions ${chromosome}:${start_pos}-${end_pos} ${imputed_vcf} \
        -S ${common_strains} |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 |\\
        bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
        plink --vcf finemap.vcf.gz \\
            --threads 5 \\
            --snps-only \\
            --maf ${params.maf} \\
            --biallelic-only \\
            --allow-extra-chr \\
            --set-missing-var-ids @:# \\
            --geno \\
            --make-bed \\
            --recode vcf-iid bgz \\
            --extract ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.txt \\
            --out ${trait_pair}.${chromosome}.${start_pos}.${end_pos}
        nsnps=`wc -l ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.txt | cut -f1 -d' '`
        chrom_num=`cat ${num_chroms} | grep -w ${chromosome} | cut -f2 -d' '`
        plink --r2 with-freqs \\
            --threads 5 \\
            --allow-extra-chr \\
            --snps-only \\
            --ld-window-r2 0 \\
            --ld-snp \$chrom_num:${peak_pos} \\
            --ld-window \$nsnps \\
            --ld-window-kb 10000 \\
            --chr \$chrom_num \\
            --out ${trait_pair}.${chromosome}:${start_pos}-${end_pos}.QTL \\
            --set-missing-var-ids @:# \\
            --vcf ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.vcf.gz
        cut ${trait_pair}.${chromosome}:${start_pos}-${end_pos}.QTL.ld -f2-10 > ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.LD_${algorithm}.tsv
        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' finemap.vcf.gz |\\
            sed 's/[[# 0-9]*]//g' |\\
            sed 's/:GT//g' |\\
            sed 's/0|0/-1/g' |\\
            sed 's/1|1/1/g' |\\
            sed 's/0|1/NA/g' |\\
            sed 's/1|0/NA/g' |\\
            sed 's/.|./NA/g'  |\\
            sed 's/0\\/0/-1/g' |\\
            sed 's/1\\/1/1/g'  |\\
            sed 's/0\\/1/NA/g' |\\
            sed 's/1\\/0/NA/g' |\\
            sed 's/.\\/./NA/g' |\\
            sed 's/^23/X/g' > ${trait_pair}.${chromosome}:${start_pos}-${end_pos}.ROI_Genotype_Matrix_${algorithm}.tsv
    """
}

process prepare_gcta_files {

    // machineType 'n1-standard-4'
    label "large"
    container 'mckeowr1/prep_sims:1.1'
    executor 'slurm'

    cpus 5

    input:
        tuple val(trait_pair), val(peak), val(peak_pos), file(traits), file(common_strains), file(vcf), file(index), file(num_chroms)

    output:
        tuple val(trait_pair), val(peak), val(peak_pos), file(traits), file("${trait_pair}.bed"), file("${trait_pair}.bim"), file("${trait_pair}.fam"), file("${trait_pair}.map"), file("${trait_pair}.nosex"), file("${trait_pair}.ped"), file("${trait_pair}.log")

    """
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools view -S ${common_strains} -Ou |\\
    bcftools filter -i N_MISSING=0 -Oz --threads 5 -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz
    plink --vcf renamed_chroms.vcf.gz \\
          --threads 5 \\
          --snps-only \\
          --biallelic-only \\
          --maf ${params.maf} \\
          --set-missing-var-ids @:# \\
          --indep-pairwise 50 10 0.8 \\
          --geno \\
          --allow-extra-chr
    tail -n +2 ${traits} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_formated_trats.tsv
    plink --vcf renamed_chroms.vcf.gz \\
          --threads 5 \\
          --make-bed \\
          --snps-only \\
          --biallelic-only \\
          --maf ${params.maf} \\
          --set-missing-var-ids @:# \\
          --extract plink.prune.in \\
          --geno \\
          --recode \\
          --out ${trait_pair} \\
          --allow-extra-chr \\
          --pheno plink_formated_trats.tsv
    """
}

process gcta_fine_maps {
    // machineType 'n2-highmem-8'
    container 'andersenlab/nemascan:20220407173056db3227'
    executor 'slurm'
    cpus 9

    publishDir "${params.out}/${trait_pair}/${peak}/Data", mode: 'copy', pattern: "*inbred.fastGWA"
    publishDir "${params.out}/${trait_pair}/${peak}/Data", mode: 'copy', pattern: "*_genes_inbred.tsv"
    publishDir "${params.out}/${trait_pair}/${peak}/Plots", mode: 'copy', pattern: "*inbred.pdf"

    publishDir "${params.out}/${trait_pair}/${peak}/Data", mode: 'copy', pattern: "*loco.fastGWA"
    publishDir "${params.out}/${trait_pair}/${peak}/Data", mode: 'copy', pattern: "*_genes_loco.tsv"
    publishDir "${params.out}/${trait_pair}/${peak}/Plots", mode: 'copy', pattern: "*loco.pdf"

    input:
        tuple val(trait_pair), val(peak), val(peak_pos), val(chromosome), val(start_pos), val(end_pos), file(roi_geno), file(roi_ld), file(roi_bim), file(roi_bed), file(roi_fam), \
        file(traits), file(trait_pair_bed), file(trait_pair_bim), file(trait_pair_fam), file(trait_pair_map), file(trait_pair_nosex), file(trait_pair_ped), file(trait_pair_log), \
        val(algorithm),file(annotation), file(genefile), file(finemap_qtl_intervals), file(plot_genes)

    output:
        tuple file("*.fastGWA"), val(trait_pair), file("*.prLD_df*.tsv"), file("*.pdf"), file("*_genes*.tsv"), val(algorithm)
        //tuple file("*_genes*.tsv"), val(TRAIT), val(algorithm), emit: finemap_done
        //tuple val(TRAIT), file("*inbred.fastGWA"), file("*.prLD_df_inbred.tsv"), file("*_genes_inbred.tsv"), emit: finemap_html_inbred, optional: true
        //tuple val(TRAIT), file("*loco.fastGWA"), file("*.prLD_df_loco.tsv"), file("*_genes_loco.tsv"), emit: finemap_html_loco, optional: true

    """
    tail -n +2 ${traits} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_finemap_traits.tsv
    gcta64 --bfile ${trait_pair} \\
            --autosome \\
            --maf ${params.maf} \\
            --make-grm-inbred \\
            --out ${trait_pair}.FM_grm_inbred.${algorithm} \\
            --thread-num 9
    gcta64 --grm ${trait_pair}.FM_grm_inbred.${algorithm} \\
            --make-bK-sparse ${params.sparse_cut} \\
            --out ${trait_pair}.sparse_FM_grm_inbred.${algorithm}  \\
            --thread-num 9
    gcta64 --grm ${trait_pair}.FM_grm_inbred.${algorithm} \\
            --pca 1 \\
            --out ${trait_pair}.sparse_FM_grm_inbred.${algorithm}  \\
            --thread-num 9
    gcta64 --fastGWA-lmm-exact \\
            --grm-sparse ${trait_pair}.sparse_FM_grm_inbred.${algorithm} \\
            --bfile ${trait_pair}.${chromosome}.${start_pos}.${end_pos} \\
            --qcovar ${trait_pair}.sparse_FM_grm_inbred.${algorithm}.eigenvec \\
            --out ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.finemap_inbred.${algorithm} \\
            --pheno plink_finemap_traits.tsv \\
            --maf ${params.maf} \\
            --thread-num 9
    
    Rscript --vanilla ${finemap_qtl_intervals} ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.finemap_inbred.${algorithm}.fastGWA ${roi_geno} ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.LD_${algorithm}.tsv ${algorithm}
    
    Rscript --vanilla ${plot_genes} ${trait_pair}.${chromosome}.${start_pos}.${end_pos}.prLD_df_${algorithm}.tsv ${traits} ${genefile} ${annotation} ${algorithm}

    """
}

workflow{

Channel.fromPath(qtl_overlap)
    .splitCsv(header: true, sep: '\t')
    .map { row -> [row."trait_pair", file(row."traitA_file"), file(row."traitB_file"), file(params.common_strains_script), row."peak"] } | common_strains

//prepare the plink files for the gcta analysis
// Use all variants rather than just variants in the reigon
// Allows us to construcut the GRM
Channel.fromPath(qtl_overlap)
    .splitCsv(header: true, sep: '\t')
    .map { row -> [row."trait_pair", row."peak", row."peak_pos", file(row."trait_file")]}
    //join to the common strains output using the trait_pair data
    .join(common_strains.out, by:[0, 1])
    .combine(Channel.fromPath(params.vcf))
    .combine(Channel.fromPath(params.vcf_index))
    .combine(Channel.fromPath(params.rename_chroms)) | prepare_gcta_files


//common_strains.out.view()
// get the common strains between the two phenotypesChannel
Channel.fromPath(qtl_overlap)
    .splitCsv(header: true, sep: '\t')// input order for finemap LD files
    .map { row -> [row."trait_pair", row."peak", row."CHROM_A", row."leftmost", row."rightmost", row."peak_pos" ]}
    // join the common strains with the phenotype files using the first two values of the common strains output and the first two values from the map output (trait pair name and the peakPOS either peakPOS_A or peakPOS_B)
    .join(common_strains.out, by:[0, 1]) 
    .combine(Channel.fromPath(params.rename_chroms))
    .combine(Channel.fromPath(params.vcf))
    .combine(Channel.fromPath(params.vcf_index))
    .combine(Channel.from(params.maf))
    .combine(Channel.from("inbred"))| prep_ld_files
//tuple file(overlaps), file(common_strains), file(phenotype_A), file(pheontype_B), file(num_chroms), file(imputed_vcf), val(algorithm), val(maf)

// Join the outputs from the prepare_gcta_files and the prep_ld_files
// This is the input for the gcta_fine_maps process

prep_ld_files.out.finemap_preps
    .join(prepare_gcta_files.out, by: [0, 1, 2]) 
    .combine(Channel.from("inbred")) // *** THIS IS A SPACE FILLER current inputs don't have detection algorithm to pass 
    .combine(ann_file)
    .combine(Channel.fromPath("${params.genes}"))
    .combine(Channel.fromPath("${params.bin_dir}/Finemap_QTL_Intervals.R"))
    .combine(Channel.fromPath("${params.bin_dir}/plot_genes.R"))| gcta_fine_maps
    
}