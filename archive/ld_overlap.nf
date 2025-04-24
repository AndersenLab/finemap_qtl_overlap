#! usr/bin/env nextflow
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2
// nextflow.enable.dsl=2


date = new Date().format( 'yyyyMMdd' )

/*
~ ~ ~ > * Parameters setup - GENERAL
*/

//params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.rename_chroms = "${workflow.projectDir}/bin/rename_chromosomes"
params.out = "Analysis_Results-${date}"
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
//params.vcf = "${workflow.projectDir}/data/test_data/c_elegans/genotypes/c_elegans.test.rename.vcf.gz" // must have an index in the same directory
//params.pheno_1 = "${workflow.projectDir}/data/test_data/c_elegans/phenotypes/test_trait_1.tsv"
//params.pheno_2 = "${workflow.projectDir}/data/test_data/c_elegans/phenotypes/test_trait_2.tsv" 
//params.interval_chrom = "I"
//params.interval_start = "1000"
//params.interval_end = "100000"
//params.trait_id = "test_trait_1-test_trait_2"
params.common_strains_script = "${workflow.projectDir}/bin/common_strains.py"
params.maf = "0.05"
//params.peak_pos = "10000"

include {common_strains; make_finemap_vcf; interval_ld} from './scripts/overlaps.nf'


workflow {
    // Create an input channel for each trait
    trait_ch1 = Channel.fromPath(params.pheno_1)
    trait_ch2 = Channel.fromPath(params.pheno_2)

    traits = trait_ch1.combine(trait_ch2) 
    
    (Channel.fromPath(params.common_strains_script))
        .combine(traits) | common_strains
    

    // create trait_id by combining the filenames of trait 1 and trait 2
    //trait_id = "${trait_ch1.baseName}-${trait_ch2.baseName}"
    // print the trait_id to the console
    //println trait_id

    Channel.fromPath(params.vcf)
        .combine(common_strains.out)
        .combine(Channel.of(params.interval_chrom))
        .combine(Channel.of(params.interval_start))
        .combine(Channel.of(params.interval_end))
        .combine(Channel.fromPath(params.rename_chroms))
        .combine(Channel.of(params.trait_id)) | make_finemap_vcf
    
    make_finemap_vcf.out
        .combine(Channel.of(params.interval_chrom))
        .combine(Channel.of(params.interval_start))
        .combine(Channel.of(params.interval_end))
        .combine(Channel.of(params.peak_pos))
        .combine(Channel.of(params.maf))
       .combine(Channel.of(params.trait_id)) | interval_ld


}