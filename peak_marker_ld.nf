#! usr/bin/env nextflow

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


// define the workflow
workflow {
 if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}


date = new Date().format( 'yyyyMMdd' )
if (params.out == null) {
    out = "Analysis_Results-${date}"
}

 /* check if the params are set */

if (params.vcf == null) {
    println "vcf is not set"
    exit 1
}

if (params.vcf_index == null) {
    println "vcf_index is not set"
    exit 1
}

if (params.qtl_overlap == null) {
    println "qtl_overlap is not set"
    exit 1
}

 


}



