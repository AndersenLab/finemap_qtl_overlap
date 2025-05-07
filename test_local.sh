#!/bin/bash

# export path to the newest NF version locally
export PATH="/Users/ryanmckeown/Desktop/finemap_qtl_overlap:$PATH"

#source the zshell profile
source ~/.zshrc

# run the workflow
nextflow run main.nf -profile local \
-resume \
--vcf /Users/ryanmckeown/Desktop/finemap_qtl_overlap/WI.20220216.hard-filter.isotype.vcf.gz