#!/bin/bash

# export path to the newest NF version locally
export PATH="/Users/ryanmckeown/Desktop/finemap_qtl_overlap:$PATH"

#source the zshell profile
source ~/.zshrc

# run the workflow
nextflow run main.nf -profile local \
-resume \
-output-dir /Users/ryanmckeown/2021_GWA_manuscript/data/processed/qtl_overlaps/20250507_peakmarker_ld_results \
--vcf /Users/ryanmckeown/2021_GWA_manuscript/data/processed/qtl_overlaps/WI.20220216.hard-filter.isotype.vcf.gz \
--strains /Users/ryanmckeown/2021_GWA_manuscript/data/processed/phenotypes/strains.txt \
--marker_pairs /Users/ryanmckeown/2021_GWA_manuscript/data/processed/qtl_overlaps/ld_input_peak_pairs.csv
