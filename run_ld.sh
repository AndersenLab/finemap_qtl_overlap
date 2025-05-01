# script to run the peak_marker_ld.nf workflow on RF
module load screen
screen -S test_pipeline

# soruce NF settings from bash profile
source ~/.bash_profile

# load the nextflow conda 
conda activate /data/eande106/software/conda_envs/nf24_env

nextflow run andersenlab/finemap_qtl_overlap \
    -r <commit-hash> \
    --vcf /vast/eande106/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.vcf.gz \
    --peak_a "II:10203" \
    --peak_b "II:10204"