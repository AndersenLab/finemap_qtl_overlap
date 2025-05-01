# script to run the peak_marker_ld.nf workflow on RF
module load screen
screen -S test_pipeline

# soruce NF settings from bash profile
source ~/.bash_profile

# load the nextflow conda 
conda activate /data/eande106/software/conda_envs/nf24_env

nextflow run peak_marker_ld.nf \
    --vcf /path/to/your.vcf \
    --peak_a "II:10203" \
    --peak_b "II:10204"