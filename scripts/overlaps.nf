process common_strains {
    publishDir "${params.out}/strain_data", pattern: "common_strains.txt", overwrite: true

    input:
        tuple file(common_strains_script), file(pheno_1), file(pheno_2)

    output:
        file('common_strains.txt')

    """
    module load python-miniconda3
    source activate /home/rjm6024/.conda/envs/vcf_stats_1.0
    python ${common_strains_script} ${pheno_1} ${pheno_2} common_strains.txt
    """
}

process make_finemap_vcf {
    publishDir "${params.out}/vcf", pattern: "finemap.vcf.gz", overwrite: true
    publishDir "${params.out}/vcf", pattern: "*.txt", overwrite: true

    input:
        tuple file(vcf), file(common_strains), val(interval_chrom), val(interval_start), val(interval_end), file(rename_chroms), val(trait_id)

    output:
        tuple file("finemap.vcf.gz"), file("${trait_id}.${interval_chrom}.${interval_start}.${interval_end}.txt")
    
    """
    module load python-miniconda3
    source activate /projects/b1059/software/conda_envs/bcftools

    bcftools index ${vcf}

    bcftools view --regions ${interval_chrom}:${interval_start}-${interval_end} -S ${common_strains} -Ou ${vcf} |\\
    bcftools filter -i N_MISSING=0 -o finemap.vcf.gz   
    
    bcftools filter -e 'GT="het"' finemap.vcf.gz | \\
    awk '\$0 !~ "#" {print \$1":"\$2}' > ${trait_id}.${interval_chrom}.${interval_start}.${interval_end}.txt
    """
}
process interval_ld {
    
    container 'mckeowr1/prep_sims:1.1'

    publishDir "${params.out}/ld", overwrite: true
    
    input:
        tuple file(finemap_vcf), file(passing_vars), val(interval_chrom), val(interval_start), val(interval_end), val(peak_pos), val(maf), val(trait_id)

    output:
        path("${trait_id}.${interval_chrom}:${interval_start}-${interval_end}.QTL.ld") 
    """
        plink --vcf ${finemap_vcf} \\
            --threads 5 \\
            --snps-only \\
            --maf ${maf} \\
            --biallelic-only \\
            --allow-extra-chr \\
            --set-missing-var-ids @:# \\
            --geno \\
            --make-bed \\
            --recode vcf-iid bgz \\
            --extract ${passing_vars} \\
            --out ${trait_id}.${interval_chrom}.${interval_start}.${interval_end}
        
        nsnps=`wc -l ${passing_vars} | cut -f1 -d' '`

        plink --r2 with-freqs \\
            --threads 5 \\
            --allow-extra-chr \\
            --snps-only \\
            --ld-window-r2 0 \\
            --ld-snp ${interval_chrom}:${peak_pos} \\
            --ld-window \$nsnps \\
            --ld-window-kb 6000 \\
            --chr ${interval_chrom} \\
            --out ${trait_id}.${interval_chrom}:${interval_start}-${interval_end}.QTL \\
            --set-missing-var-ids @:# \\
            --vcf ${trait_id}.${interval_chrom}.${interval_start}.${interval_end}.vcf.gz
    """
}

process prep_ld_files {

    // machineType 'n1-standard-4'
    label 'med'

    tag {TRAIT}

    publishDir "${params.out}/INBRED/Fine_Mappings/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix_inbred.tsv"
    publishDir "${params.out}/INBRED/Fine_Mappings/Data", mode: 'copy', pattern: "*LD_inbred.tsv"

    publishDir "${params.out}/LOCO/Fine_Mappings/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix_loco.tsv"
    publishDir "${params.out}/LOCO/Fine_Mappings/Data", mode: 'copy', pattern: "*LD_loco.tsv"

    input:
        tuple file(overlaps), file(common_strains), file(phenotype_A), file(pheontype_B), file(num_chroms), file(imputed_vcf), val(algorithm), val(maf)

    output:
        tuple val(TRAIT), file(pheno), file("*ROI_Genotype_Matrix*.tsv"), file("*LD*.tsv"), file("*.bim"), file("*.bed"), file("*.fam"), val(algorithm), emit: finemap_preps
        tuple val(TRAIT), file("*ROI_Genotype_Matrix_inbred.tsv"), file("*LD_inbred.tsv"), emit: finemap_LD_inbred, optional: true
        tuple val(TRAIT), file("*ROI_Genotype_Matrix_loco.tsv"), file("*LD_loco.tsv"), emit: finemap_LD_loco, optional: true

    """
        filename=${overlaps}
        echo Start
        while read p; do
        #if it's the first line, skip it
        if [[ \$p == *'CHROM_A'* ]]; then
            continue
        fi
            chromosome=`echo \$p | cut -f1 -d' '`
            trait=`echo \$p | cut -f2 -d' '`
            start_pos=`echo \$p | cut -f3 -d' '`
            peak_id=`echo \$p | cut -f5 -d' '` #Denotes if the peakmarker positon is trait A (first name before "." delim in trait or trait B peak
            peak_pos=`echo \$p | cut -f6 -d' '`
            end_pos=`echo \$p | cut -f4 -d' '`
      if [ \$chromosome == "MtDNA" ]; then
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S ${common_strains} |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools annotate --rename-chrs ${num_chroms} |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S ${common_strains} |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 |\\
        bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
      else
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S ${common_strains} |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools filter -e 'GT="het"' |\\
        bcftools annotate --rename-chrs ${num_chroms} |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S ${common_strains} |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 |\\
        bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
      fi
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
            --extract \$trait.\$chromosome.\$start_pos.\$end_pos.txt \\
            --out \$trait.\$chromosome.\$start_pos.\$end_pos
        nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d' '`
        chrom_num=`cat ${num_chroms} | grep -w \$chromosome | cut -f2 -d' '`
        plink --r2 with-freqs \\
            --threads 5 \\
            --allow-extra-chr \\
            --snps-only \\
            --ld-window-r2 0 \\
            --ld-snp \$chrom_num:\$peak_pos \\
            --ld-window \$nsnps \\
            --ld-window-kb 6000 \\
            --chr \$chrom_num \\
            --out \$trait.\$chromosome:\$start_pos-\$end_pos.QTL \\
            --set-missing-var-ids @:# \\
            --vcf \$trait.\$chromosome.\$start_pos.\$end_pos.vcf.gz
        cut \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld -f2-10 > \$trait.\$chromosome.\$start_pos.\$end_pos.LD_${algorithm}.tsv
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
            sed 's/^23/X/g' > \$trait.\$chromosome:\$start_pos-\$end_pos.ROI_Genotype_Matrix_${algorithm}.tsv
        done < \$filename
    """
}
