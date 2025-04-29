Code to look at the overlap between the peak markers of two different traits.


# `finemap_overlaps.nf`

Nextflow workflow to calculate LD between peak markers for each trait pair and finemap the overlapping intervals.

## Inputs

The main input is the output from `prep_ld_input.R` which is a file with the following columns:

To-do:
- [ ] Add code to create the input file. 
- [ ] Update the `nextflow.config` file to match RF specific params. 
- [ ] Remove junk data from the script. 
- [ ] Optimize the output datastructure. (Low priority given the parsing scripts are designed to work with the current output)


### Parameters
- `qtl_overlap` - This is the main input to the pipeline. It is a path to the `.tsv`input file for overlaping interval pairs.
    The pipeline expects the following columns:
    - `trait_pair` - name of the trait pair.
    - `traitA_file` - path to the `.tsv` trait file for trait A.
    - `traitB_file` - path to the `.tsv` trait file for trait B.
    - `peak` - name of the peak marker of interest.
    - `peak_pos` - position of the peak marker of interest. Either peakPOS_A or peakPOS_B.
    - `CHROM_A` - chromosome of the peak marker of trait A. (Should always be the same chromosome for trait B because inputs are assumued to be overlapping)
    - `leftmost` - leftmost position of the interval of interest.
    - `rightmost` - rightmost position of the interval of interest.
    Each overlapping interval should have two rows of data, one where each trait is the trait of interst and one where each trait is the other trait.
    

- `out` - name of the output directory.
- `maf` - minor allele frequency threshold
- `sparse_cut` - sparse cut off for the LD calculation (Default = 0.05)
- `vcf` & `vcf_index` - path to WI VCF file. [ ] Check if this is supposed to be the imputed VCF or the hardfiltered vcf. `vcf_index` is the `.tbi` file.
- `ann_file` - path to the BCSQ annotation file used in the finemapping. 




# Preparing the input data file for Nextflow script to calculate LD between overlapping peaks

`code/qtl_overlaps/ld_between_peak_markers/20240108_add_phenos_test.Rmd` is currently used to prep inputs for finemapping NF pipeline or the LD between overlapping peaks NF pipeline.

[ ] - Remove date and convert from a markdownfile

# `ld_overlap.nf`

 Nextflow Script to calculate LD between overlapping peaks and perform finemapping of the overlapping intervals
This code is based on the this NF process:
```
process prep_ld_files {

    // machineType 'n1-standard-4'
    label 'med'

    tag {TRAIT}

    publishDir "${params.out}/INBRED/Fine_Mappings/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix_inbred.tsv"
    publishDir "${params.out}/INBRED/Fine_Mappings/Data", mode: 'copy', pattern: "*LD_inbred.tsv"

    publishDir "${params.out}/LOCO/Fine_Mappings/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix_loco.tsv"
    publishDir "${params.out}/LOCO/Fine_Mappings/Data", mode: 'copy', pattern: "*LD_loco.tsv"

    input:
        tuple val(TRAIT), val(CHROM), val(marker), val(log10p), val(start_pos), val(peak_pos), val(end_pos), val(peak_id), val(h2), val(algorithm), file(geno), file(pheno), file(aggregate_mapping), file(imputed_vcf), file(imputed_index), file(phenotype), file(num_chroms)

    output:
        tuple val(TRAIT), file(pheno), file("*ROI_Genotype_Matrix*.tsv"), file("*LD*.tsv"), file("*.bim"), file("*.bed"), file("*.fam"), val(algorithm), emit: finemap_preps
        tuple val(TRAIT), file("*ROI_Genotype_Matrix_inbred.tsv"), file("*LD_inbred.tsv"), emit: finemap_LD_inbred, optional: true
        tuple val(TRAIT), file("*ROI_Genotype_Matrix_loco.tsv"), file("*LD_loco.tsv"), emit: finemap_LD_loco, optional: true

    """
        echo "HELLO"
        cat ${aggregate_mapping} |\\
        awk '\$0 !~ "\\tNA\\t" {print}' |\\
        awk '!seen[\$1,\$12,\$19,\$20,\$21]++' |\\
        awk 'NR>1{print \$1, \$11, \$18, \$19, \$20}' OFS="\\t" > ${TRAIT}_QTL_peaks.tsv
        filename='${TRAIT}_QTL_peaks.tsv'
        echo Start
        while read p; do 
            chromosome=`echo \$p | cut -f1 -d' '`
            trait=`echo \$p | cut -f2 -d' '`
            start_pos=`echo \$p | cut -f3 -d' '`
            peak_pos=`echo \$p | cut -f4 -d' '`
            end_pos=`echo \$p | cut -f5 -d' '`
      if [ \$chromosome == "MtDNA" ]; then
        cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools annotate --rename-chrs ${num_chroms} |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 |\\
        bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
      else
        cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools filter -e 'GT="het"' |\\
        bcftools annotate --rename-chrs ${num_chroms} |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
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

```

The first bit of the function looks like it sets this process to iterate across the significant peak markers from the raw mapping data. - This is not something we will need to look at LD between peak markers
```
        echo "HELLO"
        cat ${aggregate_mapping} |\\
        awk '\$0 !~ "\\tNA\\t" {print}' |\\
        awk '!seen[\$1,\$12,\$19,\$20,\$21]++' |\\
        awk 'NR>1{print \$1, \$11, \$18, \$19, \$20}' OFS="\\t" > ${TRAIT}_QTL_peaks.tsv
        filename='${TRAIT}_QTL_peaks.tsv'
        echo Start
        while read p; do 
            chromosome=`echo \$p | cut -f1 -d' '`
            trait=`echo \$p | cut -f2 -d' '`
            start_pos=`echo \$p | cut -f3 -d' '`
            peak_pos=`echo \$p | cut -f4 -d' '`
            end_pos=`echo \$p | cut -f5 -d' '`
```

Next, a series of `bcftools` commands are run to generate the `finemap.vcf` file which is used in downstream processes. * we will certainly need this * 

These bcftools commands require a few inputs 
1) phenotype data
2) chromosome variable
3) start position variable
4) end position variable
5) path to chromosome rename key
6) path to the imputed vcf file * We may want to use imputed vcf files here * 
```{bash}
        cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools filter -e 'GT="het"' |\\
        bcftools annotate --rename-chrs ${num_chroms} |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 |\\
        bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
```

Following the bcftools filtering steps, the `finemap.vcf` file is converted to a plink binary file using the following command:
```{bash}
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
```

one of the plink command output files is the `n_snps` variable 
```{bash}
        nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d' '`
```
The `n_snps` variable is used to set the LD window parameter of a subsequent `plink` command.

This next step runs another plink command. Which is used to calculate the LD between the peak marker and all other markers in the region. 
```{bash}
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
```
The `--r2 with freqs` is used to calculate the LD between the peak marker and all other markers in the region.
The `ld-window` parameter is set to the `n_snps` variable calculated in the previous step. This is defined in the `plink` documentation as the number of SNPs to include in the LD window. --ld-window <max variant ct + 1>
The `ld-snp` parameter is set to the peak marker position. This parameter sets the reference SNP for the LD calculation (correlation with the peak marker for all snps in the region). 

If the interval contains both peak markers, it doesn't matter which one is the reference SNP. Since there is not a directionality to the LD calculation. So this command will work!

The output of this command is a file called `*.QTL.ld` which contains the LD between the peak marker and all other markers in the region. If we run this on the full interval overlap, it will contain the LD between the two peak markers for the trait.

To mirror this process we will need these inputs:

1) phenotype data for both traits
2) path to the chromosome rename key
3) path to the vcf file
4) Interval chrom
5) Interval start
6) Interval end
7) peak marker position for trait 1
8) peak marker position for trait 2

The first step of the ld_overlap workflow will be to find the set of strains in that are common across the traits. I will use a python script that reads in the phenotype data for both traits and writes out a list of common strains.

```{python}
import pandas as pd
import sys

pheno1 = pd.read_csv(sys.argv[1], sep='\t')
pheno2 = pd.read_csv(sys.argv[2], sep='\t')

common_strains = list(set(pheno1['strain']).intersection(set(pheno2['strain'])))

with open(sys.argv[3], 'w') as f:
    for strain in common_strains:
        f.write(strain + '\n')
```

The next step will be to run the `bcftools` commands to filter the vcf file to the common strains and generate the `finemap.vcf` file. We also will need the chromosome rename key, the interval chrom, start, and end positions, and the peak marker positions for both traits. 

```{bash}
bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S common_strains.txt |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools filter -e 'GT="het"' |\\
        bcftools annotate --rename-chrs ${num_chroms} |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S common_strains.txt |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 |\\
        bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
```

This method is ungodly inefficent.. But we will ignore that for now.

## Process enviornments
To enable the python script to run we need to load an enviornment with the pandas package. 

```{bash}
srun --account=b1042 \
--partition=genomicsguest \
--time=01:00:00  \
--mem=30G \
--pty bash -l

conda create -n ld_overlap
conda activate ld_overlap
pip install pandas
```
Until I containerize the path to this conda env will be in my home directory on QUEST.
`environment location: /home/rjm6024/.conda/envs/ld_overlap`
To enable the use of conda envs in NF we need to make sure that the `nextflow.config` file has the following line:
```{config}
conda.enabled = true
```

If this is contained in the file then we should be able to specify the conda env in the process ex) `conda 'bwa samtools multiqc'`
`
## `bcftools`
bcftools/1.10.1

## Testing with real QTL overlap data
- Example Overlap
    Trait_1: Carboxin
        - Start: 2330097
        - End: 2677403
    Trait_2: Chlorothalonil
        - Start: 1997851
        - End: 2735727
    - 
```{bash}
interval_chrom=3 #using the renamed chrom
interval_start=1997851
interval_end=2735727
trait_id="Carboxin_Chlorothalonil"

peak_marker=2527992

vcf=/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.rename.vcf.gz

pheno_1=/projects/b1059/projects/Ryan/Toxin_GWAS/qtl_overlaps/ld_between_peak_markers/data/test_data/real_overlaps/Carboxin_medianlength48hours.tsv
pheno_2=/projects/b1059/projects/Ryan/Toxin_GWAS/qtl_overlaps/ld_between_peak_markers/data/test_data/real_overlaps/Chlorothalonil_medianlength48hours.tsv

nextflow run ld_overlap.nf \
    --interval_chrom $interval_chrom \
    --interval_start $interval_start \
    --interval_end $interval_end \
    --peak_pos $peak_marker \
    --vcf $vcf \
    --pheno_1 $pheno_1 \
    --pheno_2 $pheno_2 \
    --trait_id $trait_id
```

## Running on other example QTL overlaps

Trait 1: Methyl Mercury
    - Start: 2628141
    - End: 2835866
    - Peak: 2730209
Trait 2: Silver Nitrate 7.8
    - Start: 2622088
    - End: 3571413
    - Peak: 3016298
Chrom = II

```{bash}
interval_chrom=2 #using the renamed chrom
interval_start=2622088
interval_end=3571413
trait_id="Mercury_Silver.nitrate.7.8"

peak_marker=2730209

vcf=/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.rename.vcf.gz

pheno_1=/projects/b1059/projects/Ryan/Toxin_GWAS/qtl_overlaps/ld_between_peak_markers/data/test_data/real_overlaps/Methyl.mercury_medianlength48hours.tsv

pheno_2=/projects/b1059/projects/Ryan/Toxin_GWAS/qtl_overlaps/ld_between_peak_markers/data/test_data/real_overlaps/Silver.nitrate.7.8_medianlength48hours.tsv

nextflow run ld_overlap.nf \
    --interval_chrom $interval_chrom \
    --interval_start $interval_start \
    --interval_end $interval_end \
    --peak_pos $peak_marker \
    --vcf $vcf \
    --pheno_1 $pheno_1 \
    --pheno_2 $pheno_2 \
    --trait_id $trait_id \
    --out 20231212_Silver.nitrate.7.8_Methylmercury
```

** MAKE sure 

## To - Do list 
- [X] Write python script to find common strains between two traits
- [X] Env for python script to find common strains. 
    - [X] create env on quest
    - [X] specify on the process definition
- [X] Add the chromosome rename key to the bin directory
- [X] Add N_snps to the LD_Process
- [X] Publish to github
- [ ] Add code to publish the LD files to the output directory
- [ ] Check if the common strains code is working properly
- [X] Version control for bcftools used
- [X] Fix the atrocity that is the bcftools code
  - [ ] figure out what the awk code is doing 
- [ ] 

[ ] - put copy of the input file in the NF output directory
# Extract LD between overlapping peak markers from `ld_overlap.nf` output with `pull_overlap_LD.R`

The script `code/qtl_overlaps/ld_between_peak_markers/pull_overlap_LD.R` is used to extract the output from the `ld_overlap.nf` pipeline output folder stored in `data/processed/qtl_overlaps/ld_between_peak_markers/`.The script will extract the LD between the peak markers for the overlapping intervals and write it to a file. 

## Inputs
- relative path to the output directory of the `ld_overlap.nf` pipeline
- nf pipeline input file `toxin_length_finemapping_input_with_pheno_files.tsv`

## Outputs
- LD between the peak markers for the overlapping intervals 1`data/processed/qtl_overlaps/ld_between_peak_markers/{today}_INBRED_QTL_overlap_peak_LD.tsv`