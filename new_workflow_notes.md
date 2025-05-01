
The nextflow script `peak_marker_ld.nf` is used to calculate the LD between a pair of peak markers from overlapping intervals

Mandatory inputs to the pipeline:
- `params.vcf` - path to the vcf file on RF (path to index file is generated so a .csi index is also present)
- `strains` - Input to bcftools `-S` File of sample names to include or exclude if prefixed with "^". One sample per line. 
- `peak_a` - peakid a (e.g. "II:10203")
- `peak_b` - peakid b

# `filter_vcf`

First is to filter a vcf file specified by the vcf parameter of the pipeline

The process will run the commnd which will filter the VCF file to the overlapping region and to include data for only the strains provided in the list supplied to the `-S` argument.
```
bcftools view  -S ${common_strains} -Ou ${vcf} |\\
bcftools filter -i N_MISSING=0 -o finemap.vcf.gz   
```
We will use the bcftools container used in the NemaScan simulation pipeline. We will define the container for the process as well as the execution parameters for RF in a config file `conf/rockfish.config` within the `process` definition

process{
    executor = "slurm"
    clusterOptions = '-A eande106 -e errlog.txt -N 1'
    time = "1.hour"
    cpus = 1
    memory = "4G"
    partition = "parallel"

    withLabel: "bcftools_filter_vcf" {
        container = "docker://quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"
        time = "3.hour"
        cpus = 4
        memory = "40G"
    }
}

# `calculate_ld`

Second we take that VCF file and use plink to calculate the LD

The process will take the following inputs:
file(vcf)
val(peak_a) val(peak_b)

Process will output a log file with the ld values for the pairs of SNPs which will be saved to the output directory.

Calculate the r2 for a pair of SNPs in the VCF file.
```
plink --vcf ${vcf} \\
    --threads 5 \\
    --snps-only \\
    --maf 0.05 \\
    --biallelic-only \\
    --allow-extra-chr \\
    --set-missing-var-ids @:# \\
    --ld ${peak_a} ${peak_b} \\
    --out ${peak_a}_${peak_b}
```

Plink parameter notes:
- Calculate LD between a pair of SNPs
  --ld rs2840528 rs7545940

The process will use this container definiton:
```
withLabel: "plink_recode_vcf" {
    container = "docker://andersenlab/plink:1.9"
    time = "3.hour"
    cpus = 4
    memory = "40G"
}
```