#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)
library(data.table)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

# arguments
# 1 = raw fine mapping
# 2 - ROI geno matrix
# 3 - ld file

save_name <- gsub(glue::glue(".LD_{args[4]}.tsv"), "", args[3])
finemap <- data.table::fread(args[1]) %>%
  dplyr::mutate(CHR = as.factor(CHR),
                SNP = as.factor(SNP),
                POS = as.numeric(POS),
                AF1 = as.numeric(AF1),
                P = as.numeric(P))
ROI.geno.matrix <- data.table::fread(args[2])
ROI.LD <- data.table::fread(args[3])

peakp <- unique(ROI.LD$BP_A)

finemap_peaks <- na.omit(finemap) %>%
  dplyr::filter(POS == peakp)

finemap_peaks$POS <- as.numeric(finemap_peaks$POS)

# add strain genotype to ld file for snpeff later
roi_genotype <- ROI.geno.matrix %>%
  tidyr::gather(strain, allele, -c(CHROM:ALT)) %>%
  dplyr::mutate(allele = ifelse(allele == -1, "REF", ifelse(allele == 1, "ALT", NA))) %>%
  dplyr::group_by(CHROM, POS, REF, ALT, allele) %>%
  dplyr::summarize(strains = paste(strain, collapse = ",")) 

pr_roi_ld <- ROI.LD %>%
 # dplyr::mutate(peak_marker = gsub("_", ":", unique(finemap_peaks$SNP))) %>%
  dplyr::mutate(peak_marker = ifelse(identical(gsub("_", ":", unique(finemap_peaks$SNP)), character(0)), NA, 
                                     gsub("_", ":", unique(finemap_peaks$SNP)))) %>% 
  dplyr::mutate(SNP = SNP_B) %>%
  dplyr::select(peak_marker, peak_maf = MAF_A, SNP, maf_marker_b = MAF_B, ld_r2 = R2) %>%
  dplyr::mutate(ld_r2 = as.numeric(ld_r2)) %>%
  dplyr::filter(!is.na(ld_r2)) %>%
  dplyr::left_join(finemap,., by = "SNP") %>%
  dplyr::left_join(.,roi_genotype)

readr::write_tsv(pr_roi_ld, path = glue::glue("{save_name}.prLD_df_{args[4]}.tsv"))

peak_roi_marker <- dplyr::filter(pr_roi_ld, POS == peakp)


ld_plot <- ggplot(pr_roi_ld, mapping = aes(x = POS/1000000, y = as.numeric(-log(P)), fill = ld_r2)) +
  theme_bw(15) +
  geom_point(shape = 23, size = 3) +
  geom_point(aes(y = -log(P)), shape = 23, size = 3, fill = "red",
             data = peak_roi_marker) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = bquote(r^2)) +
  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p))))

ggsave(ld_plot, filename = glue::glue("{save_name}_finemap_plot_{args[4]}.pdf"),
       height = 4,
       width = 12)