library(tidyverse)

output_directory <- "data/processed/qtl_overlaps/ld_between_peak_markers/Analysis_Results-20240111"

#list the folders in the output directory
pair_folders <- dir(output_directory, full.names = TRUE)
head(pair_folders)
#remove the .html files
pair_folders <- pair_folders[!grepl(".html", pair_folders)]

# Get the name of the trait pairs with QTL overlap analysis data from the folder names
pair_names <- basename(pair_folders)


# Load the input file to the nextflow pipeline
overlap_input_file <- "code/qtl_overlaps/ld_between_peak_markers/toxin_length_finemapping_input_with_pheno_files.tsv"

overlaps <- data.table::fread(overlap_input_file)

n_overlaps <- nrow(overlaps)

print(
  glue::glue("There are {n_overlaps} overlap entries in the input file.")
)

# For all trait pair with a least one QTL overlap the overlaping QTL intervals 
# are stored as two entries in the overlaps dataframe.

n_trait_paris <- overlaps %>% 
  distinct(trait_pair) %>% 
  nrow()

print(
  glue::glue("There are {n_trait_paris} trait pairs with at least one QTL overlap.")
)
# Reduce the overlaps data to just inputs 

overlaps_intervals_df <- overlaps  %>% 
    select(trait_pair, CHROM_A, leftmost, rightmost, peak, peak_pos)  %>% 
    #pivot the peak column (peakPOS_A or peak_POS_B) to columns with peak_pos as the value
    pivot_wider(names_from = peak, values_from = peak_pos) 

# Write a function to pull ld for each row 

get_overlap_ld <- function(output_directory, overlap_df){
    # Input path to the NF output data and an overlap intervals df filtered to a single overlap
    # get the LD file for the trait pair folder
    ## Pull the relevant info 
    trait_pair_id <- overlap_df$trait_pair

    chrom_id <- overlap_df$CHROM_A

    leftmost <- overlap_df$leftmost

    rightmost <- overlap_df$rightmost

    peakA_id <- overlap_df$peakPOS_A

    peakB_id <- overlap_df$peakPOS_B
    
    ld_file_name <- glue::glue("{output_directory}/{trait_pair_id}/peakPOS_A/Data/{trait_pair_id}.{chrom_id}.{leftmost}.{rightmost}.LD_inbred.tsv")

    #Load the LD file for the trait pair folder
    ld_df <- data.table::fread(ld_file_name)

    print(head(ld_df))
}

# Test the function on the first row of the overlaps df
get_overlap_ld(output_directory, overlaps_intervals_df[1,])



get_peak_marker_ld <- function(trait_pair_folder, overlaps_df) {
  # function that acts on a trait_pair folder and
  # returns the LD between the peak markers for each trait pair.
  # The function will use the `overlaps` data.frame to get the peak marker IDs for each trait.
  # Then the function will read in the LD files stored in the trait_pair folder in the <trait_pair>/peakPOS_A/Data/_*

  # Get the basename of the trait pair folder
  pair_name <- basename(trait_pair_folder)

  # Get the peak marker IDs for each trait
  peak_marker_A <- overlaps_df %>%
    filter(trait_pair == pair_name) %>%
    filter(peak == "peakPOS_A") %>%
    select(CHROM_A, peak_pos) %>%
    # mutate(peak_marker = paste0(CHROM_A, ":", peak_pos)) %>%
    # pull(peak_marker)%>%
    pull(peak_pos)

  print(peak_marker_A)

  peak_marker_B <- overlaps_df %>%
    filter(trait_pair == pair_name) %>%
    filter(peak == "peakPOS_B") %>%
    select(CHROM_A, peak_pos) %>%
    # mutate(peak_marker = paste0(CHROM_A, ":", peak_pos)) %>%
    # pull(peak_marker)
    pull(peak_pos)

  print(peak_marker_B)

  if (length(peak_marker_A) > 1 && length(peak_marker_B) > 1) {
    print("More than one peak marker found for trait pair")

    # find all unique overlap intervals
    for (peak_b in peak_marker_B) {
      ld_dir <- glue::glue("{trait_pair_folder}/peakPOS_A/Data")
      # create a list of all *.LD_inbred.tsv files in the dir
      ld_files <- list.files(ld_dir, pattern = "*.LD_inbred.tsv", full.names = TRUE)

      # Load the LD file for the trait pair folder
      trait_pair_overlaps <- data.frame()
      # iterate over the LD files
      for (i in 1:length(ld_files)) {
        ld_df <- data.table::fread(ld_files[1])

        # Filter the LD file for the peak marker by looking for the the `peak_marker_B` in the `BP_B` column
        ld_df <- ld_df %>%
          filter(BP_B == peak_marker_B)

        # Create a dataframe output with the columns:
        #-`trait_pair`
        #- `SNP_A`
        #- `MAF_A`
        #- `SNP_B`
        #- `MAF_B`
        #- `R2
        output_df <- data.frame(
          trait_pair = pair_name,
          SNP_A = ld_df$SNP_A,
          MAF_A = ld_df$MAF_A,
          SNP_B = ld_df$SNP_B,
          MAF_B = ld_df$MAF_B,
          R2 = ld_df$R2
        )

        trait_pair_overlaps <- rbind(trait_pair_overlaps, output_df)
      }
    }
  } else {



    # Get the LD file for the trait pair folder
    ld_dir <- glue::glue("{trait_pair_folder}/peakPOS_A/Data")
  }
  # create a list of all *.LD_inbred.tsv files in the dir
  ld_files <- list.files(ld_dir, pattern = "*.LD_inbred.tsv", full.names = TRUE)

  # Load the LD file for the trait pair folder
  trait_pair_overlaps <- data.frame()
  # iterate over the LD files
  for (i in 1:length(ld_files)) {
    ld_df <- data.table::fread(ld_files[1])

    # Filter the LD file for the peak marker by looking for the the `peak_marker_B` in the `BP_B` column
    ld_df <- ld_df %>%
      filter(BP_B == peak_marker_B)

    # Create a dataframe output with the columns:
    #-`trait_pair`
    #- `SNP_A`
    #- `MAF_A`
    #- `SNP_B`
    #- `MAF_B`
    #- `R2
    output_df <- data.frame(
      trait_pair = pair_name,
      SNP_A = ld_df$SNP_A,
      MAF_A = ld_df$MAF_A,
      SNP_B = ld_df$SNP_B,
      MAF_B = ld_df$MAF_B,
      R2 = ld_df$R2
    )

    trait_pair_overlaps <- rbind(trait_pair_overlaps, output_df)
  }
  return(trait_pair_overlaps)
}


#Now lets test the function on the first trait pair folder

test_file <- pair_folders[2]
print(basename(test_file))

test_output <- get_peak_marker_ld(test_file, overlaps)
head(test_output)

## Checking error with file number 4
test_file <- pair_folders[4]
print(test_file)

test_output <- get_peak_marker_ld(test_file, overlaps)
head(test_output)


# Now lets run the function on all trait pair folders and output the results to a single dataframe

ld_df <- data.frame()
for(i in 1:length(pair_folders)){
  print(i)
  print(pair_folders[i])
  ld_df <- rbind(ld_df, get_peak_marker_ld(pair_folders[i], overlaps))
}

head(ld_df)
