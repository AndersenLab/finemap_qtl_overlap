library(tidyverse)

source("bin/outs.R")

#### Inputs ####
# Directory containing trait pair folders with LD data
tp_folders_dir <- "data/processed/qtl_overlaps/ld_between_peak_markers/Analysis_Results-20240111"

# Input file containing overlap information
overlap_input_file <- "data/processed/qtl_overlaps/ld_between_peak_markers/Analysis_Results-20240111/toxin_length_finemapping_input_with_pheno_files.tsv"

#### Outputs ####
# Compressed TSV file containing LD between peak markers
output_tsv <- "data/processed/qtl_overlaps/inbred_qtl_overlap_peak_LD.tsv.gz"

#### Functions #### ----------
# Function to pull LD for an overlap interval
pull_ld <- function(overlap_int_df, tp_folder){
    # Get the name of the trait pair from the folder
    tp_name <- basename(tp_folder)

    # Get the path to the LD files using the trait pair folder and data from the overlaps_df
    overlap_ld_file_A <- glue::glue(
        "{tp_folder}/peakPOS_A/Data/{tp_name}.{overlap_int_df$CHROM_A}.{overlap_int_df$leftmost}.{overlap_int_df$rightmost}.LD_inbred.tsv"
    )

    overlap_ld_file_B <- glue::glue(
        "{tp_folder}/peakPOS_B/Data/{tp_name}.{overlap_int_df$CHROM_A}.{overlap_int_df$leftmost}.{overlap_int_df$rightmost}.LD_inbred.tsv"
    )

    # Check if the files exist
    if(file.exists(overlap_ld_file_A) & file.exists(overlap_ld_file_B)){
        print("Files exist")
        
        ld_data_A <- data.table::fread(overlap_ld_file_A)

        ld_data_B <- data.table::fread(overlap_ld_file_B)

        # combine the two data sets
        interval_ld <- bind_rows(ld_data_A, ld_data_B)

        # Get the peak ids from the overlap_int_df
        peak_ids <- overlap_int_df$peak_pos %>% strsplit(",") %>% unlist()

        peak_1 <- peak_ids[1]
        peak_2 <- peak_ids[2]

        peak_ld_intervals <- interval_ld %>% 
            filter(
                (BP_A == peak_1 & BP_B == peak_2) | 
                (BP_B == peak_2 & BP_A == peak_1)
            )
        
        # add the LD (R2) to the overlap_int_df
        peak_ld_intervals <- peak_ld_intervals %>% 
            mutate(
                trait_pair = tp_name,
                CHROM_A = overlap_int_df$CHROM_A,
                leftmost = overlap_int_df$leftmost,
                rightmost = overlap_int_df$rightmost
            )
            
        return(peak_ld_intervals)

    }
}

# Function to pull LD for peak markers for a trait pair 
tp_ld <- function(tp_folder, overlaps_df){
    # Get the name of the trait pair from the folder
    tp_name <- basename(tp_folder)

    # Filter the overlaps_df for the trait pair
    tp_overlaps_df <- overlaps_df %>% 
        filter(trait_pair == tp_name)
    
    # Print the number of overlaps for the trait pair
    n_overlaps <- nrow(tp_overlaps_df)
    cat("Number of overlaps:", n_overlaps, "\n")

    # Get the unique overlap intervals for the trait pair
    tp_overlaps_int <- tp_overlaps_df %>% 
        dplyr::distinct() %>%
        dplyr::group_by(CHROM_A, leftmost, rightmost) %>%
        dplyr::summarise(peak_pos = paste(peak_pos, collapse = ",")) %>% 
        dplyr::ungroup()

    # Get trait A and trait B names
    trait_A <- tp_overlaps_df$traitA[1]
    trait_B <- tp_overlaps_df$traitB[1]

    # Pull the LD for each overlap interval by iterating 
    # pull_ld function over the rows of tp_overlaps_int
    tp_overlaps_ld <- map_dfr(1:nrow(tp_overlaps_int), ~{
        overlap_int_df <- tp_overlaps_int[.x, ]
        pull_ld(overlap_int_df, tp_folder)
    })

    # Add the Trait A and Trait B IDs to the data frame
    tp_overlaps_ld <- tp_overlaps_ld %>% 
        mutate(
            traitA = trait_A,
            traitB = trait_B
        )

    return(tp_overlaps_ld)

} 

#### Testing #### ----------
# tp_folder <- "data/processed/qtl_overlaps/ld_between_peak_markers/Analysis_Results-20240111/length_Chlorothalonil.length_Arsenic_trioxide"

# overlap_input_file <- "code/qtl_overlaps/ld_between_peak_markers/toxin_length_finemapping_input_with_pheno_files.tsv"
# overlaps <- data.table::fread(overlap_input_file)

# tp_ld(tp_folder, overlaps) %>% 
#     write_tsv("test_peak_LD.tsv")


#### Run for all trait pairs #### ----------

# Get the list of trait pair folders
tp_folders <- list.dirs(tp_folders_dir, full.names = TRUE, recursive = FALSE)

## return the number of trait pair folders
cat("Number of trait pair folders:", length(tp_folders), "\n")

# Read the input file
overlaps <- data.table::fread(overlap_input_file)

# Run the tp_ld function for each trait pair folder and return a list of data frames
tp_ld_list <- map(tp_folders, ~tp_ld(.x, overlaps))

print("Done")

# Combine the list of data frames into a single data frame
all_peak_ld <- bind_rows(tp_ld_list)

# Write the data frame to a file
overlap_peak_ld <- all_peak_ld %>% 
    # create a peakidA and peakidB column from SNP_A and SNP_B
    mutate(
        peakidA = paste(CHROM_A, BP_A, sep = ":"),
        peakidB = paste(CHROM_A, BP_B, sep = ":")
    ) %>%
    select(traitA, traitB, peakidA, peakidB, R2)  %>% 
    # order by highest R2
    arrange(desc(R2))%>%
    # add column to identify trait pairs that are outliers > 1 SD
    mutate(
        outlier = ifelse(R2 > mean(R2) + sd(R2), TRUE, FALSE)
    )
    
    


save_compressed_tsv(
    overlap_peak_ld,
    output_tsv
)


# #### Upload the supplementary data table to gsheet ####
# # Link to the supplementary data spreadsheet
# tl <- "https://docs.google.com/spreadsheets/d/1qJbh2a1JgOGAKwWlc6WOSMuG-aVh-vI_rCNDYxnIL20/edit?usp=sharing"

# # # Authorize sheet access
# # gs4_auth(
# # )

# # write the as a new sheet in the spreadsheet
# googlesheets4::sheet_write(
#     data = overlap_peak_ld,
#     ss = tl,
#     sheet = "Inbred QTL Overlap Peak LD"
# )
