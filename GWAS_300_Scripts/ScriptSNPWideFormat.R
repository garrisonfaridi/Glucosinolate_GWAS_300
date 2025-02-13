# SNP WIDE FORMAT TABLES

# ALIPHATIC INDIVIDUAL GLUC GWAS WIDE FORMAT
library(tidyverse)

# File paths and corresponding names from Aliphatic_List
file_paths <- list(
  "output_300/300_gluc1_maf0.03_miss0.1GF.assoc.txt", #1
  "output_300/300_gluc2_maf0.03_miss0.1GF.assoc.txt", #2
  "output_300/300_gluc3_maf0.03_miss0.1GF.assoc.txt", #3
  "output_300/300_gluc4_maf0.03_miss0.1GF.assoc.txt", #4
  "output_300/300_gluc5_maf0.03_miss0.1GF.assoc.txt", #5
  "output_300/300_gluc6_maf0.03_miss0.1GF.assoc.txt", #6
  "output_300/300_gluc7_maf0.03_miss0.1GF.assoc.txt", #7
  "output_300/300_gluc8_maf0.03_miss0.1GF.assoc.txt", #8
  "output_300/300_gluc9_maf0.03_miss0.1GF.assoc.txt", #9
  "output_300/300_gluc10_maf0.03_miss0.1GF.assoc.txt", #10
  "output_300/300_gluc11_maf0.03_miss0.1GF.assoc.txt", #11
  "output_300/300_gluc12_maf0.03_miss0.1GF.assoc.txt", #12
  "output_300/300_gluc13_maf0.03_miss0.1GF.assoc.txt" #13
)

Aliphatic_List <- list(
  "7mSh",
  "8MTO",
  "3mSOp",
  "4mSOb",
  "5mSOp",
  "6mSOh",
  "7mSOh",
  "8mSOo",
  "Pren",
  "Buen",
  "Peen",
  "S2hBuen",
  "2hPeen"
)

# Initialize an empty data frame to hold the merged data
gwas_all_ali <- NULL

# Loop through file paths and aliphatic list
for (i in seq_along(file_paths)) {
  # Read the current file
  current_file <- read.delim(file_paths[[i]], sep = "\t", stringsAsFactors = FALSE)
  
  # Create the rs column by combining chr and ps
  current_file$rs <- paste(current_file$chr, current_file$ps, sep = "_")
  
  # Keep only the required columns
  current_file <- current_file[, c("rs", "chr", "ps", "af", "p_wald")]
  
  # Rename the p_wald column to the corresponding name in Aliphatic_List
  colnames(current_file)[5] <- Aliphatic_List[[i]]
  
  # Merge with the main data frame
  if (is.null(gwas_all_ali)) {
    # First file initializes the merged data
    gwas_all_ali <- current_file
  } else {
    # Subsequent files are merged by the rs column
    gwas_all_ali <- merge(gwas_all_ali, current_file, by = c("rs", "chr", "ps", "af"), all = TRUE)
  }
}

# Add min p value column
gwas_all_ali$p_min <- do.call(pmin, gwas_all_ali[, 5:ncol(gwas_all_ali)])

# Write the final merged data to a CSV file
write.csv(gwas_all_ali, "output_300/snp_table_all_aliphatic_pvals.csv", row.names = FALSE)

print("Merged CSV has been created at output_300/snp_table_all_aliphatic_pvals.csv")


# INDOLIC INDIVIDUAL GLUC GWAS WIDE FORMAT

# Create a list of Indolic file paths
file_paths <- list(
  "output_300/300_gluc14_maf0.03_miss0.1GF.assoc.txt",
  "output_300/300_gluc15_maf0.03_miss0.1GF.assoc.txt",
  "output_300/300_gluc16_maf0.03_miss0.1GF.assoc.txt",
  "output_300/300_gluc17_maf0.03_miss0.1GF.assoc.txt"
)

Indolic_List <- list(
  "IM",
  "1moIM",
  "1hIM",
  "4moIM"
)

# Initialize an empty data frame to hold the merged data
gwas_all_ind <- NULL

# Loop through file paths and aliphatic list
for (i in seq_along(file_paths)) {
  # Read the current file
  current_file <- read.delim(file_paths[[i]], sep = "\t", stringsAsFactors = FALSE)
  
  # Create the rs column by combining chr and ps
  current_file$rs <- paste(current_file$chr, current_file$ps, sep = "_")
  
  # Keep only the required columns
  current_file <- current_file[, c("rs", "chr", "ps", "af", "p_wald")]
  
  # Rename the p_wald column to the corresponding name in Aliphatic_List
  colnames(current_file)[5] <- Indolic_List[[i]]
  
  # Merge with the main data frame
  if (is.null(gwas_all_ind)) {
    # First file initializes the merged data
    gwas_all_ind <- current_file
  } else {
    # Subsequent files are merged by the rs column
    gwas_all_ind <- merge(gwas_all_ind, current_file, by = c("rs", "chr", "ps", "af"), all = TRUE)
  }
}

# Add min p value column
gwas_all_ind$p_min <- do.call(pmin, gwas_all_ind[, 5:ncol(gwas_all_ind)])

# Write the final merged data to a CSV file
write.csv(gwas_all_ind, "output_300/snp_table_all_indolic_pvals.csv", row.names = FALSE)

print("Merged CSV has been created at output_300/snp_table_all_indolic_pvals.csv")

