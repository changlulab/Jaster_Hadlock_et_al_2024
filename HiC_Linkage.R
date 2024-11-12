############preparation##############
library(dplyr)
library(ggplot2)


setwd("C:/ChIP_seq")
# Read the two input files
file1 <- read.table("Differential_B1A1_NACMF_out.bed", header = FALSE)
file2 <- read.csv("Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}

# Write the matched_rows data frame to a BED file
write.table(matched_rows, "B1A1_NAC_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "B1A1_NAC_non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")