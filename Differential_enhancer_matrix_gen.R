B1A1_NAC_combined_Annotated.bed

############preparation##############

library(ggplot2)
library(dplyr)
library(stats)
library(gplots)
library(tidyverse)
library(pheatmap)


setwd("C:/data/LSD_Op/ChIP_seq")
dataframe_A <- read.csv("counts_NAC_MF.csv")

bed_file1 <- read.table("C1A1_NAC_combined_Annotated.bed", header = FALSE, stringsAsFactors = FALSE)
bed_file2 <- read.table("D1A1_NAC_combined_Annotated.bed", header = FALSE, stringsAsFactors = FALSE)
bed_file3 <- read.table("C1D1_NAC_combined_Annotated.bed", header = FALSE, stringsAsFactors = FALSE)
bed_file4 <- read.table("G1E1_NAC_combined_Annotated.bed", header = FALSE, stringsAsFactors = FALSE)
bed_file5 <- read.table("H1E1_NAC_combined_Annotated.bed", header = FALSE, stringsAsFactors = FALSE)
bed_file6 <- read.table("G1H1_NAC_combined_Annotated.bed", header = FALSE, stringsAsFactors = FALSE)

# Combine the six data frames into one
dataframe_B <- rbind(bed_file1, bed_file2, bed_file3, bed_file4, bed_file5, bed_file6)

# Get unique values from column 2 of dataframe B
values_to_match <- unique(dataframe_B[, 2])

# Filter rows from dataframe A where column 2 matches any value from dataframe B
matching_rows <- dataframe_A[dataframe_A[, 2] %in% values_to_match, ]

write.csv(matching_rows, "NACMF_3comp_DP_Count.csv", row.names = FALSE)

#check for duplicate start values
# Check for duplicates in "START" column
duplicate_rows <- matching_rows[duplicated(matching_rows$START) | duplicated(matching_rows$START, fromLast = TRUE), ]

# Print rows with duplicates in "START" column
if (nrow(duplicate_rows) > 0) {
  cat("Duplicate values in START column:\n")
  print(duplicate_rows)
} else {
  cat("No duplicate values found in START column.\n")
}
##May have to remove duplicates manually, raise start position by 1 then save and read-in file

matching_rows <- read.csv("NACMF_3comp_DP_Count.csv")
# Create a new dataframe with columns from the 5th column onwards
new_dataframe <- matching_rows %>% select(5:ncol(matching_rows))

# Normalize the data using L2 norm
normalize_l2 <- function(x) {
  norm <- sqrt(sum(x^2))
  return(x / norm)
}

# Apply the normalization function to each row
normalized_rows <- t(apply(new_dataframe, 1, normalize_l2))

# Convert the result back to a dataframe
normalized_dataframe <- as.data.frame(normalized_rows)

# Extract columns containing Group.A, Group.B, Group.C, and Group.D
GroupA_columns <- normalized_dataframe %>%
  select(matches("A."))

GroupC_columns <- normalized_dataframe %>%
  select(matches("C."))

GroupD_columns <- normalized_dataframe %>%
  select(matches("D."))

GroupE_columns <- normalized_dataframe %>%
  select(matches("E."))

GroupG_columns <- normalized_dataframe %>%
  select(matches("G."))

GroupH_columns <- normalized_dataframe %>%
  select(matches("H."))


# Calculate the mean for each Group
GroupA_mean <- rowMeans(GroupA_columns, na.rm = TRUE)
GroupC_mean <- rowMeans(GroupC_columns, na.rm = TRUE)
GroupD_mean <- rowMeans(GroupD_columns, na.rm = TRUE)
GroupE_mean <- rowMeans(GroupE_columns, na.rm = TRUE)
GroupG_mean <- rowMeans(GroupG_columns, na.rm = TRUE)
GroupH_mean <- rowMeans(GroupH_columns, na.rm = TRUE)


row_means <- rowMeans(cbind(GroupA_mean, GroupC_mean, GroupD_mean, GroupE_mean, GroupG_mean, GroupH_mean), na.rm = TRUE)

# Create the new dataframe
input_1 <- data.frame(
  start = matching_rows[,2],
  GroupA = GroupA_mean,
  GroupC = GroupC_mean,
  GroupD = GroupD_mean,
  GroupE = GroupE_mean,
  GroupG = GroupG_mean,
  GroupH = GroupH_mean,
  Average = row_means # Column name changed to "Average"
)

input_2 <- data.frame(
  CHR = matching_rows[,1],
  START = matching_rows[,2],
  END = matching_rows[,3],
  GroupA = GroupA_mean,
  GroupC = GroupC_mean,
  GroupD = GroupD_mean,
  GroupE = GroupE_mean,
  GroupG = GroupG_mean,
  GroupH = GroupH_mean,
  Average = row_means # Column name changed to "Average"
)
input_3 <- data.frame(START = matching_rows[,2], Average = row_means, normalized_dataframe)
input_4 <- data.frame(matching_rows[,1:3], Average = row_means, normalized_dataframe)
# Write the new dataframe to a CSV file
write.csv(input_1, "NACMF_3comp_6col_Grouped_input_norm.csv", row.names = FALSE)
write.csv(input_2, "NACMF_3comp_6col_Grouped_input_norm_2.csv", row.names = FALSE)
write.csv(input_3, "NACMF_3comp_6col_input_norm.csv", row.names = FALSE)
write.csv(input_4, "NACMF_3comp_6col_input_norm_2.csv", row.names = FALSE)
