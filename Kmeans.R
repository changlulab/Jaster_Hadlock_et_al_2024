############preparation##############
library(ggplot2)
library(dplyr)
library(stats)
library(gplots)
library(tidyverse)
library(pheatmap)

###K-Means

my_data1 <- read.csv("NACMF_male_3comp_6col_Grouped_input_norm.csv")
rownames(my_data1) <- my_data1$start

# Select the column you want to use for clustering
my_data <- my_data1[, 2:4]  # Assuming you want to use the average column for clustering

# Perform k-means clustering, number of clusters determined by silhouette
k <- 4
set.seed(123)
kmeans_result <- kmeans(my_data, centers = k, nstart = 25)

# Get cluster assignments
cluster_labels <- kmeans_result$cluster

# Initialize list to store genes in each cluster
cluster_gene_lists <- vector("list", length = k)

# Fill the list with genes in each cluster
for (cluster in 1:k) {
  genes_in_cluster <- rownames(my_data1)[cluster_labels == cluster]
  cluster_gene_lists[[cluster]] <- genes_in_cluster
}

# Export gene lists to CSV
for (cluster in 1:k) {
  cluster_filename <- paste("cluster_", cluster, "_genes_PFCMF_6col_grouped_norm4.csv", sep = "")
  write.csv(cluster_gene_lists[[cluster]], file = cluster_filename, row.names = FALSE)
}


my_data2 <- read.csv("NACMF_male_3comp_6col_Grouped_input_norm_2.csv")
my_data2$Cluster <- cluster_labels
MDT <- my_data2[,-3:-7]
MDT1 <- MDT[,-1]
write.csv(MDT1, "NACMF_male_3comp_6col_Grouped_Clusters.csv", row.names = FALSE)

###Reorder clusters to 1, 2, 3, 4
data <- read_csv("NACMF_male_3comp_6col_Grouped_Clusters.csv", col_names = TRUE)

# Map values in the second column according to your transformation
data$Cluster <- case_when(
  data$Cluster == 2 ~ 1,
  data$Cluster == 1 ~ 2,
  data$Cluster == 4 ~ 3,
  data$Cluster == 3 ~ 4,
  TRUE ~ data$Cluster  # handle any other values if necessary
)

# Write the modified data back to a CSV file
write_csv(data, "modified_PFCMF_male_3comp_6col_Grouped_Clusters.csv", col_names = FALSE)









# Load necessary libraries
library(dplyr)

# File paths of your two CSV files
file1_path <- "modified_PFCMF_male_3comp_6col_Grouped_Clusters.csv"
file2_path <- "modified_PFCMF_female_3comp_6col_Grouped_Clusters.csv"

# Read the first CSV file
data1 <- read.csv(file1_path, header = FALSE, stringsAsFactors = FALSE)
colnames(data1) <- c("Column1", "Column2")

# Read the second CSV file
data2 <- read.csv(file2_path, header = FALSE, stringsAsFactors = FALSE)
colnames(data2) <- c("Column1", "Column2")

# Combine the two datasets
combined_data <- bind_rows(data1, data2)

# Remove duplicates based on Column1
combined_data_unique <- combined_data %>%
  distinct(Column1, .keep_all = TRUE)

# Write the combined and de-duplicated data to a new CSV file
write.csv(combined_data_unique, "combined_PFCMF_3comp_6col_Grouped_Cluster_Numbers.csv", row.names = FALSE)





# Load necessary libraries (step6)
library(dplyr)

# Read file 1 (assuming 8 columns)
file1 <- read.csv("PFCMF_3comp_6col_Grouped_input_norm.csv", header = TRUE)

# Read file 2 (assuming 2 columns)
file2 <- read.csv("combined_PFCMF_3comp_6col_Grouped_Cluster_Numbers.csv", header = TRUE)

# Rename columns for clarity
colnames(file1) <- paste0("V", 1:8)
colnames(file2) <- c("Column1", "Column2")

# Perform left join to add Column2 values from file2 to file1
combined_data <- left_join(file1, file2, by = c("V1" = "Column1"))

# Mutate to update V9 based on matching values
combined_data <- combined_data %>%
  mutate(V9 = ifelse(!is.na(Column2), Column2, V9)) %>%
  select(-Column1, -Column2)


# Write the combined data with updated column 9 to a new CSV file
write.csv(combined_data, "combined_PFCMF_3comp_6col_Grouped_input.csv", row.names = FALSE)

my_data2 <- read.csv("combined_PFCMF_3comp_6col_Grouped_input.csv")
colnames(my_data2) <- c("start", "VEH-VEH, Male","OXY-VEH, Male", "OXY-PSI, Male", "VEH-VEH, Female","OXY-VEH, Female", "OXY-PSI, Female", "Average", "Cluster")
custom_order <- c(1, 2, 3, 4)

# Use match function to create a temporary ordering column based on custom_order
my_data2$order <- match(my_data2$Cluster, custom_order)

# Sort my_data2 by the new order column
my_data2 <- my_data2[order(my_data2$order), ]

# Remove the temporary order column
my_data2$order <- NULL


# rownames(my_data2) <- my_data2$START

mydata3 <- my_data2[, 2:7]


# Convert to matrix for heatmap plotting
mydata3 <- as.matrix(mydata3)


library(pheatmap)

phtmap <- pheatmap::pheatmap(mydata3)
unique_labels <- unique(my_data2$Cluster)
labels_row <- rep("", nrow(mydata3))
for (i in 1:length(unique_labels)) {
  indices <- which(my_data2$Cluster == unique_labels[i])
  labels_row[indices[1]] <- paste("Cluster", as.character(unique_labels[i]))
}
breaks <- seq(0, 0.3, by = 0.001172)
# Plot heatmap with custom column clustering

figure <- pheatmap(mydata3, 
                   scale = "none",
                   breaks = breaks,
                   cluster_rows = FALSE,       # Turn off row clustering
                   cluster_cols = FALSE,
                   #labels_row = labels_row,
                   border_color = NA,          # No border
                   fontsize_row = 8,           # Adjust row label font size
                   fontsize = 8,               # Adjust font size for cell values
                   color = colorRampPalette(c("white", "#ADD8E6" ,"#000080" ))(256))  # Color palette