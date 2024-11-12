############Pearson Correlation 2##############

my_data2 <- read.csv("combined_NACMF_3comp_6col_Grouped_input.csv")
colnames(my_data2) <- c("start", "VEH-VEH, Male","OXY-VEH, Male", "OXY-PSI, Male", "VEH-VEH, Female","OXY-VEH, Female", "OXY-PSI, Female", "Average", "Cluster")
custom_order <- c(1, 2, 3, 4)

# Use match function to create a temporary ordering column based on custom_order
my_data2$order <- match(my_data2$Cluster, custom_order)

# Sort my_data2 by the new order column
my_data2 <- my_data2[order(my_data2$order), ]

# Remove the temporary order column
my_data2$order <- NULL

cluster1_data <- my_data2[my_data2$Cluster == 1, ]
cluster2_data <- my_data2[my_data2$Cluster == 2, ]
cluster3_data <- my_data2[my_data2$Cluster == 3, ]
cluster4_data <- my_data2[my_data2$Cluster == 4, ]
C1 <- cluster1_data[,2:7]
C2 <- cluster2_data[,2:7]
C3 <- cluster3_data[,2:7]
C4 <- cluster4_data[,2:7]


calculate_row_correlations <- function(df) {
  # Initialize an empty list to store correlations
  correlations <- numeric(nrow(df))
  
  # Loop through each row
  for (i in 1:nrow(df)) {
    # Extract the vectors for the current row
    group1 <- c(df[i, 1], df[i, 2], df[i, 3])
    group2 <- c(df[i, 4], df[i, 5], df[i, 6])
    
    # Compute the Pearson correlation coefficient
    correlation <- cor(group1, group2, method = "pearson")
    
    # Store the result in the list
    correlations[i] <- correlation
  }
  
  # Return the list of correlations
  return(correlations)
}


# Call the function with the dataframe
correlations <- calculate_row_correlations(C1)

correlations_df <- data.frame(Correlation = correlations)

# Append the new column to the original dataframe
C1_with_correlations <- cbind(C1, Correlation = correlations_df$Correlation)
# Print the result
median_correlation <- median(C1_with_correlations$Correlation, na.rm = TRUE)
average_correlation <- mean(C1_with_correlations$Correlation, na.rm = TRUE)

sorted_df <- C1_with_correlations %>%
  arrange(Correlation)

# Select the top 100 rows with the lowest values
lowest_100_rows <- head(sorted_df, 100)

# Create a new dataframe with these rows
new_dataframe <- lowest_100_rows


calculate_row_correlations2 <- function(df) {
  # Initialize an empty list to store correlations
  correlations <- numeric(nrow(df))
  
  # Loop through each row
  for (i in 1:nrow(df)) {
    # Extract the vectors for the current row
    group1 <- c(df[i, 1], df[i, 2], df[i, 3])
    group2 <- c(df[i, 4], df[i, 5], df[i, 6])
    
    # Compute the Pearson correlation coefficient
    correlation <- cor(group1, group2, method = "pearson")
    
    # Store the result in the list
    correlations[i] <- correlation
  }
  
  # Return the list of correlations
  return(correlations)
}

# Call the function with the dataframe
correlations2 <- calculate_row_correlations2(C2)

correlations_df2 <- data.frame(Correlation = correlations2)

# Append the new column to the original dataframe
C2_with_correlations <- cbind(C2, Correlation = correlations_df2$Correlation)
# Print the result
median_correlation2 <- median(C2_with_correlations$Correlation, na.rm = TRUE)
average_correlation2 <- mean(C2_with_correlations$Correlation, na.rm = TRUE)



calculate_row_correlations3 <- function(df) {
  # Initialize an empty list to store correlations
  correlations <- numeric(nrow(df))
  
  # Loop through each row
  for (i in 1:nrow(df)) {
    # Extract the vectors for the current row
    group1 <- c(df[i, 1], df[i, 2], df[i, 3])
    group2 <- c(df[i, 4], df[i, 5], df[i, 6])
    
    # Compute the Pearson correlation coefficient
    correlation <- cor(group1, group2, method = "pearson")
    
    # Store the result in the list
    correlations[i] <- correlation
  }
  
  # Return the list of correlations
  return(correlations)
}


# Call the function with the dataframe
correlations3 <- calculate_row_correlations3(C3)

correlations_df3 <- data.frame(Correlation = correlations3)

# Example dataframe C1 (Replace this with your actual dataframe)

# Append the new column to the original dataframe
C3_with_correlations <- cbind(C3, Correlation = correlations_df3$Correlation)
# Print the result
median_correlation3 <- median(C3_with_correlations$Correlation, na.rm = TRUE)
average_correlation3 <- mean(C3_with_correlations$Correlation, na.rm = TRUE)


calculate_row_correlations4 <- function(df) {
  # Initialize an empty list to store correlations
  correlations <- numeric(nrow(df))
  
  # Loop through each row
  for (i in 1:nrow(df)) {
    # Extract the vectors for the current row
    group1 <- c(df[i, 1], df[i, 2], df[i, 3])
    group2 <- c(df[i, 4], df[i, 5], df[i, 6])
    
    # Compute the Pearson correlation coefficient
    correlation <- cor(group1, group2, method = "pearson")
    
    # Store the result in the list
    correlations[i] <- correlation
  }
  
  # Return the list of correlations
  return(correlations)
}


# Call the function with the dataframe
correlations4 <- calculate_row_correlations4(C4)

correlations_df4 <- data.frame(Correlation = correlations4)

# Example dataframe C1 (Replace this with your actual dataframe)

# Append the new column to the original dataframe
C4_with_correlations <- cbind(C4, Correlation = correlations_df4$Correlation)
# Print the result
median_correlation4 <- median(C4_with_correlations$Correlation, na.rm = TRUE)
average_correlation4 <- mean(C4_with_correlations$Correlation, na.rm = TRUE)


C1_hold <-C1_with_correlations[,7]
C2_hold <-C2_with_correlations[,7]
C3_hold <-C3_with_correlations[,7]
C4_hold <-C4_with_correlations[,7]



data_list <- list(Cluster1 = C1_hold, Cluster2 = C2_hold, Cluster3 = C3_hold, Cluster4 = C4_hold)

# Create the box plot
boxplot(data_list, 
        main = "NAc Pairwise Pearson Correlation Coefficient",  # Title of the plot
        #xlab = "Cluster",                          # X-axis label
        #ylab = "PPCC",                         # Y-axis label
        col = c("#87CEEB", "#99B2FF", "#6699FF", "#000080"),  # Colors for each box
        border = "black",                        # Border color of the boxes
        notch = TRUE,                            # Add notches to the boxes
        outline = FALSE)                         # Do not show outliers
