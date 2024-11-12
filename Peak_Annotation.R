setwd("C:/ChIP_seq")
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

library(clusterProfiler)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)

rm(list = ls())
options(stringsAsFactors = F)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
setwd("C:/ChIP_seq")

bedPeaksFile = "B1A1_NAC_non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('There are ', length(peak), 'peaks for this data'))
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)

## Remove Promoters

# Define a function to filter rows based on column 14
remove_outliers <- function(data) {
  # Filter rows where column 14 is not within the range (-2000, 2000)
  filtered_data <- data[!(data[, 14] > -2000 & data[, 14] < 2000), ]
  return(filtered_data)
}

# Example usage:
# Assuming df is your data frame, call the function to remove outliers
# Replace df with the name of your actual data frame
new_peakAnno_df <- remove_outliers(peakAnno_df)

# Now, filtered_df contains the data frame with rows removed where the value of column 14 is within the range (-2000, 2000)

#### Only use below if you are cutting post HiC comp non-matched peaks
new_peakAnno_df = new_peakAnno_df[,-17]
new_peakAnno_df = new_peakAnno_df[,-4:-15]

####

##Combine the HiC & nearest neighbor enhancer annotations
HiC_Match <- read.table("B1A1_NAC_matched_rows.bed", header = FALSE)
colnames(HiC_Match)=c('seqnames', 'start','end', 'SYMBOL')
combined_bed <- rbind(new_peakAnno_df, HiC_Match)

write.table(combined_bed, "B1A1_NAC_combined_Annotated.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
