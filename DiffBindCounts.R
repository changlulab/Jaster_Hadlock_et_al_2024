####setup
library(DiffBind)
library(dplyr)
library(GenomicRanges)
library(ChIPpeakAnno)
library(ChIPseeker)

# # See DiffBind documentation for how to put together sample sheet

setwd("C:/ChIP_seq")

infiles = "C:/ChIP_seq/infiles"
outfiles = "C:/ChIP_seq/outfiles"

gene_link = "C:/ChIP_seq/mm10Promoter.bed"

# 
groups=c('V/V_F','V/P_F','O/V_F','O/P_F', 'V/V_M','V/P_M','O/V_M','O/P_M')
fdr=0.05
run_fcs=c(1) #if you want to check multiple fold-changes, you can add them here

data <- dba(sampleSheet = "C:/ChIP_seq/PSI_Op_DiffAn_Samplesheet_MaleFemaleNAC_Dup95.csv")
data
data_blacklist_remove = dba.blacklist(data, blacklist = DBA_BLACKLIST_MM10, greylist = TRUE)
gl = read.table(gene_link, header = FALSE, stringsAsFactors = FALSE, sep="\t", col.names=c("chr","start","stop"))
gl_gr = makeGRangesFromDataFrame(gl, keep.extra.columns = TRUE)

# Consensus should be the column that has the groupings that you are comparing
rep_consensus_set <- dba.peakset(data_blacklist_remove, consensus = DBA_TREATMENT, minOverlap = 2)
rep_consensus_set
rep_consensus <- dba(rep_consensus_set, mask=rep_consensus_set$masks$Consensus, minOverlap=1)
rep_consensus
ConsensusPeaks <- dba.peakset(rep_consensus, bRetrieve=TRUE)
ConsensusPeaks

#Generate counts matrix
data_analysis <- dba.count(data_blacklist_remove,summits = FALSE, peaks = ConsensusPeaks,filter=1,
                           bScaleControl = TMM_MINUS_FULL)

counts <- dba.peakset(data_analysis, bRetrieve=T, DataType=DBA_DATA_FRAME)
write.csv(counts, 'counts_NAC_MF.csv')

