############preparation##############
library(dplyr)
library(DESeq2)
library(umap)
library(ggplot2)

rm(list = ls())
options(stringsAsFactors = F)
setwd("C:/ChIP_seq")
factors = read.csv('cofactorsNACMF.csv')
row.names(factors) = factors[,1]
factors = factors[,-1]
factors$Peaks = round(factors$Peaks,2)
factors$NSC = round(factors$NSC,2)
factors$sequencing_depth = round(factors$sequencing_depth,2)
factors$unique_reads = round(factors$unique_reads,2)

raw_def = read.csv('counts_NAC_MF.csv')
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
def = round(def)
peak_info = raw_def[,1:4]


####strong correlation: cor>0.7################
###############################B_1 vs. A_1###############################################################
setwd("C:/ChIP_seq")
meta = bind_rows(factors[11:18,1:7],factors[1:10,1:7])
# might need to change below so the second numbers are equivalent
meta$condition_num = c(rep(1,8),rep(2,10))
group_B1A1 = data.frame(def[,11:18],def[,1:10])
miss_1 <- c()
for (i in 1:nrow(group_B1A1)) {
  if(length(which(group_B1A1[i,1:10]<10)) > 0.5*ncol(group_B1A1[,1:10])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_B1A1)) {
  if(length(which(group_B1A1[i,11:18]<10)) > 0.5*ncol(group_B1A1[,11:18])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_B1A1 = group_B1A1[-miss,]
t_group_B1A1 = t(group_B1A1)
umap_results <- umap(t_group_B1A1,n_neighbors=8,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')

corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')
# 
corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')
# 
corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')
# 
corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')
# 
corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')
# 
corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')
# 
group_list1 = meta$condition
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC

group_list5 = meta$FRiP
group_list6 = meta$sequencing_depth
group_list7 = meta$batch
# 
colData=data.frame(row.names = colnames(group_B1A1),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   FRiP = group_list5,
                   sequencing_depth = group_list6,
                   batch = group_list7)

colData$condition = factor(colData$condition,levels = c('B','A'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=factor(colData$batch,)
# 
dds=DESeqDataSetFromMatrix(countData = group_B1A1,
                           colData = colData,
                           design = ~ condition)
# #tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)
# 
res_NACM=results(dds,
                 contrast = c("condition","B","A"))
resOrdered_NACM=res_NACM[order(res_NACM$padj),]
DEG_NACM=as.data.frame(resOrdered_NACM)
DEG_NACM=na.omit(DEG_NACM)
nrDEG_NACM=DEG_NACM[which(DEG_NACM$padj<0.05 & abs(DEG_NACM$log2FoldChange) > 0.6),]
choose_gene_neuron=rownames(nrDEG_NACM)
choose_matrix=group_B1A1[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_NACM,choose_matrix)
up = nrDEG_NACM[which(nrDEG_NACM$log2FoldChange > 0.6),]
down = nrDEG_NACM[which(nrDEG_NACM$log2FoldChange < -0.6),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/data/LSD_Op/ChIP_seq")
all_NACM = group_B1A1[rownames(DEG_NACM),]
all_matrix = data.frame(DEG_NACM,all_NACM)
write.csv(all_matrix,'All_B1A1_NACMF.csv')
write.csv(new_DEG_matrix,'Differential_B1A1_NACMF.csv')

dat <- read.csv("Differential_B1A1_NACMF.csv")
dat1 <- dat[,2:4]
write.table(dat1,'Differential_B1A1_NACMF_out.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)