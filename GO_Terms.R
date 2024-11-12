###Prep
library(ggplot2)
library(ChIPseeker)
library(clusterProfiler)

###Gene Ontology
dat5 <- read.table("C1D1_PFC_combined_Annotated.bed")
dat6 <- dat5[,4]
Go_result_up <- enrichGO(dat6, 'org.Mm.eg.db', ont = "BP", keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH")
ggo <- Go_result_up
ggo@result$logP = -log10(ggo@result$p.adjust)
print(max(ggo@result$logP))
ggo@result$Description <- factor(ggo@result$Description, levels = rev(ggo@result$Description), ordered=TRUE)

sig_ggo <- ggo@result[ggo@result$p.adjust < 0.1,]
fig<-ggplot(sig_ggo[1:min(nrow(sig_ggo),15),], aes(x=logP,y=Description)) + 
  xlab("-logP") +
  ylab("Biological Processes") +
  coord_cartesian(xlim=c(0,10)) +
  theme_classic()
print(fig)
dotplot(ggo)