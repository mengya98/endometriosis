setwd('G:/endometriosis/Analysis/p40/EPI/remove_cilia/all/')

library(Seurat)
library(ggplot2) 

#F5A epi immune dotplot#
epi_immu<-read.table('epi_imm_genes.txt',stringsAsFactors = F)
feature<-list('IFNr stimulated genes'=c('STAT1','STAT2','IRF1','IFIT3','IFITM1','IFITM2','IFITM3','GBP1','GBP2','GBP4'),
              'MHCII genes'=c(epi_immu$V1[c(1:4)],'HLA-DQA1',epi_immu$V1[c(6:7,5)],'HLA-DRB6'), 'Complement genes'=epi_immu$V1[8:13],
              'Chemokines'=epi_immu$V1[c(14:15,19,21,17:18)],
              'Chronic inflammatory genes'=c('SAA1','SAA2','LY6E','PTGS2','RARRES2','S100A8','S100A9'))
pdf('F5A, EPI_immu_dotplot.pdf',width=14,height=2.5)
DotPlot(epi,group.by = 'origin',dot.scale = 6,features = feature)+
  theme_bw()+
  theme(strip.background = element_blank())+
  scale_colour_distiller(palette='Reds',direction = 1)+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
  theme(axis.text = element_text(color = 'black'))
dev.off()