setwd('G:/endometriosis/Analysis/EPI')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(ggsci)

#EV3A#
P36_GO<-read_xlsx('./DEG/EC/GEO/P36_ec_up_meta.xlsx',sheet = 3)
P36_GO$Description<-factor(P36_GO$Description,levels = rev(P36_GO$Description))
pdf("EV3A, P36_ECvsEU_upregulation_GOterms.pdf",width = 7.5,height = 4)
ggplot(P36_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#D2AF81FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-LogP", title = "Go terms of EC vs EU in P36")
dev.off()

#EV3B epi immune dotplot#
epi<-readRDS('./epi.rds')
epi_immu<-read.table('epi_imm_genes.txt',stringsAsFactors = F)
feature<-list('IFNr stimulated genes'=c('STAT1','STAT2','IRF1','IFIT3','IFITM1','IFITM2','IFITM3','GBP1','GBP2','GBP4'),
              'MHCII genes'=c(epi_immu$V1[c(1:4)],'HLA-DQA1',epi_immu$V1[c(6:7,5)],'HLA-DRB6'), 'Complement genes'=epi_immu$V1[8:13],
              'Chemokines'=epi_immu$V1[c(14:15,19,21,17:18)],
              'Chronic inflammatory genes'=c('SAA1','SAA2','LY6E','PTGS2','RARRES2','S100A8','S100A9'))
pdf('EV3B, EPI_immu_dotplot.pdf',width=14,height=2.5)
DotPlot(epi,group.by = 'origin',dot.scale = 6,features = feature)+
  theme_bw()+
  theme(strip.background = element_blank())+
  scale_colour_distiller(palette='Reds',direction = 1)+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
  theme(axis.text = element_text(color = 'black'))
dev.off()
