setwd('G:/endometriosis/Analysis/EPI')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(ggsci)
library(readxl)

epi<-readRDS('./epi.rds')

#F5A P15/P23 EC GO Terms#
P15_GO<-read_xlsx('./DEG/EC/GEO/p15_ec_up_meta.xlsx',sheet = 3)
P15_GO$Description<-factor(P15_GO$Description,levels = rev(P15_GO$Description))
pdf("F5A, P15_ECvsEU_up_GOterms.pdf",width = 10.8,height = 4)
ggplot(P15_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#D2AF81FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-LogP", title = "Go terms of EC vs EU in P15")
dev.off()

P23_GO<-read_xlsx('./DEG/EC/GEO/P23_ec_up_meta.xlsx',sheet = 3)
P23_GO$Description<-factor(P23_GO$Description,levels = rev(P23_GO$Description))
pdf("F5A, P23_ECvsEU_upregulation_GOterms.pdf",width = 6.7,height = 4)
ggplot(P23_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#D2AF81FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-LogP", title = "Go terms of EC vs EU in P23")
dev.off()


#F5B EPI Immune StackedVlnPlot#
epi_immu<-read.table('epi_imm_genes.txt',stringsAsFactors = F)
feature<-list('IFNr stimulated genes'=c('STAT1','STAT2','IRF1','IFIT3','IFITM1','IFITM2','IFITM3','GBP1','GBP2','GBP4'),
              'MHCII genes'=c(epi_immu$V1[c(1:4)],'HLA-DQA1',epi_immu$V1[c(6:7,5)],'HLA-DRB6'), 'Complement genes'=epi_immu$V1[8:13],
              'Chemokines'=epi_immu$V1[c(14:15,19,21,17:18)],
              'Chronic inflammatory genes'=c('SAA1','SAA2','LY6E','PTGS2','RARRES2','S100A8','S100A9'))
patient<-subset(epi,patient%in%c('P15','P23','P36'))
patient$patient_origin<-paste0(patient$patient,'_',patient$origin)
patient$patient_origin<-factor(patient$patient_origin,levels = c('P15_Eutopic','P15_Ectopic','P23_Eutopic','P23_Ectopic','P36_Eutopic','P36_Ectopic'))
pdf('F5B, stackvlnplot_EPI_immu_P15P23P36.pdf',width=3,height=16)
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = 'patient_origin',... )  + 
    xlab("") + ylab(feature) + ggtitle(feature) + theme_classic()+
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5,color ='black'),
          plot.margin = plot.margin, title = element_blank())
  return(p)
}
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,color ='black'), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
StackedVlnPlot(patient, features = c(feature$`IFNr stimulated genes`[c(2:3,5:7,9:10)],feature$`MHCII genes`[c(1:2,4:6,8:9)],feature$`Complement genes`[c(1:2,4:6)],feature$Chemokines[c(1:3,5:6)],feature$`Chronic inflammatory genes`), 
               cols = pal_simpsons()(14)[c(13,4,13,4,13,4)],pt.size = 0)
dev.off()
