setwd('G:/endometriosis/Analysis/p40/Epi_immune')
library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(patchwork)
library(ggsci)

epi<-readRDS('G:/endometriosis/Analysis/p40/EPI/remove_cilia/all/epi_removeCilia_all.rds')
mac<-readRDS('G:/endometriosis/Analysis/p40/Immune/Macro/myeloid.rds')

combined <- merge(epi, y = mac, add.cell.ids = c("EPI", "Macro"), project = "epi_macro")
celltype<-sapply(strsplit(as.character(rownames(combined@meta.data)),'_'),'[[',1)
combined[['celltype']]<-sapply(strsplit(as.character(paste(celltype,'0',sep = '_')),'_'),'[[',1)

meta<-data.frame(combined@meta.data)
##epi_macro 6+1##
combined$epi_macro<-''
combined$epi_macro<-paste0(combined$celltype,'_',combined$state_loci_phase2)
combined$epi_macro[combined$epi_macro=='Macro_NA']<-'Macro'
pal<-c('#BC3C29','#0072B5','#E18727','#20854E','#7876B1','#6F99AD','#71D0F5FF')
combined$epi_macro<-factor(combined$epi_macro,levels = c('EPI_Normal_proliferative','EPI_Eutopic_proliferative','EPI_Ectopic_proliferative',
                                                         'EPI_Normal_secretory','EPI_Eutopic_secretory','EPI_Ectopic_secretory','Macro'))

epi_immu<-read.table('G:/endometriosis/gene_list/epi_imm_genes.txt',stringsAsFactors = F)
feature<-list('IFNr stimulated genes'=c('STAT1','STAT2','IRF1','IFIT3','IFITM1','IFITM2','IFITM3','GBP1','GBP2','GBP4'),
              'MHCII genes'=c(epi_immu$V1[c(1:4)],'HLA-DQA1',epi_immu$V1[c(6:7,5)],'HLA-DRB6'), 'Complement genes'=epi_immu$V1[8:13],
              'Chemokines'=epi_immu$V1[c(14:15,19,21,17:18)],
              'Chronic inflammatory genes'=c('SAA1','SAA2','LY6E','PTGS2','RARRES2','S100A8','S100A9'))
#S4A#
pdf('S4A, stackvlnplot_epi_macro_immu_genes.pdf',width=16,height=6)
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = 'epi_macro',... )  + 
    xlab("") + ylab('') + ggtitle(feature) + theme_classic()+
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = rel(1),color ='black'), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_line(color ='black'),
          axis.title.x = element_blank(),
          plot.margin = plot.margin,plot.title = element_text(hjust = 0.5,color ='black'))
  return(p)
}
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,color ='black'), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 7)
  return(p)
}
StackedVlnPlot(combined, features = c(feature$`MHCII genes`,feature$`Complement genes`,feature$Chemokines,feature$`Chronic inflammatory genes`), 
               pt.size=0, cols=pal)
dev.off()
