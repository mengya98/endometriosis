setwd('G:/endometriosis/Analysis/EPI')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(ggsci)

epi<-readRDS('./epi.rds')
meta<-data.frame(epi@meta.data)

epi_sec_cilia<-subset(epi,bigphase=='Secretory' & celltype%in%c('Ciliated') & origin%in%c('Normal','Eutopic'))
#F3D#
pdf('F3D, EPI_cilia_sec_SULT1E1.pdf',width=4.8,height=4)
VlnPlot(epi_sec_cilia,'SULT1E1',group.by = 'origin',cols = pal_simpsons()(14)[c(7,13)],pt.size = 0)
dev.off()

#F3G#
pdf('F3G, EPI_ciliagenes.pdf',width=1.5,height=3)
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = 'origin',... )  + 
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
StackedVlnPlot(epi_sec_cilia,c('MEF2C','RAF1','BCL2'),
               cols = pal_simpsons()(14)[c(7,13)],pt.size = 0)
dev.off()
