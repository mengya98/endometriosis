setwd('G:/endometriosis/Analysis/T')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(ggsci)
library(ggrepel)
library(ggsignif)
library(ggalluvial)
library(patchwork)

tcell<-readRDS('./tcell.rds')

#EV4A QC#
col<-rev(c("#fcc2f1",#P25
  "#ed7647",#P28
  "#a2bcf2",#P16
  "#d8589f",#P27
  "#91c3d8",#P29
  "#f9c7b1",#P3
  "#edb39c",#P7
  "#c69c29",#P8
  "#5bc8cc",#P11
  "#210b66",#P15
  "#db3262",#P20
  "#d374e0",#P22
  "#a1a5e0",#P23
  "#95a6f4",#P30
  "#c16d30",#P34
  "#5bef58",#P26
  "#edb0a1",#P32
  "#91f97f",#P33
  "#74e072",#P35
  "#83d8db"))#P36


tcell<-subset(tcell,patient%in%c('P25','P28','P16','P27','P29','P3','P7','P8','P11','P15','P20','P22','P23','P30','P34','P26','P32','P33','P35','P36'))
tcell$patient<-factor(tcell$patient,levels = rev(c('P25','P28','P16','P27','P29','P3','P7','P8','P11','P15','P20','P22','P23','P30','P34','P26','P32','P33','P35','P36')))
meta2<-data.frame(tcell@meta.data)
meta2$patient<-factor(meta2$patient,levels = rev(c('P25','P28','P16','P27','P29','P3','P7','P8','P11','P15','P20','P22','P23','P30','P34','P26','P32','P33','P35','P36')))
v1<-VlnPlot(object = tcell, features = c("nFeature_RNA"),group.by = 'patient',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = col)+
  scale_y_log10()+
  coord_flip()
v2<-VlnPlot(object = tcell, features = c("nCount_RNA"),group.by = 'patient',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = col)+
  scale_y_log10()+
  coord_flip()
tcell_cellnumber<-data.frame(number=summary(meta2$patient))
tcell_cellnumber$patient<-rownames(tcell_cellnumber)
tcell_cellnumber$patient<-factor(tcell_cellnumber$patient,levels = rev(c('P25','P28','P16','P27','P29','P3','P7','P8','P11','P15','P20','P22','P23','P30','P34','P26','P32','P33','P35','P36')))
v3<-ggplot(tcell_cellnumber, aes(x = patient,y=number,fill=patient))+
  geom_col(aes(fill=patient),color='black')+
  scale_fill_manual(values = col)+
  theme_classic()+
  theme(legend.position="none")+
  labs(x='',y='Cell number')+
  theme(axis.text = element_text(color='black'),axis.ticks = element_line(color='black'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+coord_flip()+
  scale_x_discrete(position = "top")+
  scale_y_reverse()

meta2$origin<-factor(meta2$origin,levels = c('Normal','Eutopic','Ectopic'))
origin_cellnumber<-data.frame(number=summary(meta2$origin))
origin_cellnumber$origin<-rownames(origin_cellnumber)
origin_cellnumber$origin<-factor(origin_cellnumber$origin,levels = rev(c('Normal','Eutopic','Ectopic')))
v5<-ggplot(origin_cellnumber, aes(x = origin,y=number,fill=origin))+
  geom_col(aes(fill=origin),color='black',width = 0.6)+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  theme_classic()+
  theme(legend.position="none")+
  labs(x='',y='Cell number')+
  theme(axis.text = element_text(color='black'),axis.ticks = element_line(color='black'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+coord_flip()#+

tcell$origin<-factor(tcell$origin,levels = rev(c('Normal','Eutopic','Ectopic')))
v6<-VlnPlot(object = tcell, features = c("nFeature_RNA"),group.by = 'origin',pt.size = 0)+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_y_log10()+
  coord_flip()
v7<-VlnPlot(object = tcell, features = c("nCount_RNA"),group.by = 'origin',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  scale_y_log10()+
  coord_flip()

#sankey#
tcell$patient<-factor(tcell$patient,levels = (c('P25','P28','P16','P27','P29','P3','P7','P8','P11','P15','P20','P22','P23','P30','P34','P26','P32','P33','P35','P36')))
tcell$origin<-factor(tcell$origin,levels = (c('Normal','Eutopic','Ectopic')))
sankey<-meta2[,c(5,8)]
sankey<-group_by(sankey, patient,origin) %>% summarise(., count = n())
sankey<-sankey[c(1:19,21:24,26:30),]
sankey$patient<-factor(sankey$patient,levels = (c('P25','P28','P16','P27','P29','P3','P7','P8','P11','P15','P20','P22','P23','P30','P34','P26','P32','P33','P35','P36')))
v4<-ggplot(data = sankey,aes(axis1 = patient, axis2 = origin,y = count,label=F))+#fill=patient,
  scale_x_discrete(limits = c("patient", "origin"), expand = c(.1, .05)) +#
  geom_alluvium(aes(fill = patient),width = 1/8,color='black') +
  geom_stratum(width = 1/8,size=0.1) + 
  scale_fill_manual(values = rev(col))+
  theme_minimal() +
  theme(panel.border = element_blank(),panel.grid=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank())+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  theme(legend.position="none")+
  ggtitle("Patient ID -> Locations")

pdf('EV4A, tcell_QCplot.pdf',width=16,height=6)
v2+v1+v3+v4+v5+v6+v7+plot_layout(widths = c(1,1,1,3,1,1,1),ncol = 7)
dev.off()
