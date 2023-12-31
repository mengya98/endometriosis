setwd('G:/endometriosis/Analysis/All/')

library(Seurat)
library(ggplot2)
library(readxl)
library(RColorBrewer)
library(pheatmap)
library(ggsci)

all<-readRDS('all_seurat.rds')

#col<-randomColor(34)
col<-rev(c("#97d35b",#P2
           "#0e51d6",#P9
           "#eddb1a",#P17
           "#6b5ee0",#P19
           "#fcc2f1",#P25
           "#ed7647",#P28
           "#80e5b6",#P14
           "#a2bcf2",#P16
           "#d8589f",#P27
           "#91c3d8",#P29
           "#7ec635",#P38
           "#89b4e5",#P1
           "#40a9d6",#P31
           "#74e072",#P35
           "#f9c7b1",#P3
           "#edb39c",#P7
           "#5bc8cc",#P11
           "#210b66",#P15
           "#db3262",#P20
           "#ba2337",#P21
           "#a1a5e0",#P23
           "#95a6f4",#P30
           "#c16d30",#P34
           "#8473f4",#P40
           "#c4751b",#P12
           "#2eccb1",#P18
           "#5bef58",#P26
           "#edb0a1",#P32
           "#83d8db",#P36
           "#c69c29",#P8
           "#f7c309",#P10
           "#d374e0",#P22
           "#167f0c",#P24
           "#91f97f"))#P33

#Figure S1A QC#
all$patient<-factor(all$patient,levels = rev(c('P2','P9','P17','P19','P25','P28','P14','P16','P27','P29','P38',
                                               'P1','P31','P35',
                                               'P3','P7','P11','P15','P20','P21','P23','P30','P34','P40','P12','P18','P26','P32','P36',
                                               'P8','P10','P22','P24','P33')))
meta2<-data.frame(all@meta.data)
meta2$patient<-factor(meta2$patient,levels = rev(c('P2','P9','P17','P19','P25','P28','P14','P16','P27','P29','P38','P1','P31','P35','P3','P7','P11','P15','P20','P21','P23','P30','P34','P40','P12','P18','P26','P32','P36','P8','P10','P22','P24','P33')))

t1<-VlnPlot(object = all, features = c("nFeature_RNA"),group.by = 'patient',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = col)+
  scale_y_log10()+
  coord_flip()
t2<-VlnPlot(object = all, features = c("nCount_RNA"),group.by = 'patient',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = col)+
  scale_y_log10()+
  coord_flip()
all_cellnumber_new<-data.frame(number=summary(meta2$patient))
all_cellnumber_new$patient<-rownames(all_cellnumber_new)
all_cellnumber_new$patient<-factor(all_cellnumber_new$patient,levels = rev(c('P2','P9','P17','P19','P25','P28','P14','P16','P27','P29','P38','P1','P31','P35','P3','P7','P11','P15','P20','P21','P23','P30','P34','P40','P12','P18','P26','P32','P36','P8','P10','P22','P24','P33')))
t3<-ggplot(all_cellnumber_new, aes(x = patient,y=number,fill=patient))+
  geom_col(aes(fill=patient),color='black')+
  scale_fill_manual(values = col)+
  theme_classic()+
  theme(legend.position="none")+
  labs(x='',y='Cell number')+
  theme(axis.text = element_text(color='black'),axis.ticks = element_line(color='black'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+coord_flip()+
  scale_x_discrete(position = "top")+
  scale_y_reverse()
meta2<-data.frame(all@meta.data)
meta2$origin<-factor(meta2$origin,levels = c('Normal','Eutopic','Ectopic'))
origin_cellnumber_new<-data.frame(number=summary(meta2$origin))
origin_cellnumber_new$origin<-rownames(origin_cellnumber_new)
origin_cellnumber_new$origin<-factor(origin_cellnumber_new$origin,levels = rev(c('Normal','Eutopic','Ectopic')))
t5<-ggplot(origin_cellnumber_new, aes(x = origin,y=number,fill=origin))+
  geom_col(aes(fill=origin),color='black',width = 0.6)+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  theme_classic()+
  theme(legend.position="none")+
  labs(x='',y='Cell number')+
  theme(axis.text = element_text(color='black'),axis.ticks = element_line(color='black'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+coord_flip()
all$origin<-factor(all$origin,levels = rev(c('Normal','Eutopic','Ectopic')))
t6<-VlnPlot(object = all, features = c("nFeature_RNA"),group.by = 'origin',pt.size = 0)+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_y_log10()+
  coord_flip()
t7<-VlnPlot(object = all, features = c("nCount_RNA"),group.by = 'origin',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  scale_y_log10()+
  coord_flip()
##sankey#
all$patient<-factor(all$patient,levels = (c('P2','P9','P17','P19','P25','P28','P14','P16','P27','P29','P38','P1','P31','P35','P3','P7','P11','P15','P20','P21','P23','P30','P34','P40','P12','P18','P26','P32','P36','P8','P10','P22','P24','P33')))
all$origin<-factor(all$origin,levels = (c('Normal','Eutopic','Ectopic')))
sankey2<-meta2[,c(5,11)]
sankey2<-group_by(sankey2, patient,origin) %>% summarise(., count = n())
sankey2$patient<-factor(sankey2$patient,levels = (c('P2','P9','P17','P19','P25','P28','P14','P16','P27','P29','P38','P1','P31','P35','P3','P7','P11','P15','P20','P21','P23','P30','P34','P40','P12','P18','P26','P32','P36','P8','P10','P22','P24','P33')))
t4<-ggplot(data = sankey2,aes(axis1 = patient, axis2 = origin,y = count,label=F))+#fill=patient,
  scale_x_discrete(limits = c("patient", "origin"), expand = c(.1, .05)) +#
  geom_alluvium(aes(fill = patient),width = 1/8,color='black') +
  geom_stratum(width = 1/8,size=0.1)+ 
  scale_fill_manual(values = rev(col))+
  theme_minimal() +
  theme(panel.border = element_blank(),panel.grid=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank())+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  theme(legend.position="none")+
  ggtitle("Participant ID->Tissue type")

pdf('Figure S1A, all_QCplot.pdf',width=16,height=5)
t2+t1+t3+t4+t5+t6+t7+plot_layout(widths = c(1,1,1,3,1,1,1),ncol = 7)
dev.off()


#find DEGs#
all$cell_type<-factor(all$cell_type,levels=c('Epithelium','Stromal fibroblasts','T cells','Macrophages','NK cells','Endothelium'))
all.markers<-FindAllMarkers(t1,logfc.threshold = 0.25,min.pct = 0.5)
all.markers_order<-all.markers[order(all.markers$cluster,all.markers$avg_log2FC,decreasing = T),]
all.markers_order<-all.markers_order[order(all.markers_order$cluster,decreasing = F),]
all.markers_order_positive<-all.markers_order[all.markers_order$avg_log2FC>=0&all.markers_order$p_val_adj<=0.05,]
all.markers_order_positive %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)-> top500.markers
write.table(top500.markers,file = "./DEG/p1_p40_celltype_top500markers.xls", row.names = F, col.names = T, quote = F, sep = "\t")
all.markers_order_positive %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)-> top50.markers
write.table(top50.markers,file = "./DEG/p1_p40_celltype_top50markers.xls", row.names = F, col.names = T, quote = F, sep = "\t")

#Figure S1B###GO terms##
epi_GO<-read_xlsx('./GEO/EPI_GO.xlsx',sheet = 2)
epi_GO<-epi_GO[c(1,3,22,58,162,242,202,141),]
epi_GO$Description<-factor(epi_GO$Description,levels = rev(epi_GO$Description))
pdf("Figure S1B, EPI_top500DEGs_GOterms.pdf",width = 6.9,height = 4)
ggplot(epi_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#91331FFF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of Epithelium")
dev.off()

str_GO<-read_xlsx('./GEO/STR_GO.xlsx',sheet = 2)
str_GO<-str_GO[c(6,13,29,34,61,159,196,218),]
str_GO$Description<-factor(str_GO$Description,levels = rev(str_GO$Description))
pdf("Figure S1B, STR_top500DEGs_GOterms.pdf",width = 6.4,height = 4)
ggplot(str_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#709AE1FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =1))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of Stromal fibroblasts")
dev.off()

T_GO<-read_xlsx('./GEO/T_GO.xlsx',sheet = 2)
T_GO<-T_GO[c(1,14,49,84,165,178,450,464),]
T_GO$Description<-factor(T_GO$Description,levels = rev(T_GO$Description))
pdf("Figure S1B, T_top500DEGs_GOterms.pdf",width = 7.0,height = 4)
ggplot(T_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#7876B1FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of T cells")
dev.off()

mac_GO<-read_xlsx('./GEO/Macro_GO.xlsx',sheet = 2)
mac_GO<-mac_GO[c(3,57,65,95,111,185,218,272),]
mac_GO$Description<-factor(mac_GO$Description,levels = rev(mac_GO$Description))
pdf("Figure S1B, Mac_top500DEGs_GOterms.pdf",width = 7,height = 4)
ggplot(mac_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#71D0F5FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of Macrophages")
dev.off()

nk_GO<-read_xlsx('./GEO/NK_GO.xlsx',sheet = 2)
nk_GO<-nk_GO[c(1,17,21,66,76,87,169,171),]
nk_GO$Description<-factor(nk_GO$Description,levels = rev(nk_GO$Description))
pdf("Figure S1B, NK_top500DEGs_GOterms.pdf",width = 7.1,height = 4)
ggplot(nk_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#8A9197FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of NK cells")
dev.off()

endo_GO<-read_xlsx('./GEO/Endo_GO.xlsx',sheet = 2)
endo_GO<-endo_GO[c(1,7,24,50,55,114,139,169),]
endo_GO$Description<-factor(endo_GO$Description,levels = rev(endo_GO$Description))
pdf("Figure S1B, Endo_top500DEGs_GOterms.pdf",width = 6.4,height = 4)
ggplot(endo_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#FED439FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of Endothelium")
dev.off()

#Figure S1C#
pdf('Figure S1C, all_patient.pdf',width=7.3,height=6)
DimPlot(all, reduction = "umap", pt.size = 0.5, group.by = 'patient',label = F)+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  scale_color_manual(values = rev(col))+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()
