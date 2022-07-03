setwd('G:/endometriosis/Analysis/p40/All/remove_FLU_BLOOD_OV/')

library(Seurat)
library(ggplot2)
library(readxl)
library(RColorBrewer)
library(pheatmap)
library(ggsci)

all<-readRDS('all_seurat.rds')

#S1A QC#
all$QC<-''
v1<-VlnPlot(object = all, features = c("nFeature_RNA"),group.by = 'QC',pt.size = 0,cols = '#91331FFF')+
  theme(axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  theme(legend.position="none")
v2<-VlnPlot(object = all, features = c("nCount_RNA"),group.by = 'QC',pt.size = 0,cols = '#91331FFF')+
  theme(axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  theme(legend.position="none")
pdf('S1A, all_QCvlnplot.pdf',width=5,height=4)
patchwork::wrap_plots(v1,v2,ncol = 2)
dev.off()

#S1B# cell number of samples
meta<-data.frame(all@meta.data)
pdf('S1B, all_patient_cellcount.pdf',width=5,height=4)
meta$patient<-factor(meta$patient,levels = c('P1','P2','P3','P7','P8','P9','P10','P11','P12','P14','P15','P16','P17','P18','P19','P20','P22','P21','P23','P24','P25','P26','P27','P28','P29','P30','P31','P32','P33','P34','P35','P36','P38','P40'))
ggplot(meta, aes(x = patient), stat = "count")+
  geom_bar(fill='black',width = 0.8)+
  theme_classic()+
  labs(x='Patient',y='Cell number')+
  theme(axis.text = element_text(color='black'),axis.ticks = element_line(color='black'))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
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

#S1C#
pdf('S1C, all_heatmap.pdf',width=6,height=4)
all$cell_type<-factor(all$cell_type,levels=(c('Epithelium','Stromal fibroblasts','T cells','Macrophages','NK cells','Endothelium')))
DoHeatmap(all,features=top50.markers$gene,
          group.by = 'cell_type',group.colors = c(pal_simpsons()(15)[c(14,2)],'#7876B1FF',pal_simpsons()(15)[c(10,3,1)]),
          label = F,draw.lines = F)+
  scale_fill_distiller(palette = "RdYlBu",direction = -1)
dev.off()

##GO terms##
#S1D#
epi_GO<-read_xlsx('./GEO/EPI_GO.xlsx',sheet = 2)
epi_GO<-epi_GO[c(1,3,22,58,162,242,202,141),]
epi_GO$Description<-factor(epi_GO$Description,levels = rev(epi_GO$Description))
pdf("S1D, EPI_top500DEGs_GOterms.pdf",width = 6.9,height = 4)
ggplot(epi_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#91331FFF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of Epithelium")
dev.off()
#S1E#
str_GO<-read_xlsx('./GEO/STR_GO.xlsx',sheet = 2)
str_GO<-str_GO[c(6,13,29,34,61,159,196,218),]
str_GO$Description<-factor(str_GO$Description,levels = rev(str_GO$Description))
pdf("S1E, STR_top500DEGs_GOterms.pdf",width = 6.4,height = 4)
ggplot(str_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#709AE1FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =1))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of Stromal fibroblasts")
dev.off()
#S1F#
T_GO<-read_xlsx('./GEO/T_GO.xlsx',sheet = 2)
T_GO<-T_GO[c(1,14,49,84,165,178,450,464),]
T_GO$Description<-factor(T_GO$Description,levels = rev(T_GO$Description))
pdf("S1F, T_top500DEGs_GOterms.pdf",width = 7.0,height = 4)
ggplot(T_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#7876B1FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of T cells")
dev.off()
#S1G#
mac_GO<-read_xlsx('./GEO/Macro_GO.xlsx',sheet = 2)
mac_GO<-mac_GO[c(3,57,65,95,111,185,218,272),]
mac_GO$Description<-factor(mac_GO$Description,levels = rev(mac_GO$Description))
pdf("S1G, Mac_top500DEGs_GOterms.pdf",width = 7,height = 4)
ggplot(mac_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#71D0F5FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of Macrophages")
dev.off()
#S1H#
nk_GO<-read_xlsx('./GEO/NK_GO.xlsx',sheet = 2)
nk_GO<-nk_GO[c(1,17,21,66,76,87,169,171),]
nk_GO$Description<-factor(nk_GO$Description,levels = rev(nk_GO$Description))
pdf("S1H, NK_top500DEGs_GOterms.pdf",width = 7.1,height = 4)
ggplot(nk_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#8A9197FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of NK cells")
dev.off()
#S1I#
endo_GO<-read_xlsx('./GEO/Endo_GO.xlsx',sheet = 2)
endo_GO<-endo_GO[c(1,7,24,50,55,114,139,169),]
endo_GO$Description<-factor(endo_GO$Description,levels = rev(endo_GO$Description))
pdf("S1I, Endo_top500DEGs_GOterms.pdf",width = 6.4,height = 4)
ggplot(endo_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#FED439FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of Endothelium")
dev.off()
