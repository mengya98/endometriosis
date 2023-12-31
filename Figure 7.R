setwd('G:/endometriosis/Analysis/STR')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(ggsci)
library(readxl)

str<-readRDS('str.rds')

meta_str_old<-data.frame(str@meta.data)

str.data<-data.frame(str@assays$RNA@counts,check.names = F)
str<-CreateSeuratObject(counts = str.data,Project="str",min.cells = 3,min.features = 500)
meta_str_new<-data.frame(str@meta.data)
meta_str_new<-cbind(meta_str_new,meta_str_old[,c(4:7,11,12)])
str@meta.data<-meta_str_new
str<-NormalizeData(object = str, normalization.method = 'LogNormalize',scale.factor = 1e6)
str<-FindVariableFeatures(str,selection.method = 'vst',nfeatures = 3000)
str.genes<-rownames(str)
str<-ScaleData(str,features = str.genes)
str<-RunPCA(str,features = VariableFeatures(object = str))
str<-FindNeighbors(str,dims = 1:10)
str<-FindClusters(str, resolution = 0.8)
str<-RunUMAP(str,dims = 1:10)

#F7A#
str$celltype<-''
str$celltype[str$seurat_clusters%in%c(0,1,8,9,10)]<-'Endometrial stroma proliferative'
str$celltype[str$seurat_clusters%in%c(5,7,12,13,14,16)]<-'Endometrial stroma secretory'
str$celltype[str$origin=='Ovary']<-'Ovarian stroma'

pdf('F7A, STR_celltype.pdf',width=6.2,height=4)
str$celltype<-factor(str$celltype,levels = (c('Endometrial stroma proliferative','Endometrial stroma secretory','Ovarian stroma')))
DimPlot(str, reduction = "umap", pt.size = 1,label = F,group.by = 'celltype',cols=c('#CD9B1D',pal_simpsons()(16)[c(2,5)]))+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme_bw()+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

#F7B#
pdf('F7B, STR_celltype_dotplot.pdf',width=9,height=3)
str$celltype<-factor(str$celltype,levels = rev(c('Endometrial stroma proliferative','Endometrial stroma secretory','Ovarian stroma')))
DotPlot(str,features = c('MMP11','PMAIP1','SFRP1','WNT5A','IGF1',
                         'S100A4','DKK1','CRYAB','FOXO1','IL15',
                         'ARX','STAR','MYH11','DCN','LMOD1'),
        group.by = 'celltype',dot.scale = 6)+
  scale_colour_distiller(palette = "Reds",direction = 1)+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
dev.off()

#F7C#
pdf('F7C, str_tissuetypes.pdf',width=4.8,height=4)
str$origin<-factor(str$origin,levels = c('Normal','Eutopic','Ectopic'))
DimPlot(str, reduction = "umap", pt.size = 1,label = F,group.by = 'origin',cols=c(pal_simpsons()(16)[c(7,13,4,5)]))+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme_bw()+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

#F7D#
goterm<-read_xlsx('STR_early_pro_EU_PvsN_up_GOterms.xlsx',sheet = 7)
goterm$Description<-factor(goterm$Description,levels = rev(goterm$Description))
pdf("F7D, STR_pro_GOterms.pdf",width = 6.6,height = 4)
ggplot(goterm,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#8B4500")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-log(Pvalue)", title = "Go terms of eutopic stroma")
dev.off()

#F7E earlypro DEGs volcano#
str<-subset(str,celltype!='Ovarian stroma')
meta<-data.frame(str@meta.data)
early_pro_EU_N<-rownames(meta[meta$phase=='early_proliferative'& meta$origin%in%c('Normal'),])
early_pro_EU_P<-rownames(meta[meta$phase=='early_proliferative'& meta$origin%in%c('Eutopic'),])

volcano_earlypro<-FindMarkers(str,ident.1 = early_pro_EU_P, ident.2 = early_pro_EU_N,only.pos = F,min.pct = 0, logfc.threshold = 0)
volcano_earlypro<-volcano_earlypro[volcano_earlypro$p_val_adj<1,]
volcano_earlypro$gene<-rownames(volcano_earlypro)
volcano_earlypro<-volcano_earlypro[order(volcano_earlypro$avg_log2FC,decreasing = T),]
volcano_earlypro$label_genes<-''
volcano_earlypro$label_genes_sig<-''
volcano_earlypro$label_genes_sig<-'con'
volcano_earlypro$label_genes_sig[volcano_earlypro$p_val_adj <= 0.05 & (volcano_earlypro$avg_log2FC >= 0.25 | volcano_earlypro$avg_log2FC <= -0.25)]<-'mark'
for (i in 1:nrow(volcano_earlypro)) {
  if (volcano_earlypro[i,6]%in%c(showgenes)) {
    volcano_earlypro[i,7]<-volcano_earlypro[i,6]
    volcano_earlypro[i,8]<-'label'
  }
}
pdf("F7E, STR_volcano_earlypro_PvsN_DEGs.pdf",width=4,height=4)
ggplot(volcano_earlypro,aes(x=avg_log2FC,y=-log10(p_val_adj),color=label_genes_sig,size=label_genes_sig))+
  geom_point()+
  theme_classic()+
  labs(x = 'log2(FoldChange)', y = '-log10(FDR)', title = "Volcano plot based on DEGs between EU and N")+
  geom_text_repel(data = volcano_earlypro, aes(x = avg_log2FC, 
                                               y = -log10(p_val_adj), 
                                               label = label_genes),
                  size = 2,box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.3, "lines"), 
                  max.overlaps = 600,
                  segment.color = "black", 
                  segment.size=0.2,
                  color="black",
                  show.legend = T)+
  scale_color_manual(name='',values =c("con"="grey",'mark'='gray24',"label"="red"))+
  scale_size_manual(values = c(0.00001,1,0.00001))+
  scale_y_continuous(limits = c(0,22))+
  geom_vline(xintercept = c(0.25,-0.25),size=0.1,color='grey')+
  geom_hline(yintercept = c(1.301),size=0.1,color='grey')+
  theme(legend.position = 'none')
dev.off() 
