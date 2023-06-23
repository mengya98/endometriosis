setwd('G:/endometriosis/Analysis/EPI')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(ggsci)
library(monocle)

epi<-readRDS('./epi.rds')
epi.data<-data.frame(epi@assays$RNA@counts,check.names = F)
epi<-CreateSeuratObject(counts = epi.data,Project="epi",min.cells = 3,min.features = 500)
epi<-NormalizeData(object = epi, normalization.method = 'LogNormalize',scale.factor = 1e6)
epi<-FindVariableFeatures(epi,selection.method = 'vst',nfeatures = 1500)

epi.genes<-rownames(epi)
epi<-ScaleData(epi,features = epi.genes)
epi<-RunPCA(epi,features = VariableFeatures(object = epi))

epi<-JackStraw(epi, num.replicate = 100,dims = 20)
epi<-ScoreJackStraw(epi,dims = 1:20)
epi<-FindNeighbors(epi,dims = 1:8)
epi<-FindClusters(epi, resolution = 0.4)
epi<-RunUMAP(epi,dims = 1:8)

epi$celltype<-''
epi$celltype[epi$phase%in%c('early_proliferative','late_proliferative')]<-'Glandular proliferative'
epi$celltype[epi$phase%in%c('early_secretory','mid_secretory')]<-'Glandular secretory'
epi$celltype[epi$origin%in%c('Ectopic')]<-'NNMT+'
epi$celltype[epi$seurat_clusters%in%c(3)]<-'Ciliated'

#F2A epi QC#
col<-rev(c("#97d35b",#P2
            "#eddb1a",#P17
            "#6b5ee0",#P19
            "#ed7647",#P28
            "#80e5b6",#P14
            "#a2bcf2",#P16
            "#d8589f",#P27
            "#91c3d8",#P29
            "#7ec635",#P38
            "#f9c7b1",#P3
            "#5bc8cc",#P11
            "#210b66",#P15
            "#ba2337",#P21
            "#a1a5e0",#P23
            "#c16d30",#P34
            "#8473f4",#P40
            "#c4751b",#P12
            "#2eccb1",#P18
            "#40a9d6",#P31
            "#74e072",#P35
            "#83d8db"))#P36

epi<-subset(epi,patient%in%c('P2','P17','P19','P28','P14','P16','P27','P29','P38','P3','P11','P15','P21','P23','P34','P40','P12','P18','P31','P35','P36'))

epi$patient<-factor(epi$patient,levels = rev(c('P2','P17','P19','P28','P14','P16','P27','P29','P38','P3','P11','P15','P21','P23','P34','P40','P12','P18','P31','P35','P36')))
meta2<-data.frame(epi@meta.data)
meta2$patient<-factor(meta2$patient,levels = rev(c('P2','P17','P19','P28','P14','P16','P27','P29','P38','P3','P11','P15','P21','P23','P34','P40','P12','P18','P31','P35','P36')))
t1<-VlnPlot(object = epi, features = c("nFeature_RNA"),group.by = 'patient',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = col)+
  scale_y_log10()+
  coord_flip()
t2<-VlnPlot(object = epi, features = c("nCount_RNA"),group.by = 'patient',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = col)+
  scale_y_log10()+
  coord_flip()
epi_cellnumber<-data.frame(number=summary(meta2$patient))
epi_cellnumber$patient<-rownames(epi_cellnumber)
epi_cellnumber$patient<-factor(epi_cellnumber$patient,levels = rev(c('P2','P17','P19','P28','P14','P16','P27','P29','P38','P3','P11','P15','P21','P23','P34','P40','P12','P18','P31','P35','P36')))
t3<-ggplot(epi_cellnumber, aes(x = patient,y=number,fill=patient))+
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
origin_cellnumber_new<-data.frame(number=summary(meta2$origin))
origin_cellnumber_new$origin<-rownames(origin_cellnumber_new)
origin_cellnumber_new$origin<-factor(origin_cellnumber_new$origin,levels = rev(c('Normal','Eutopic','Ectopic')))
t5<-ggplot(origin_cellnumber_new, aes(x = origin,y=number,fill=origin))+
  geom_col(aes(fill=origin),width = 0.6,color='black')+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  theme_classic()+
  theme(legend.position="none")+
  labs(x='',y='Cell number')+
  theme(axis.text = element_text(color='black'),axis.ticks = element_line(color='black'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+coord_flip()#+
epi$origin<-factor(epi$origin,levels = rev(c('Normal','Eutopic','Ectopic')))
t6<-VlnPlot(object = epi, features = c("nFeature_RNA"),group.by = 'origin',pt.size = 0)+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_y_log10()+
  coord_flip()
t7<-VlnPlot(object = epi, features = c("nCount_RNA"),group.by = 'origin',pt.size = 0)+
  theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values = pal_simpsons()(14)[c(4,13,7)])+
  scale_y_log10()+
  coord_flip()

epi$patient<-factor(epi$patient,levels = (c('P2','P17','P19','P28','P14','P16','P27','P29','P38','P3','P11','P15','P21','P23','P34','P40','P12','P18','P31','P35','P36')))
meta2<-data.frame(epi@meta.data)
sankey2<-meta2[,c(5,8)]

sankey2<-group_by(sankey2, patient,origin) %>% summarise(., count = n())
sankey2$patient<-factor(sankey2$patient,levels = (c('P2','P17','P19','P28','P14','P16','P27','P29','P38','P3','P11','P15','P21','P23','P34','P40','P12','P18','P31','P35','P36')))
t4<-ggplot(data = sankey2,aes(axis1 = patient, axis2 = origin,y = count,label=F))+
  scale_x_discrete(limits = c("patient", "origin"), expand = c(.1, .05)) +#
  geom_alluvium(aes(fill = patient),width = 1/8,color='black') +
  geom_stratum(width = 1/8,size=0.1) + 
  scale_fill_manual(values = rev(col))+
  theme_minimal() +
  theme(panel.border = element_blank(),panel.grid=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank())+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  theme(legend.position="none")+
  ggtitle("Patient ID -> Locations")


pdf('F2A, epi_QCplot.pdf',width=16,height=6)
t2+t1+t3+t4+t5+t6+t7+plot_layout(widths = c(1,1,1,4,1,1,1),ncol = 7)
dev.off()

#F2B epi celltype#
epi$celltype<-factor(epi$celltype,levels = c('Ciliated','Glandular proliferative','Glandular secretory','NNMT+'))
pdf('F2B, epi_celltype.pdf',width=6.6,height=5)
DimPlot(epi, reduction = "umap", group.by = 'celltype',pt.size = 1,label = F,
        cols = c('#7876B1FF','#CD9B1D',pal_simpsons()(5)[c(2,5)]))+
  theme_bw()+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme(axis.text= element_blank(), axis.ticks=element_blank())+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))
dev.off()

#F2C epi celltype dotplot#
epi$celltype<-factor(epi$celltype,levels = rev(c('Ciliated','Glandular proliferative','Glandular secretory','NNMT')))
pdf('F2C, epi_dotplot.pdf',width=9,height=3)
DotPlot(epi,features = c('PIFO','FOXJ1','AGR3','SPAG17','SNTN',
                         'MMP7','IHH','EMID1','NPAS3','CPM',
                         'MT1E','SPP1','C2CD4A','PAEP','CXCL14',
                         'NNMT','CXCL1','HLA-DRB1','MUC16','FBLN2'),
        group.by = 'celltype',dot.scale = 6)+
  scale_colour_distiller(palette = "Reds",direction = 1)+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
dev.off()

#F2D epi tissue types#
pdf('F2D, EPI_tissue types.pdf',width=4.8,height=4)
epi$origin<-factor(epi$origin,levels = c('Normal','Eutopic','Ectopic'))
DimPlot(epi, reduction = "umap", pt.size = 1, group.by = 'origin',label = F,
        cols = pal_simpsons()(14)[c(7,13,4)])+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme_bw()+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

#F2E#
meta<-data.frame(epi@meta.data)
pro_EU_N<-rownames(meta[meta$origin=='Normal' & meta$celltype%in%c('Glandular proliferative'),])
pro_EU_P<-rownames(meta[meta$origin=='Eutopic' & meta$celltype%in%c('Glandular proliferative'),])
volcano_pro<-FindMarkers(epi,ident.1 = pro_EU_P, ident.2 = pro_EU_N,only.pos = F,min.pct = 0.25, logfc.threshold = 0)
volcano_pro<-volcano_pro[volcano_pro$p_val_adj<1,]
volcano_pro$gene<-rownames(volcano_pro)
volcano_pro<-volcano_pro[order(volcano_pro$avg_log2FC,decreasing = T),]

showgenes<-read.table('./volcano_DEGs.txt')
showgenes<-showgenes$V1
volcano_pro$label_genes<-''
volcano_pro$label_genes_sig<-''
volcano_pro$label_genes_sig<-'con'
volcano_pro$label_genes_sig[volcano_pro$p_val_adj <= 0.05 & (volcano_pro$avg_log2FC >= 0.25 | volcano_pro$avg_log2FC <= -0.25)]<-'mark'
for (i in 1:nrow(volcano_pro)) {
  if (volcano_pro[i,6]%in%c(showgenes)) {
    volcano_pro[i,7]<-volcano_pro[i,6]
    volcano_pro[i,8]<-'label'
  }
}

pdf("F2E, volcano_pro_PvsN_DEGs.pdf",width=4,height=4)
ggplot(volcano_pro,aes(x=avg_log2FC,y=-log10(p_val_adj),color=label_genes_sig,size=label_genes_sig))+geom_point()+
  theme_classic()+
  labs(x = 'log2(FoldChange)', y = '-log10(FDR)', title = "Volcano plot based on DEGs between EU and N")+
  geom_text_repel(data = volcano_pro, aes(x = avg_log2FC, 
                                          y = -log10(p_val_adj), 
                                          label = label_genes),
                  size = 2,box.padding = unit(0.1, "lines"),#shape=18,
                  point.padding = unit(0.3, "lines"), 
                  max.overlaps = 100,
                  segment.color = "black", 
                  segment.size=0.2,
                  color="black",
                  show.legend = FALSE)+
  scale_color_manual(name='',values =c("con"="grey",'mark'='gray24',"label"="red"))+
  scale_size_manual(values = c(0.1,1,0.1))+
  geom_vline(xintercept = c(0.25,-0.25),size=0.1,color='grey')+
  geom_hline(yintercept = c(1.301),size=0.1,color='grey')+
  theme(legend.position = 'none')
dev.off() 

#F2F-F2I#
data <- as(as.matrix(epi@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = epi@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
var.genes <- FindVariableFeatures(epi,nfeatures = 1000)
var.genes <- VariableFeatures(var.genes)
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, var.genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 3,method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "State")
cds <- orderCells(cds, root_state = 2)

#F2F#
pdf("F2F, trainmonocle_celltype.pdf",width =4.5,height = 5)
plot_cell_trajectory(cds_1, color_by = "celltype", cell_size = 1,show_backbone=TRUE)+
  scale_color_manual(values=c('#7876B1FF','#CD9B1D',pal_simpsons()(5)[c(2)]))
dev.off()
#F2G#
pdf("F2G, trainmonocle_origin.pdf",width = 4.5,height = 5)
plot_cell_trajectory(cds_1, color_by = "origin", cell_size = 1,show_backbone=T)+
  scale_color_manual(values=c(pal_simpsons()(14)[c(7,13,4)]))
dev.off()
#F2H#
pdf("F2H, trainmonocle_Pseudotime.pdf",width = 4.3,height = 5)
plot_cell_trajectory(cds_1, color_by = "Pseudotime", cell_size = 1,show_backbone=TRUE)+
  scale_colour_gradientn(colours = c("#9E0142","#D53E4F","#F46D43",'orange',"#FDAE61","#FEE08B","#E6F598","#BABABA","#878787","#4D4D4D","#1A1A1A"),
                         values = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
dev.off()

exprData <- data.frame(t(epi@assays$RNA@data),check.rows = F,check.names = F)
cds$FOXJ1 <- exprData[,"FOXJ1"] 
cds$TPPP3 <- exprData[,"TPPP3"]
cds$PIFO <- exprData[,"PIFO"]
cds$RSPH1 <- exprData[,"RSPH1"]

p1<-plot_cell_trajectory(cds, color_by = "FOXJ1", cell_size = 1,show_backbone=T)+
  scale_colour_gradientn(colours = c("#9E0142","#9E0142","#D53E4F","#D53E4F","#F46D43",'orange',"#FEE08B","#E6F598","#E6F598","#BABABA","#878787"),
                         values = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),limit=c(-0.1,9.2))+
  theme(axis.line.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  labs(title = "FOXJ1")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")
p2<-plot_cell_trajectory(cds, color_by = "TPPP3", cell_size = 1,show_backbone=T)+
  scale_colour_gradientn(colours = c("#9E0142","#9E0142","#D53E4F","#D53E4F","#F46D43",'orange',"#FEE08B","#E6F598","#E6F598","#BABABA","#878787"),
                         values = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),limit=c(-0.1,9.2))+
  theme_classic()+
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())+
  labs(title = "TPPP3")+
  theme(plot.title = element_text(hjust = 0.5))
p3<-plot_cell_trajectory(cds, color_by = "PIFO", cell_size = 1,show_backbone=T)+
  scale_colour_gradientn(colours = c("#9E0142","#9E0142","#D53E4F","#D53E4F","#F46D43",'orange',"#FEE08B","#E6F598","#E6F598","#BABABA","#878787"),
                         values = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),limit=c(-0.1,9.2))+
  labs(title = "PIFO")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")
p4<-plot_cell_trajectory(cds, color_by = "RSPH1", cell_size = 1,show_backbone=T)+
  scale_colour_gradientn(colours = c("#9E0142","#9E0142","#D53E4F","#D53E4F","#F46D43",'orange',"#FEE08B","#E6F598","#E6F598","#BABABA","#878787"),
                         values = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),limit=c(-0.1,9.2))+
  labs(title = "RSPH1")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line.y = element_blank(),axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.title.y = element_blank())+
  theme(legend.position="none")

#F2I#
pdf("F2I, trainmonocle_ciliagenes.pdf",width = 8.5,height = 10)
patchwork::wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()

