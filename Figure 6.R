setwd('G:/endometriosis/Analysis/T')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(ggsci)

tcell<-readRDS('../tcell.rds')
meta_old<-data.frame(tcell@meta.data)
meta_old<-meta_old[,c(1:7,11,12)]
tcell.data<-data.frame(tcell@assays$RNA@counts,check.names = F)
tcell<- CreateSeuratObject(counts = tcell.data, project = "tcell",min.cells = 3, min.features =  500)
tcell@meta.data<-meta_old
tcell <- NormalizeData(object = tcell, normalization.method = "LogNormalize", scale.factor = 1e6)
tcell <- FindVariableFeatures(tcell,selection.method = 'vst',nfeatures = 2000)
tcell.genes <- rownames(x = tcell)
tcell <- ScaleData(object = tcell, features = tcell.genes,vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"))
tcell <- RunPCA(object = tcell, features = VariableFeatures(object = tcell))
tcell<-JackStraw(tcell, num.replicate = 100,dims = 20)
tcell<-ScoreJackStraw(tcell,dims = 1:20)
JackStrawPlot(tcell,dims = 1:20)
ElbowPlot(tcell,ndims = 20)

#harmony
tcell <- tcell %>%
  RunHarmony("patient", plot_convergence = TRUE)

tcell <- FindNeighbors(object = tcell, reduction="harmony",dims = 1:9)
tcell <- FindClusters(object = tcell, resolution = 2.0)
tcell <- RunUMAP(object = tcell,reduction="harmony",dims = 1:9)

tcell$newname<-''
tcell$newname[tcell$seurat_clusters%in%c(0,6,11,12,15)]<-'CD4_TRM' 
tcell$newname[tcell$seurat_clusters%in%c(1,2,7,10,14)]<-'CD8_TRM'
tcell$newname[tcell$seurat_clusters%in%c(3,5)]<-'MAIT'
tcell$newname[tcell$seurat_clusters%in%c(4,13)]<-'TCM'
tcell$newname[tcell$seurat_clusters%in%c(8)]<-'TEM'
tcell$newname[tcell$seurat_clusters%in%c(9)]<-'CTL'
tcell$newname[tcell$seurat_clusters%in%c(16)]<-'Treg'

#F6A T celltype#
pdf('F6A, Tcell_celltype.pdf',width=5.1,height=4)
tcell$newname<-factor(tcell$newname,levels=(c('CTL','TEM','TCM',
                                              'CD4_TRM','CD8_TRM','Treg','MAIT')))
DimPlot(tcell, reduction = "umap", pt.size = 1,label = F,group.by = 'newname',
        cols = c("#A6CEE3","#1F78B4","#FB9A99","#B2DF8A","#FDBF6F","#E31A1C","#33A02C"))+
  labs(x = "UMAP1", y = "UMAP2",title = 'Cell Type')+
  theme_bw()+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

#F6B#
all_markers<-c('CD3E','CD4','CD8A','CD8B','FCGR3A','KLRD1','GZMH',
               'STMN1','HMGB2','MKI67','IL7R','TCF7','CCR7','SELL',
               'CXCR6','ZNF683','ITGAE','XCL1','FOXP3','CTLA4','IL2RA','SLC4A10','ZBTB16','DPP4','NCR3')
pdf('F6B, Tcell_celltype_dotplot.pdf',width=9,height=3)
tcell$newname<-factor(tcell$newname,levels=rev(c('CTL','TEM','TCM',
                                                 'CD4_TRM','CD8_TRM','Treg','MAIT')))
DotPlot(tcell,dot.scale = 6,features = all_markers,group.by = 'newname')+
  theme(strip.background = element_blank())+
  scale_colour_gradientn(colours = c("#9E0142","#D53E4F","#F46D43",'orange',"#FDAE61","#FEE08B","#E6F598","#BABABA","#878787","#4D4D4D","#1A1A1A"),
                         values = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
  theme(axis.text = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))
dev.off()

#F6C IFNr score#
tcell_N_EC<-subset(tcell,origin%in%c('Normal','Ectopic'))
pdf('F6C, tcell_ifnr_score(N+EC).pdf',width=3.6,height=4)
VlnPlot(tcell_N_EC,'ifnr_genes1',group.by = 'origin',cols=pal_simpsons()(14)[c(7,4)],pt.size = 0)+
  labs(title = 'IFNr signature score',x='')+
  theme_bw() + theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5),axis.text = element_text(color='black'))+
  scale_y_continuous(limits = c(-2.5,3.1),breaks = c(-2,0,2))+
  theme(legend.position = 'none')+
  geom_signif(comparisons = list(c("Normal",'Ectopic')),
              map_signif_level=T, 
              tip_length=c(0.02,0.02), 
              y_position = c(2.7),
              size=0.3,
              textsize = 3, 
              test = "t.test")
dev.off()

#F6D Tcell EPI interaction#
selected_rows = read.table('rows_use.txt',sep = '\t')
selected_rows=selected_rows$V1
selected_columns = read.table('cols_use.txt',sep = '\t')
selected_columns=selected_columns$V1
all_pval = read.table('./pvalues.txt', header=T, stringsAsFactors = F, sep='\t', comment.char = '', check.names=F)
all_means = read.table('./means.txt', header=T, stringsAsFactors = F, sep='\t', comment.char = '', check.names=F)

intr_pairs = all_pval$interacting_pair
all_pval = all_pval[,-c(1:11)]
all_means = all_means[,-c(1:11)]
if(is.null(selected_rows)){
  selected_rows = intr_pairs
}
if(is.null(selected_columns)){
  selected_columns = colnames(all_pval)
}
sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
df_names = expand.grid(selected_rows, selected_columns)
pval = unlist(sel_pval)
pval[pval==0] = 0.0009
plot.data = cbind(df_names,pval)
pr = unlist(as.data.frame(sel_means))
pr[pr==0] = 1
plot.data = cbind(plot.data,log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

filter<-plot.data[(plot.data$pair%in%c(selected_rows) 
                   & plot.data$clusters%in%c('Epi_Normal|T_Normal','Epi_Eutopic|T_Eutopic','Epi_Ectopic|T_Ectopic')), ]
my_palette <- colorRampPalette(rev(c("#9E0142","#D53E4F","#F46D43",'orange',"#FDAE61","#FEE08B","#E6F598","#BABABA","#878787","#4D4D4D","#1A1A1A")), alpha=TRUE)(n=399)

pdf('F6D, EPI_Tcell_interaction.pdf',height = 6,width = 8)  
filter$clusters<-factor(filter$clusters,levels = (c('Epi_Normal|T_Normal','Epi_Eutopic|T_Eutopic','Epi_Ectopic|T_Ectopic')))
ggplot(filter,aes(x=pair,y=clusters)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette, limit=c(-1.75,1.75),breaks=c(-1,0,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=16, colour = "black"),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.8, linetype = "solid", colour = "black"))
dev.off()
