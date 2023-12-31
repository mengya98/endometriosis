setwd('G:/endometriosis/Analysis/All/')

library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggsci)

a<-read.table('G:/endometriosis/matrix/all_p1_p22_symbol.txt',header=T,row.names=1,sep='\t')
b<-read.table('G:/endometriosis/matrix/all_p23_28_matrix_symb.txt',header=T,row.names=1,sep='\t')
c<-read.table('G:/endometriosis/matrix/all_p29_35_matrix_symb.txt',header=T,row.names=1,sep='\t')
d<-read.table('G:/endometriosis/matrix/all_p36_symb.txt',header=T,row.names=1,sep='\t')
e<-read.table('G:/endometriosis/matrix/all_p37_40_matrix_symb.txt',header=T,row.names=1,sep='\t')

endome=cbind(a,b,c,d,e)
all<-CreateSeuratObject(counts = endome, min.cells = 3, min.features =  500,project = 'Endometriosis')
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^MT")
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
patient<-sapply(strsplit(as.character(rownames(all@meta.data)),'_'),'[[',1)
all[['patient']]<-sapply(strsplit(as.character(paste(patient,'0',sep = '_')),'_'),'[[',1)
control<-sapply(strsplit(as.character(rownames(all@meta.data)),'_'),'[[',2)
all[['control']]<-sapply(strsplit(as.character(paste(control,'0',sep = '_')),'_'),'[[',1)
loci<-sapply(strsplit(as.character(rownames(all@meta.data)),'_'),'[[',3)
all[['loci']]<-sapply(strsplit(as.character(paste(loci,'0',sep = '_')),'_'),'[[',1)

all<-subset(all, subset = nCount_RNA > 10000 & nCount_RNA < 1500000 & percent.mt<25)
all<-subset(all,loci!='FLU' & loci!= 'BLOOD' & loci!= 'OV' & patient!='AD1' & patient!='P37' & patient!='P39')

all<-NormalizeData(object = all, normalization.method = 'LogNormalize',scale.factor = 1e6)
all<-FindVariableFeatures(all,selection.method = 'vst',nfeatures = 3500)

all.genes<-rownames(all)
all<-ScaleData(all,features = all.genes)
all<-RunPCA(all,features = VariableFeatures(object = all))

all<-JackStraw(all, num.replicate = 100,dims = 30)
all<-ScoreJackStraw(all,dims = 1:30)
pdf("all.pc.pdf",width = 10,height = 6)
JackStrawPlot(all,dims = 1:30)
ElbowPlot(all,ndims = 30)
dev.off()

all<-FindNeighbors(all,dims = 1:12)
all<-FindClusters(all, resolution = 0.4)
all<-RunUMAP(all,dims = 1:12)

#celltype#
all$cell_type<-''
all$cell_type[all$seurat_clusters%in%c(4,14,15,16,18,21,22)]='Epithelium'
all$cell_type[all$seurat_clusters%in%c(1,13,20)]='Macrophages'
all$cell_type[all$seurat_clusters%in%c(2,5,6,7,8,9,10,17,19)]='Stromal fibroblasts'
all$cell_type[all$seurat_clusters%in%c(0,12)]='T cells'
all$cell_type[all$seurat_clusters%in%c(3)]='NK cells'
all$cell_type[all$seurat_clusters%in%c(11)]='Endothelium'

#F1C#
all$cell_type<-factor(all$cell_type,levels=c('Epithelium','Stromal fibroblasts','T cells','Macrophages','NK cells','Endothelium'))
pdf('F1C, all_celltype.pdf',width=7.7,height=6)
DimPlot(all,reduction = "umap", pt.size = 0.5, group.by = 'cell_type',label = F,
        cols = c(pal_simpsons()(15)[c(14,2)],'#7876B1FF',pal_simpsons()(15)[c(10,3,1)]))+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()


T_markers<-c('CD3D','CD3G','CD3E','CD2')
NK_markers<-c('NCAM1','KLRD1','KLRC1','NKG7')
Macro_markers<-c('CD14','CD163','FCGR3A','CSF1R')
EPI_markers<-c('EPCAM','MUC1','CLDN4','KRT18')
STR_markers<-c('THY1','COL5A1','COL1A1','PDGFRB')
Endo_markers<-c('PECAM1','VWF','CLDN5','CDH5')

#F1D#
all$cell_type<-factor(all$cell_type,levels=rev(c('Epithelium','Stromal fibroblasts','T cells','Macrophages','NK cells','Endothelium')))
pdf('F1D, all_celltype_dotplot.pdf',width=10,height=4)
DotPlot(all,features = c(EPI_markers,STR_markers,T_markers,Macro_markers,NK_markers,Endo_markers),
        group.by = 'cell_type',dot.scale = 6)+
  scale_colour_distiller(palette = "Reds",direction = 1)+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+
  theme(axis.text= element_text(color = 'black'), axis.ticks=element_line(color = 'black'))
dev.off()

#F1E tissue types#
all$origin<-''
all$origin[all$control%in%c('N')&all$loci%in%c('EU')]<-'Normal'
all$origin[all$control%in%c('P')&all$loci%in%c('EU')]<-'Eutopic'
all$origin[all$control%in%c('P')&all$loci%in%c('EC')]<-'Ectopic'
pdf('F1E, all_tissue types.pdf',width=7,height=6)
DimPlot(all, reduction = "umap", pt.size = 0.5, group.by = 'origin',label = F,
        cols = pal_simpsons()(14)[c(4,13,7)])+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

#F1F Menstrual cycle#
all$bigphase=''
all$bigphase[which(all$patient%in%c('P1','P2','P3','P7','P11','P15','P19','P20','P22','P21','P25','P30','P34','P40','P17','P23','P28','P9','P10'))]='Proliferative'
all$bigphase[which(all$patient%in%c('P14','P16','P27','P36','P38','P12','P18','P24','P26','P29','P31','P32','P35','P8','P33'))]='Secretory'
pdf('F1F, all_menstrual_cycle.pdf',width=7.2,height=6)
DimPlot(all, reduction = "umap", pt.size = 0.5, group.by = 'bigphase',label = F,cols = c("#B09C85FF","#E64B35FF"))+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

saveRDS(all,'./all_seurat.rds')

epi_includeCilia<-subset(all,cell_type=='Epithelium')
saveRDS(epi_includeCilia,'G:/endometriosis/Analysis/EPI/epi.rds')
remove(epi_includeCilia)

tcell<-subset(all,cell_type=='T cells')
saveRDS(tcell,'G:/endometriosis/Analysis/T/tcell.rds')
remove(tcell)

str<-subset(all,cell_type=='Stromal fibroblasts')
saveRDS(str,'G:/endometriosis/Analysis/STR/str.rds')
remove(str)
