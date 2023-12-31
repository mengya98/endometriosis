setwd('G:/endometriosis/Analysis/EPI')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(ggsci)

#S3D#
GSE203191<-readRDS('GSE203191.rds')
GSE203191.data<-data.frame(GSE203191@assays$RNA@counts,check.names = F)
GSE203191<-CreateSeuratObject(counts = GSE203191.data,Project="GSE203191")
GSE203191<-NormalizeData(object = GSE203191, normalization.method = 'LogNormalize',scale.factor = 1e6)
GSE203191<-FindVariableFeatures(GSE203191,selection.method = 'vst',nfeatures = 3000)
GSE203191.genes<-rownames(GSE203191)
GSE203191<-ScaleData(GSE203191,features = GSE203191.genes)
GSE203191<-RunPCA(GSE203191,features = VariableFeatures(object = GSE203191))
GSE203191<-FindNeighbors(GSE203191,dims = 1:9)
GSE203191<-FindClusters(GSE203191, resolution = 0.4)
GSE203191<-RunUMAP(GSE203191,dims = 1:9)

GSE203191$celltype<-'Non-cilia'
GSE203191$celltype[which(GSE203191$seurat_clusters%in%c('7'))]='Cilia'
meta<-data.frame(GSE203191@meta.data)
#cilia proportion#
meta$celltype<-factor(meta$celltype,levels = c('Non-cilia','Cilia'))
pdf("S3D, GSE203191.cilia.proportion.pdf",width = 3.5, height = 3)
ggplot(meta,aes(x=pheno,fill=celltype))+
  geom_bar(stat = 'count',
           position = 'fill',
           width = 0.5)+
  labs(x = "",y = "Proportion", title = "Data from GSE203191",fill='')+
  theme_classic()+
  theme(axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        legend.text = element_text(color="black", size = 10),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  theme(axis.text = element_text(colour = 'black'))
dev.off()

#S3E#
epi<-readRDS('./epi.rds')
epi_rmCilia<-subset(epi,celltype!=c('Ciliated'))
pdf('S3E, EPI_rmCilia_SULT1E1.pdf',width=4.8,height=4)
VlnPlot(epi_rmCilia,'SULT1E1',group.by = 'celltype',cols = c('#CD9B1D',pal_simpsons()(5)[c(2,5)]),pt.size = 0)
dev.off()
