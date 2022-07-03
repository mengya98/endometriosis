setwd('G:/endometriosis/Analysis/p40/EU')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(ggsci)
library(readxl)
library(uwot)

epi<-readRDS('./epi_removeEC.rds')
meta_epi_old<-data.frame(epi@meta.data)
epi.data<-data.frame(epi@assays$RNA@counts,check.names = F)
epi<-CreateSeuratObject(counts = epi.data,Project="epi",min.cells = 3,min.features = 500)
meta_epi_new<-data.frame(epi@meta.data)
meta_epi_new<-cbind(meta_epi_new,meta_epi_old[,c(4:9,12)])
epi@meta.data<-meta_epi_new
epi
epi<-NormalizeData(object = epi, normalization.method = 'LogNormalize',scale.factor = 1e6)
epi<-FindVariableFeatures(epi,selection.method = 'vst',nfeatures = 2000)

epi.genes<-rownames(epi)
epi<-ScaleData(epi,features = epi.genes)
epi<-RunPCA(epi,features = VariableFeatures(object = epi))

epi<-JackStraw(epi, num.replicate = 100,dims = 20)
epi<-ScoreJackStraw(epi,dims = 1:20)
pdf("epi.pc.pdf",width = 10,height = 6)
JackStrawPlot(epi,dims = 1:20)
ElbowPlot(epi,ndims = 20)
dev.off()

epi<-FindNeighbors(epi,dims = 1:8)
epi<-FindClusters(epi, resolution = 0.4)
epi<-RunUMAP(epi,dims = 1:8)

##cilia cluster##
epi$origin[epi$seurat_clusters%in%c(3)]<-'Ciliated epithelium'
#remove the cells that pseudotime inconsistent with real time#
cell_remove<-c('P23_P_EU_EPI_SC55','P23_P_EU_EPI_SC34','P23_P_EU_EPI_SC70','P19_N_EU_EP_SC85',
               'P14_N_EU_M_SC78','P27_N_EU_EPI_SC26','P27_N_EU_EPI_SC87','P27_N_EU_EPI_SC92',
               'P27_N_EU_EPI_SC2','P16_N_EU_EP_SC80','P14_N_EU_M_SC81','P28_N_EU_STR_SC79',
               'P23_P_EU_EPI_SC28','P28_N_EU_STR_SC70','P28_N_EU_STR_SC80','P18_P_EU_EP_SC31',
               'P28_N_EU_STR_SC25','P28_N_EU_STR_SC38','P18_P_EU_EP_SC46','P18_P_EU_EP_SC1',
               'P18_P_EU_EP_SC40','P36_P_EU_EPI_SC92','P28_N_EU_EPI_SC31','P28_N_EU_STR_SC57',
               'P16_N_EU_IMM_SC32','P18_P_EU_EP_SC49','P19_N_EU_EP_SC96')
epi<-epi[,!((colnames(epi))%in%cell_remove)]
epi<-subset(epi,loci=='EU')

#F2A epi celltype#
epi$big_phase_cilia<-''
epi$big_phase_cilia[epi$phase%in%c('early_proliferative','late_proliferative')]<-'MMP7 epithelium'
epi$big_phase_cilia[epi$phase%in%c('early_secretory','mid_secretory')]<-'Glandular epithelium'
epi$big_phase_cilia[epi$seurat_clusters%in%c(3)]<-'Ciliated epithelium'
epi$big_phase_cilia<-factor(epi$big_phase_cilia,levels = c('MMP7 epithelium','Glandular epithelium','Ciliated epithelium'))
pdf('F2A, epi_bigphase_cilia.pdf',width=5.6,height=4)
DimPlot(epi, reduction = "umap", pt.size = 1,label = F,group.by = 'big_phase_cilia',
        cols = c('#CD9B1D',pal_simpsons()(5)[c(2)],'#7876B1FF'))+
  theme_bw()+
  theme(axis.text= element_blank(), axis.ticks=element_blank())+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  labs(x = "UMAP1", y = "UMAP2",title = '')
dev.off()

#F2B#
epi$big_phase_cilia<-factor(epi$big_phase_cilia,levels = rev(c('MMP7 epithelium','Glandular epithelium','Ciliated epithelium')))
pdf('F2B, epi_dotplot(EU+Cilia).pdf',width=8,height=3)
DotPlot(epi,features = c('MMP7','IHH','EMID1','NPAS3','CPM',
                         'MT1E','SPP1','C2CD4A','PAEP','CXCL14',
                         'PIFO','FOXJ1','AGR3','SPAG17','SNTN'),
        group.by = 'big_phase_cilia',dot.scale = 6)+
  scale_colour_distiller(palette = "Reds",direction = 1)+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
dev.off()

#Str#
str<-readRDS("../str_removeOvary_EU.rds")
meta_str_old<-data.frame(str@meta.data)
str.data<-data.frame(str@assays$RNA@counts,check.names = F)
str<-CreateSeuratObject(counts = str.data,Project="str",min.cells = 3,min.features = 500)
meta_str_new<-data.frame(str@meta.data)
meta_str_new<-cbind(meta_str_new,meta_str_old[,c(4:9)])
str@meta.data<-meta_str_new
str
str<-NormalizeData(object = str, normalization.method = 'LogNormalize',scale.factor = 1e6)
str<-FindVariableFeatures(str,selection.method = 'vst',nfeatures = 3000)
str.genes<-rownames(str)
str<-ScaleData(str,features = str.genes)
str<-RunPCA(str,features = VariableFeatures(object = str))

str<-JackStraw(str, num.replicate = 100,dims = 20)
str<-ScoreJackStraw(str,dims = 1:20)
pdf("str.pc.pdf",width = 10,height = 6)
JackStrawPlot(str,dims = 1:20)
ElbowPlot(str,ndims = 20)
dev.off()

str<-FindNeighbors(str,dims = 1:9)
str<-FindClusters(str, resolution = 0.8)
str<-RunUMAP(str,dims = 1:9)

#F2C str celltype#
str$celltype<-''
str$celltype[str$phase%in%c('early_proliferative','late_proliferative')]<-'Endometrial stroma'
str$celltype[str$phase%in%c('early_secretory','mid_secretory')]<-'Decidualized stroma'
str$celltype<-factor(str$celltype,levels = c('Endometrial stroma','Decidualized stroma'))
pdf('F2C, str_celltype.pdf',width=5.7,height=4)
DimPlot(str, reduction = "umap", pt.size = 1,label = F,group.by = 'celltype',cols = c('#CD9B1D',pal_simpsons()(5)[c(2)]))+
  labs(x = "UMAP_1", y = "UMAP_2",title = '')+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

#F2D dotplot endo deci marker#
str$celltype<-factor(str$celltype,levels = rev(c('Endometrial stroma','Decidualized stroma')))
features<-list('Proliferative phase'=c('MMP11','SFRP1','CRABP2'),'Secretory phase'=c('S100A4','DKK1','CFD'))
pdf('F2D, dotplot_endo_deci_sepa.pdf',width=6,height=2)
DotPlot(str,features = features,group.by = 'celltype',dot.scale = 6)+
  theme_bw()+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(strip.background = element_blank())+
  scale_colour_distiller(palette='Reds',direction = 1)+
  theme(axis.text = element_text(color = 'black'),axis.title = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()

epi_EU<-subset(epi,origin!='Ciliated epithelium')
#F2E epi pseudotime#
counts<-epi_EU@assays$RNA@counts
sce1<-SingleCellExperiment(assays = List(counts = counts))
genes_quake<-read_xlsx('G:/endometriosis/slingshot/quake_supp_table.xlsx')
geneFilter<-genes_quake$Gene[genes_quake$Cell_type=='Unciliated epithelia']
geneFilter<-intersect(geneFilter,row.names(epi_EU@assays$RNA@data))
sce1<-sce1[geneFilter, ]
#normalize data
FQnorm<-function(counts){
  rk<-apply(counts,2,rank,ties.method='min')
  counts.sort<-apply(counts,2,sort)
  refdist<-apply(counts.sort,1,median)
  norm<-apply(rk,2,function(r){ refdist[r] })
  rownames(norm)<-rownames(counts)
  return(norm)
}
assays(sce1)$norm<-FQnorm(assays(sce1)$counts)
#dimension reduction PCA
pca<-prcomp(t(log1p(assays(sce1)$norm)), scale. = FALSE)
rd1<-pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
#dimension reduction UMAP
rd2<-uwot::umap(t(log1p(assays(sce1)$norm)))
colnames(rd2)<-c('UMAP1', 'UMAP2')
plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sce1)<-SimpleList(PCA = rd1, UMAP = rd2)
sce1<-slingshot(sce1, reducedDim = 'UMAP')
sce1
sce1$phase<-epi_EU$phase
sce1$patient<-epi_EU$patient
sce1$origin<-epi_EU$origin

pseu1<-data.frame(colnames(sce1),sce1$patient,sce1$phase,sce1$origin,sce1$slingPseudotime_1)
colnames(pseu1)=c('Identity','patient','phase','origin','Pseudotime')
#pseu_phase#
pseu1$pseu_phase=''
pseu1$pseu_phase[which(pseu1$Pseudotime<4)]='early_proliferative'
pseu1$pseu_phase[which(pseu1$Pseudotime>=10 & pseu1$Pseudotime<=12.9)]='late_proliferative'
pseu1$pseu_phase[which(pseu1$Pseudotime>=12.9 & pseu1$Pseudotime<=15.2)]='early_secretory'
pseu1$pseu_phase[which(pseu1$Pseudotime>=26)]='mid_secretory'

pseu1$Pseudotime[which(pseu1$Pseudotime>=10)]=pseu1$Pseudotime[which(pseu1$Pseudotime>10)]-6
pseu1$Pseudotime[which(pseu1$Pseudotime>=20)]=pseu1$Pseudotime[which(pseu1$Pseudotime>20)]-10
pseu1<-pseu1[order(pseu1$Pseudotime,decreasing = F),]

cell_remove<-c('P23_P_EU_EPI_SC55','P23_P_EU_EPI_SC34','P23_P_EU_EPI_SC70','P19_N_EU_EP_SC85',
               'P14_N_EU_M_SC78','P27_N_EU_EPI_SC26','P27_N_EU_EPI_SC87','P27_N_EU_EPI_SC92',
               'P27_N_EU_EPI_SC2','P16_N_EU_EP_SC80','P14_N_EU_M_SC81','P28_N_EU_STR_SC79',
               'P23_P_EU_EPI_SC28','P28_N_EU_STR_SC70','P28_N_EU_STR_SC80','P18_P_EU_EP_SC31',
               'P28_N_EU_STR_SC25','P28_N_EU_STR_SC38','P18_P_EU_EP_SC46','P18_P_EU_EP_SC1',
               'P18_P_EU_EP_SC40','P36_P_EU_EPI_SC92','P28_N_EU_EPI_SC31','P28_N_EU_STR_SC57',
               'P16_N_EU_IMM_SC32','P18_P_EU_EP_SC49','P19_N_EU_EP_SC96')
sce1<-sce1[,!((colnames(sce1))%in%cell_remove)]
epi_EU<-epi_EU[,!((colnames(epi_EU))%in%cell_remove)]
pseu1<-pseu1[!((rownames(pseu1))%in%cell_remove),]
write.table(pseu1,'EPI_eu_pseudotime.txt',quote = F,sep = '\t')

pdf('F2E-1, EPI_eu_slingshot_phase.pdf',width = 5.5,height = 6)
plot(reducedDims(sce1)$UMAP, col = c("#B09C85FF","#E64B35FF","#00A087FF","#8491B4FF")[as.factor(sce1$phase)], pch=16, bty="o",col.axis='black',col.lab='black',asp=1,
     tck=-0.02,cex=1,xaxt="n", yaxt="n",main = "", xlab = "UMAP1",  ylab = "UMAP2")
lines(SlingshotDataSet(sce1), lwd=1, col='black',lty=2)
dev.off()

pdf('F2E-2, EPI_eu_slingshot_origin.pdf',width = 5.5,height = 6)
plot(reducedDims(sce1)$UMAP, col = pal_simpsons()(14)[c(13,7)][as.factor(sce1$origin)], pch=16,bty="o",col.axis='black',col.lab='black',asp=1,
     tck=-0.02,cex=1,xaxt="n", yaxt="n",main = "", xlab = "UMAP1",  ylab = "UMAP2")
lines(SlingshotDataSet(sce1), lwd=1, lty=2,col='black')
dev.off()

#F2F str pseudotime#
counts<-str@assays$RNA@counts
sce2<-SingleCellExperiment(assays = List(counts = counts))
genes_quake<-read_xlsx('G:/endometriosis/slingshot/quake_supp_table.xlsx')
geneFilter<-genes_quake$Gene[genes_quake$Cell_type=='Stromal fibroblasts']
geneFilter<-intersect(geneFilter,row.names(str@assays$RNA@data))
sce2<-sce2[geneFilter, ]
sce2
#normalize data
FQnorm<-function(counts){
  rk<-apply(counts,2,rank,ties.method='min')
  counts.sort<-apply(counts,2,sort)
  refdist<-apply(counts.sort,1,median)
  norm<-apply(rk,2,function(r){ refdist[r] })
  rownames(norm)<-rownames(counts)
  return(norm)
}
assays(sce2)$norm<-FQnorm(assays(sce2)$counts)
#dimension reduction PCA
pca<-prcomp(t(log1p(assays(sce2)$norm)), scale. = FALSE)
rd1<-pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
#dimension reduction UMAP
rd2<-uwot::umap(t(log1p(assays(sce2)$norm)))
colnames(rd2)<-c('UMAP1', 'UMAP2')
plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sce2)<-SimpleList(PCA = rd1, UMAP = rd2)
sce2<-slingshot(sce2, reducedDim = 'UMAP')
sce2
sce2$phase<-str$phase
sce2$patient<-str$patient
sce2$origin<-str$origin

pseu2<-data.frame(colnames(sce2),sce2$patient,sce2$phase,sce2$origin,sce2$slingPseudotime_1)
colnames(pseu)=c('Identity','patient','phase','origin','Pseudotime')
#pseu_phase#
pseu2$pseu_phase=''
pseu2$pseu_phase[which(pseu2$Pseudotime<8.65)]='early_proliferative'
pseu2$pseu_phase[which(pseu2$Pseudotime>=8.65 & pseu2$Pseudotime<=14.17)]='late_proliferative'
pseu2$pseu_phase[which(pseu2$Pseudotime>=14.17 & pseu2$Pseudotime<=14.82)]='early_secretory'
pseu2$pseu_phase[which(pseu2$Pseudotime>=14.82)]='early_proretory'
pseu2<-pseu2[order(pseu2$Pseudotime,decreasing = F),]

cell_remove<-c('P29_N_EU_STR_SC44','P29_N_EU_IMM_noCD4_SC68','P14_N_EU_STR_SC15','P31_N_EU_IMM_CD4_SC36',
               'P23_P_EU_EPI_SC8','P23_P_EU_IMM_1_SC34','P23_P_EU_STR_SC79','P23_P_EU_STR_SC13',
               'P23_P_EU_IMM_1_SC47','P23_P_EU_IMMBIG_SC95.1','P23_P_EU_STR_SC90','P17_N_EU_STR_SC24',
               'P23_P_EU_IMM_1_SC53','P23_P_EU_IMM_1_SC31','P23_P_EU_STR_SC45','P9_P_EU_SOME_2_SC7',
               'P18_P_EU_STR_SC18','P14_N_EU_STR_SC7','P18_P_EU_STR_SC42','P14_N_EU_STR_SC3',
               'P27_N_EU_STR_SC9','P14_N_EU_STR_SC11','P27_N_EU_STR_SC8','P14_N_EU_STR_SC18',
               'P14_N_EU_STR_SC30','P14_N_EU_STR_SC31','P16_N_EU_STR_SC63','P32_P_EU_STR_SC30',
               'P16_N_EU_STR_SC91.1','P14_N_EU_STR_SC20','P27_N_EU_STR_SC30','P36_P_EU_STR_SC55',
               'P27_N_EU_STR_SC27','P32_P_EU_STR_SC22','P27_N_EU_EPI_SC59','P28_N_EU_STR_SC46',
               'P28_N_EU_STR_SC55','P18_P_EU_STR_SC6','P28_N_EU_STR_SC56','P28_N_EU_STR_SC72',
               'P28_N_EU_STR_SC18','P28_N_EU_STR_SC17','P28_N_EU_STR_SC86','P28_N_EU_STR_SC93',
               'P28_N_EU_STR_SC15','P28_N_EU_STR_SC78','P28_N_EU_EPI_SC79','P28_N_EU_IMM_SC9',
               'P28_N_EU_STR_SC77','P28_N_EU_STR_SC35','P28_N_EU_STR_SC2','P28_N_EU_IMM_SC37',
               'P18_P_EU_STR_SC19','P28_N_EU_STR_SC31','P28_N_EU_STR_SC33','P28_N_EU_STR_SC61',
               'P28_N_EU_STR_SC81','P28_N_EU_STR_SC27','P28_N_EU_STR_SC45','P28_N_EU_IMM_SC45',
               'P28_N_EU_STR_SC42','P18_P_EU_STR_SC22','P18_P_EU_STR_SC43','P28_N_EU_STR_SC64',
               'P28_N_EU_STR_SC24','P18_P_EU_STR_SC49','P18_P_EU_STR_SC10','P18_P_EU_STR_SC51',
               'P28_N_EU_STR_SC40','P18_P_EU_STR_SC9','P18_P_EU_STR_SC36','P28_N_EU_STR_SC19',
               'P18_P_EU_STR_SC5','P18_P_EU_STR_SC30','P28_N_EU_STR_SC6','P28_N_EU_STR_SC54',
               'P18_P_EU_STR_SC53','P28_N_EU_STR_SC37','P28_N_EU_STR_SC36','P18_P_EU_STR_SC32',
               'P18_P_EU_STR_SC31','P18_P_EU_STR_SC15','P28_N_EU_STR_SC22','P28_N_EU_STR_SC94',
               'P28_N_EU_STR_SC1','P18_P_EU_STR_SC40','P18_P_EU_STR_SC56','P18_P_EU_STR_SC38',
               'P18_P_EU_STR_SC47','P18_P_EU_STR_SC45','P18_P_EU_STR_SC13','P18_P_EU_STR_SC8',
               'P32_P_EU_STR_SC34','P18_P_EU_STR_SC41','P28_N_EU_STR_SC90','P18_P_EU_STR_SC37',
               'P28_N_EU_STR_SC71','P3_P_EU_STR_2_SC12','P9_P_EU_1_SC7','P20_P_EU_IMM_SC54')
cell_remove1<-c(rownames(y[y$UMAP1<=-5 & y$UMAP2>=2,]),
                rownames(y[y$UMAP1>=-5 & y$UMAP1<=0 & y$UMAP2<=-3,]),
                rownames(y[y$UMAP1>=-5 & y$UMAP1<=0 & y$UMAP2>=-1 & y$UMAP2<=1 & y$pseu_phase=='late_proliferative',]))

sce2<-sce2[,!((colnames(sce2))%in%cell_remove)]
str<-str[,!((colnames(str))%in%cell_remove)]
pseu2<-pseu2[!((rownames(pseu2))%in%cell_remove),]

sce2<-sce2[,!((colnames(sce2))%in%cell_remove1)]
str<-str[,!((colnames(str))%in%cell_remove1)]
pseu2<-pseu2[!((rownames(pseu2))%in%cell_remove1),]
write.table(pseu2,'str_eu_pseudotime.txt',quote = F,sep = '\t')

pdf('F2F-1, str_eu_slingshot_phase.pdf',width = 5.5,height = 6)
plot(reducedDims(sce2)$UMAP, col = c("#B09C85FF","#E64B35FF","#00A087FF","#8491B4FF")[as.factor(sce2$phase)], pch=16, bty="o",col.axis='black',col.lab='black',asp=1,
     tck=-0.02,cex=1,xaxt="n", yaxt="n",main = "", xlab = "UMAP1",  ylab = "UMAP2")
lines(SlingshotDataSet(sce2), lwd=1, col='black',lty=2)
dev.off()

pdf('F2F-2, str_eu_slingshot_origin.pdf',width = 5.5,height = 6)
plot(reducedDims(sce2)$UMAP, col = pal_simpsons()(14)[c(13,7)][as.factor(sce2$origin)], pch=16,bty="o",col.axis='black',col.lab='black',asp=1,
     tck=-0.02,cex=1,xaxt="n", yaxt="n",main = "", xlab = "UMAP1",  ylab = "UMAP2")
lines(SlingshotDataSet(sce2), lwd=1, lty=2,col='black')
dev.off()

#F2G#
a1<-paste0('P',1:40)
a2<-c('late_proliferative','proliferative','early_proliferative','secretory','late_proliferative','secretory','mid_secretory','secretory','late_proliferative',
     'early_proliferative','early_proliferative','mid_secretory','late_proliferative','early_secretory','late_proliferative','early_secretory','early_secretory','mid_secretory',
     'early_proliferative','early_proliferative','early_proliferative','early_proliferative','late_proliferative','late_proliferative','secretory','mid_secretory','early_secretory','late_proliferative',
     'mid_secretory','late_proliferative','mid_secretory','mid_secretory','mid_secretory','early_proliferative','early_secretory','mid_secretory',
     NA,'mid_secretory','late_proliferative','early_proliferative')
a3<-c('proliferative','proliferative','proliferative','secretory','proliferative','secretory','secretory','secretory','proliferative',
     'proliferative','proliferative','secretory','proliferative','secretory','proliferative','secretory','secretory','secretory',
     'proliferative','proliferative','proliferative','proliferative','proliferative','proliferative','secretory','secretory','secretory','proliferative',
     'secretory','proliferative','secretory','secretory','secretory','proliferative','secretory','secretory',
     NA,'secretory','proliferative','proliferative')
phase_he<-data.frame(a1,a2,a3)
phase_he<-na.omit(phase_he)
epi_p<-data.frame(table(pseu1$patient))
str_p<-data.frame(table(pseu2$patient))
p_union<-union(epi_p$Var1,str_p$Var1)
phase_he_p<-phase_he[phase_he$a1%in%p_union,]
phase_he_p2<-phase_he_p
row.names(phase_he_p2)<-1:nrow(phase_he_p2)
phase_he_p2[2,2]<-'early_proliferative'
phase_he_p2[16,2]<-'early_secretory'
phase_he_p2[19,2]<-'late_proliferative'

epi_mean<-c()
epi_lower<-c()
epi_upper<-c()
for (i in p_union) {
  x=epi[epi$patient==i,]
  mean=mean(x$Pseudotime)
  lower=min(x$Pseudotime)
  upper=max(x$Pseudotime)
  epi_mean=c(epi_mean,mean)
  epi_lower=c(epi_lower,lower)
  epi_upper=c(epi_upper,upper)
}
str_mean<-c()
str_lower<-c()
str_upper<-c()
for (i in p_union) {
  x=str[str$patient==i,]
  mean=mean(x$Pseudotime)
  lower=min(x$Pseudotime)
  upper=max(x$Pseudotime)
  str_mean=c(str_mean,mean)
  str_lower=c(str_lower,lower)
  str_upper=c(str_upper,upper)
}
row.names(phase_he_p2)<-phase_he_p2$a1
phase_he_p2<-phase_he_p2[p_union,]
pseudotime<-data.frame(p_union,epi_mean,epi_lower,epi_upper,str_mean,str_lower,str_upper,phase_he_p2$a2)
pseudotime[3,5:7]<-NA
pseudotime[27:28,2:4]<-NA
pseudotime<-pseudotime[order(pseudotime$epi_mean),]
table(pseudotime$phase_he_p2.a2)
pseudotime$phase_he_p2.a2<-factor(pseudotime$phase_he_p2.a2,levels = c('early_proliferative','late_proliferative','early_secretory','mid_secretory'))
pseudotime<-pseudotime[order(pseudotime$phase_he_p2.a2),]
pseudotime<-pseudotime[c(1:18,28,19:27),]
pseudotime<-pseudotime[order(pseudotime$epi_mean),]
pseudotime<-pseudotime[order(pseudotime$phase_he_p2.a2),]
pseudotime$p_union<-factor(pseudotime$p_union,levels = pseudotime$p_union)

pseudotime$big_phase<-sapply(strsplit(as.character(pseudotime$phase_he_p2.a2),'_'),'[[',2)
pseudotime<-pseudotime[order(pseudotime$epi_mean),]
pseudotime<-pseudotime[order(pseudotime$big_phase),]
pseudotime$p_union<-factor(pseudotime$p_union,levels = pseudotime$p_union)

pdf('F2G, pseudotime_menstrual_cycle.pdf',height = 6,width = 8)
ggplot() + 
  geom_errorbar(data=pseudotime, mapping=aes(x=p_union, ymin=epi_lower, ymax=epi_upper), width=0.5, size=1, color="#91331FFF") + 
  geom_errorbar(data=pseudotime, mapping=aes(x=p_union, ymin=str_lower, ymax=str_upper), width=0.5, size=1, color="#709AE1FF") + 
  geom_point(data=pseudotime, mapping=aes(x=p_union, y=epi_mean), size=4, shape=21, fill="#91331FFF",color='#91331FFF') +
  geom_point(data=pseudotime, mapping=aes(x=p_union, y=str_mean), size=4, shape=21, fill="#709AE1FF",color='#709AE1FF')+
  theme_classic()+
  xlab(label = 'Patient')+
  ylab(label = 'Pseudotime')+
  geom_rect(aes(xmin=0,xmax=13.5,ymin=-Inf,ymax=Inf),fill="#B09C85FF",alpha=0.2)+
  geom_rect(aes(xmin=13.5,xmax=28.5,ymin=-Inf,ymax=Inf),fill='#E64B35FF',alpha=0.2)+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

