setwd('G:/endometriosis/Analysis/p40/EPI/remove_cilia/all/')
library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(ggsci)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(viridis)

epi<-readRDS('./epi_removeCilia.rds')
cell_remove<-c('P23_P_EU_EPI_SC55','P23_P_EU_EPI_SC34','P23_P_EU_EPI_SC70','P19_N_EU_EP_SC85',
               'P14_N_EU_M_SC78','P27_N_EU_EPI_SC26','P27_N_EU_EPI_SC87','P27_N_EU_EPI_SC92',
               'P27_N_EU_EPI_SC2','P16_N_EU_EP_SC80','P14_N_EU_M_SC81','P28_N_EU_STR_SC79',
               'P23_P_EU_EPI_SC28','P28_N_EU_STR_SC70','P28_N_EU_STR_SC80','P18_P_EU_EP_SC31',
               'P28_N_EU_STR_SC25','P28_N_EU_STR_SC38','P18_P_EU_EP_SC46','P18_P_EU_EP_SC1',
               'P18_P_EU_EP_SC40','P36_P_EU_EPI_SC92','P28_N_EU_EPI_SC31','P28_N_EU_STR_SC57',
               'P16_N_EU_IMM_SC32','P18_P_EU_EP_SC49','P19_N_EU_EP_SC96','P2_N_EU_SOME_2_SC28')
epi<-epi[,!((colnames(epi))%in%cell_remove)]
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

epi<-FindNeighbors(epi,dims = 1:5)
epi<-FindClusters(epi, resolution = 0.4)
epi<-RunUMAP(epi,dims = 1:5,min.dist = 0.2)

#big phase#
epi$phase2<-''
epi$phase2[which(epi$phase%in%c('early_proliferative','late_proliferative'))]<-'proliferative'
epi$phase2[which(epi$phase%in%c('early_secretory','mid_secretory'))]<-'secretory'
##state_loci_phase2##
epi$state_loci_phase2<-''
epi$state_loci_phase2<-paste0(epi$origin,'_',epi$phase2)

#color#
pal<-c('#BC3C29','#0072B5','#E18727','#20854E','#7876B1','#6F99AD')

#F4A#
pdf('F4A, EPI_location(N+EU+EC).pdf',width=4.8,height=4)
DimPlot(epi, reduction = "umap", pt.size = 1, group.by = 'origin',label = F,
        cols = pal_simpsons()(14)[c(4,13,7)])+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

#F4B#
patient<-subset(epi,patient%in%c('P15','P23','P36'))
quake<-read_xlsx('G:/endometriosis/slingshot/quake_supp_table.xlsx')
quake<-quake[quake$Cell_type=='Unciliated epithelia',]
quake<-quake[quake$Gene%in%row.names(epi@assays$RNA@data),]
rownames(quake)<-quake$Gene
#top100 quake markers#
quake_remove_late_sec<-quake[1:1279,]
quake_remove_late_sec %>% group_by(Peak_phase) %>% top_n(n = 100, wt = Peak_height)-> quake_top100
patient_quake_top100_phase_data<-data.frame(patient@assays$RNA@data,check.names = F,check.rows = F)
patient_quake_top100_phase_data<-patient_quake_top100_phase_data[quake_top100$Gene[1:300],]
patient_meta<-data.frame(patient@meta.data)
patient_quake_top100_phase_data<-scale(patient_quake_top100_phase_data)
patient_quake_top100_phase_data<-as.data.frame(t(patient_quake_top100_phase_data))
patient_quake_top100_phase_data$origin<-patient_meta$origin
patient_quake_top100_phase_data$phase<-patient_meta$phase
patient_quake_top100_phase_data<-arrange(patient_quake_top100_phase_data,patient_quake_top100_phase_data$phase)
#early_pro 88,early_sec 73,late_pro 61
patient_quake_top100_phase_data<-patient_quake_top100_phase_data[c(1:88,162:222,89:161),]
plot1_col<-data.frame(row.names=rownames(patient_quake_top100_phase_data),origin=patient_quake_top100_phase_data$origin,
                      phase=patient_quake_top100_phase_data$phase,patient_origin=patient_quake_top100_phase_data$patient_origin)
plot1_row <- data.frame(row.names=rownames(t(patient_quake_top100_phase_data[,1:300])),
                        GeneClass = factor(rep(c("Early_proliferative","Late_proliferative",'Early_secretory'), c(100,100,100))))
plot1_color<-list(origin=c('Eutopic'='#C80813FF','Ectopic'='#D2AF81FF'),
                  phase=c("early_proliferative"="#B09C85FF","late_proliferative"="#00A087FF",'early_secretory'="#E64B35FF"),
                  GeneClass=c("Early_proliferative"="#B09C85FF","Late_proliferative"="#00A087FF",'Early_secretory'="#E64B35FF"))

pdf("F4B, Epi_pheatmap_phase_related_top100genes(EU+EC).pdf",width=7,height=6)
pheatmap(t(patient_quake_top100_phase_data[,1:300]),cluster_cols = F,cluster_rows=F,show_colnames = F,show_rownames =F,
         scale='row',
         color =colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(50),
         annotation_col = plot3_col,annotation_row = plot3_row,
         annotation_colors =  plot3_color,
         breaks = unique(c(seq(-2,2,length=50))),
         main = "Phase related genes")
dev.off()

#F4C#
#proliferative phase score
phase12_genes<-read.table('G:/endometriosis/gene_list/epi_quakephase12.txt',header=T)
phase12_genes$Gene<-as.character(phase12_genes$Gene)
epi=AddModuleScore(object = epi,features = list(phase12_genes$Gene),ctrl = 5,name = 'phase12_genes')
#secretory phase score
phase34_genes<-read.table('G:/endometriosis/gene_list/epi_quakephase34.txt',header=T)
phase34_genes$MyList<-as.character(phase34_genes$MyList)
epi=AddModuleScore(object = epi,features = list(phase34_genes$MyList),ctrl = 5,name = 'phase34_genes')

epi_pro<-subset(epi,phase_abbr=='Pro')
pdf('F4C-1, epi_pro_phase_score.pdf',width=4,height=5)
VlnPlot(epi_pro,'phase12_genes1',group.by = 'state_loci_phase2',cols = pal,pt.size = 0)+
  labs(title = '')+theme(legend.position="none")+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5))
dev.off()
epi_sec<-subset(epi,phase_abbr=='Sec')
pdf('F4C-2, epi_sec_phase_score.pdf',width=4,height=5)
VlnPlot(epi_sec,'phase34_genes1',group.by = 'state_loci_phase2',cols = pal[4:6],pt.size = 0)+
  labs(title = '')+theme(legend.position="none")+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5))
dev.off()

#F4D#
p4_genes<-read.table('G:/endometriosis/gene_list/epi_p4genes.txt',header=T)
p4_genes$p4_score<-as.character(p4_genes$p4_score)
heat_gene1<-unique(c('FOXO1','HSD17B2','FKBP5','ERRFI1',
                     p4_genes$p4_score[c(8,12,15,17,19,22,23,27,29:32,34:38)]))
heatdata1<-data.frame(epi@assays$RNA@data,check.names = F,check.rows = F)
heatdata1<-heatdata1[heat_gene1,]
heatdata1<-as.data.frame(t(heatdata1))
heatdata1$state_loci_phase2<-meta$state_loci_phase2
heatdata1<-melt(heatdata1,measure.vars = colnames(heatdata1[1:22]),variable.name = "gene",value.name = "Expression")
heatdata1_mean<-aggregate(heatdata1,by=list(heatdata1$state_loci_phase2,heatdata1$gene),FUN = mean)
heatdata1_mean<-heatdata1_mean[,c(1,2,5)]
colnames(heatdata1_mean)[1:2]<-c('state_loci_phase2','gene')
heatdata1_mean$state_loci_phase2<-factor(heatdata1_mean$state_loci_phase2,levels=rev(c('Normal_proliferative','Eutopic_proliferative','Ectopic_proliferative','Normal_secretory','Eutopic_secretory','Ectopic_secretory')))
pdf("F4D, heatmap_p4.pdf",width=6.5,height=2.5)
ggplot(heatdata1_mean, aes(gene, state_loci_phase2)) + 
  geom_tile(aes(fill = Expression),colour = "white")+ 
  scale_fill_viridis(name="Expression",option = "A") +
  theme_classic()+
  theme(axis.text = element_text(color='black'))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank())+
  labs(x = " ", y = " ", title = "")
dev.off()

#F4E#
epi=AddModuleScore(object = epi,features = list(p4_genes$p4_score),ctrl = 5,name = 'p4_genes')
pdf('F4E, epi_p4_score.pdf',width=4,height=5)
epi$state_loci_phase2<-factor(epi$state_loci_phase2,levels=(c('Normal_proliferative','Eutopic_proliferative','Ectopic_proliferative','Normal_secretory','Eutopic_secretory','Ectopic_secretory')))
VlnPlot(epi,'p4_genes1',group.by = 'state_loci_phase2',cols = pal,pt.size = 0)+
  labs(title = '')+theme(legend.position="none")+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5))
dev.off()

#F4F#
epi_pro<-subset(epi,phase2=='proliferative'&origin%in%c('Eutopic','Ectopic'))
pdf('F4F, stackvlnplot_epi_pro_estrodial_sitimulated_genes(EU+EC).pdf',width=3,height=4)
epi_pro$state_loci_phase2<-factor(epi$state_loci_phase2,levels=(c('Normal_proliferative','Eutopic_proliferative','Ectopic_proliferative')))
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = 'state_loci_phase2',... )  + 
    xlab("") + ylab('') + ggtitle(feature) + theme_classic()+
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = rel(1),color ='black'), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_line(color ='black'),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5,color ='black'))
  return(p)
}
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,color ='black'), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 2)
  return(p)
}
StackedVlnPlot(epi_pro, features = c('GPER1','EGR3','CCN2','CCND1','TP53','ID3'),
               pt.size=0, cols=pal[2:3])
dev.off()

#F4H EC upregulation GO Terms##
P15_GO<-read_xlsx('./DEG/EC/GEO/p15_ec_up_meta.xlsx',sheet = 3)
P15_GO$Description<-factor(P15_GO$Description,levels = rev(P15_GO$Description))
pdf("F4H, P15_ECvsEU_upregulation_GOterms.pdf",width = 10.8,height = 4)
ggplot(P15_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#D2AF81FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-LogP", title = "Go terms of EC vs EU in P15")
dev.off()

#F4I#
heat_gene2<-unique(c('FBLN2','FN1','LAMB3','LAMA4','TIMP2',
                     'NRP1','NRP2','VEGFA','VEGFB',
                     'SOD2','ALDH3B1','GPX3','HMOX1'))
heatdata2<-data.frame(epi@assays$RNA@data,check.names = F,check.rows = F)
heatdata2<-heatdata2[heat_gene2,]
heatdata2<-as.data.frame(t(heatdata2))
heatdata2$state_loci_phase2<-meta$state_loci_phase2
heatdata2<-melt(heatdata2,measure.vars = colnames(heatdata2[1:13]),variable.name = "gene",value.name = "Expression")
heatdata2_mean<-aggregate(heatdata2,by=list(heatdata2$state_loci_phase2,heatdata2$gene),FUN = mean)
heatdata2_mean<-heatdata2_mean[,c(1,2,5)]
colnames(heatdata2_mean)[1:2]<-c('state_loci_phase2','gene')
heatdata2_mean$state_loci_phase2<-factor(heatdata2_mean$state_loci_phase2,levels=rev(c('Normal_proliferative','Eutopic_proliferative','Ectopic_proliferative','Normal_secretory','Eutopic_secretory','Ectopic_secretory')))
pdf("F4I, heatmap(ECM,blood vessel,ox stress).pdf",width=5,height=2.5)
ggplot(heatdata2_mean, aes(gene, state_loci_phase2)) + 
  geom_tile(aes(fill = Expression),colour = "white")+ 
  scale_fill_viridis(name="Expression",option = "A") +
  theme_classic()+
  theme(axis.text = element_text(color='black'))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank())+
  labs(x = " ", y = " ", title = "")
dev.off()

#F4J#
meta<-data.frame(epi@meta.data)
P15_P23_EC<-rownames(meta[meta$patient%in%c('P15','P23') & meta$phase2=='proliferative'& meta$origin%in%c('Ectopic'),])
P15_P23_EU<-rownames(meta[meta$patient%in%c('P15','P23') & meta$phase2=='proliferative'& meta$origin%in%c('Eutopic'),])
P36_EC<-rownames(meta[meta$patient=='P36' & meta$phase2=='secretory'& meta$origin%in%c('Ectopic'),])
P36_EU<-rownames(meta[meta$patient=='P36' & meta$phase2=='secretory'& meta$origin%in%c('Eutopic'),])

P15_P23_ECvsEU<-FindMarkers(epi,ident.1 = P15_P23_EC, ident.2 = P15_P23_EU,only.pos = F,min.pct = 0.25, logfc.threshold = 0)
P15_P23_ECvsEU$gene<-rownames(P15_P23_ECvsEU)
P36_ECvsEU<-FindMarkers(epi,ident.1 = P36_EC, ident.2 = P36_EU,only.pos = F,min.pct = 0.25, logfc.threshold = 0)
P36_ECvsEU$gene<-rownames(P36_ECvsEU)
P15P23_P36_ECvsEU<-merge(P15_P23_ECvsEU,P36_ECvsEU,by = 'gene')
P15P23_P36_ECvsEU<-P15P23_P36_ECvsEU[,c(1,3,8)]
colnames(P15P23_P36_ECvsEU)<-c('gene','Proliferative','Secretory')
write.table(P15P23_P36_ECvsEU,file = "./DEG/EC/P15P23_P36_ECvsEU.DEGs.xls", row.names = F, col.names = T, quote = F, sep = "\t")

label_genes<-unique(c('DCN','FBLN2','FN1','LAMB3','LAMA4','SPARC','TIMP2','NRP1','NRP2','VEGFA','VEGFB','NNMT',
                      'ALDH1A3','RDH10','STRA6','GJA1','MGP','RARRES1','RARRES2','RUNX3','CYP1B1','RBP1','RORB','RXRA'))
P15P23_P36_ECvsEU$label_genes<-''
P15P23_P36_ECvsEU$label_genes_sig<-''
P15P23_P36_ECvsEU$label_genes_sig<-'con'
for (i in 1:nrow(P15P23_P36_ECvsEU)) {
  if (P15P23_P36_ECvsEU[i,1]%in%c(label_genes)) {
    P15P23_P36_ECvsEU[i,4]<-P15P23_P36_ECvsEU[i,1]
    P15P23_P36_ECvsEU[i,5]<-'label'
  }
}
P15P23_P36_ECvsEU$label_genes_sig<-as.factor(P15P23_P36_ECvsEU$label_genes_sig)
write.table(P15P23_P36_ECvsEU,file = "./DEG/EC/P15P23_P36_ECvsEU.DEGs.label.xls", row.names = F, col.names = T, quote = F, sep = "\t")

pdf("F4J, Volcano_P15P23_P36_ECvsEU.pdf",width=4.0,height=4)
ggplot(P15P23_P36_ECvsEU,aes(x=Proliferative,y=Secretory,color=label_genes_sig,size=label_genes_sig))+geom_point()+
  theme_classic()+
  scale_color_manual(name='',values =c("label"="#91331FFF","con"="grey"))+
  scale_size_manual(values = c(0.5,1))+
  labs(x = 'EC vs EU in proliferative phase (log2)', y = 'EC vs EU in secretory phase (log2)', title = "")+
  geom_text_repel(data = P15P23_P36_ECvsEU, aes(x = Proliferative, 
                                                y = Secretory, 
                                                label = label_genes),
                  size = 2,box.padding = unit(0.1, "lines"),#shape=18,
                  point.padding = unit(0.3, "lines"), 
                  max.overlaps = 500,
                  segment.color = "black", 
                  color="black",
                  show.legend = FALSE)+
  theme(legend.position = 'none')
dev.off() 

