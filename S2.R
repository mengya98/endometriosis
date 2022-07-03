setwd('G:/endometriosis/Analysis/p40/EPI/EU')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(ggsci)
library(readxl)

epi<-readRDS('./epi_removeEC.rds')
#S2A#
epi$bigphase<-''
epi$bigphase[epi$phase%in%c('early_proliferative','late_proliferative')]<-'Proliferative'
epi$bigphase[epi$phase%in%c('early_secretory','mid_secretory')]<-'Secretory'
epi$bigphase<-factor(epi$bigphase,levels = c('Proliferative','Secretory'))
pdf('S2A, epi_menstrual_cycle.pdf',width=5.1,height=4)
DimPlot(epi, reduction = "umap", pt.size = 1,label = F,group.by = 'bigphase',
        #cols = c(pal_simpsons()(15)[c(14,2)],'#7876B1FF'),
        cols = c("#B09C85FF","#E64B35FF"))+
  theme_bw()+
  theme(axis.text= element_blank(), axis.ticks=element_blank())+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  labs(x = "UMAP1", y = "UMAP2",title = '')
dev.off()
remove(epi)

str<-readRDS('str_removeEC.rds')
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

str<-FindNeighbors(str,dims = 1:10)
str<-FindClusters(str, resolution = 0.8)
str<-RunUMAP(str,dims = 1:10)

#ovary cluster#
str$origin[str$seurat_clusters%in%c(0,3,6,10,11,14,16,17)]='Ovary'
#remove the cells that pseudotime inconsistent with real time#
cell_remove<-c('P29_N_EU_STR_SC44','P29_N_EU_IMM_noCD4_SC68','P14_N_EU_STR_SC15','P31_N_EU_IMM_CD4_SC36',
               'P23_P_EU_EPI_SC8','P23_P_EU_IMM_1_SC34','P23_P_EU_STR_SC79','P23_P_EU_STR_SC13',
               'P23_P_EU_IMM_1_SC47','P23_P_EU_IMMBIG_SC95.1','P23_P_EU_STR_SC90','P17_N_EU_STR_SC24',
               'P23_P_EU_IMM_1_SC53','P23_P_EU_IMM_1_SC31','P23_P_EU_STR_SC45',
               'P9_P_EU_1_SC7','P9_P_EU_SOME_2_SC7',
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
               'P28_N_EU_STR_SC71','P3_P_EU_STR_2_SC12','P20_P_EU_IMM_SC54',
               "P34_P_EU_EPI_SC40","P34_P_EU_EPI_SC44","P34_P_EU_STR_SC49","P34_P_EU_STR_SC51", 
               "P34_P_EU_STR_SC52","P34_P_EU_STR_SC53","P34_P_EU_STR_SC54","P34_P_EU_STR_SC55",  
               "P34_P_EU_STR_SC56","P34_P_EU_STR_SC57","P34_P_EU_STR_SC58","P34_P_EU_STR_SC59", 
               "P34_P_EU_STR_SC60","P34_P_EU_STR_SC61","P34_P_EU_STR_SC62","P34_P_EU_STR_SC63", 
               "P34_P_EU_STR_SC65","P34_P_EU_STR_SC67","P34_P_EU_STR_SC68","P34_P_EU_STR_SC69", 
               "P34_P_EU_STR_SC71","P34_P_EU_STR_SC72","P34_P_EU_STR_SC73","P34_P_EU_STR_SC74", 
               "P34_P_EU_STR_SC75","P34_P_EU_STR_SC76","P34_P_EU_STR_SC77","P34_P_EU_STR_SC78", 
               "P34_P_EU_STR_SC79","P34_P_EU_STR_SC80","P34_P_EU_STR_SC82","P34_P_EU_STR_SC84", 
               "P34_P_EU_STR_SC86","P34_P_EU_STR_SC87","P34_P_EU_STR_SC88","P34_P_EU_STR_SC89", 
               "P34_P_EU_STR_SC90","P34_P_EU_STR_SC91","P34_P_EU_STR_SC92","P34_P_EU_STR_SC93",
               "P34_P_EU_STR_SC94","P34_P_EU_STR_SC96","P11_P_EU_1_SC1","P11_P_EU_1_SC11",    
               "P11_P_EU_1_SC14","P11_P_EU_1_SC16","P11_P_EU_1_SC17","P11_P_EU_1_SC21",    
               "P11_P_EU_1_SC25","P11_P_EU_1_SC26","P11_P_EU_1_SC27","P11_P_EU_1_SC28",    
               "P11_P_EU_1_SC29","P11_P_EU_1_SC3","P11_P_EU_1_SC30","P11_P_EU_1_SC34",    
               "P11_P_EU_1_SC35","P11_P_EU_1_SC37","P11_P_EU_1_SC38","P11_P_EU_1_SC40",   
               "P11_P_EU_1_SC42","P11_P_EU_1_SC43","P11_P_EU_1_SC47","P11_P_EU_1_SC5",     
               "P11_P_EU_1_SC6","P11_P_EU_1_SC7","P23_P_EU_EPI_SC79","P23_P_EU_EPI_SC84",  
               "P23_P_EU_IMM_1_SC78","P23_P_EU_STR_SC34","P23_P_EU_STR_SC51","P23_P_EU_STR_SC53")
str<-str[,!((colnames(str))%in%cell_remove)]

#S2B#
pdf('S2B, str_ovary_location.pdf',width=4.8,height=4)
DimPlot(str, reduction = "umap", pt.size = 1,label = F,group.by = 'origin',cols=c(pal_simpsons()(16)[c(13,7,5)]))+#,'#7876B1FF'
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

##ovary marker##
umap<-data.frame(str@reductions$umap@cell.embeddings)
umap$MME<-str@assays$RNA@data['MME',]
umap$ESR1<-str@assays$RNA@data['ESR1',]
umap$ARX<-str@assays$RNA@data['ARX',]
umap$STAR<-str@assays$RNA@data['STAR',]
p1<-ggplot(umap, aes(x=UMAP_1, y=UMAP_2,color= MME))+
  geom_point(size=1) + 
  theme_classic()+
  theme(axis.title.y = element_text(size = 25, color = "black"))+
  theme(axis.line.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+
  labs(x = "UMAP1", y = "UMAP2", title = "MME")+
  theme(plot.title = element_text(hjust = 0.5,size = 25))+
  scale_colour_distiller(palette='Reds',direction = 1,limit=c(-0.1,8.3),breaks=c(0,2,4,6,8))+
  theme(legend.position="none")
p2<-ggplot(umap, aes(x=UMAP_1, y=UMAP_2,color= ESR1))+
  geom_point(size=1) + 
  theme_classic()+
  theme(axis.title.x = element_text(size = 25, color = "black"))+
  theme(axis.title.y = element_text(size = 25, color = "black"))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+
  labs(x = "UMAP1", y = "UMAP2", title = "ESR1")+
  theme(plot.title = element_text(hjust = 0.5,size = 25))+
  scale_colour_distiller(palette='Reds',direction = 1,limit=c(-0.1,8.3),breaks=c(0,2,4,6,8))+
  theme(legend.position="none")
p3<-ggplot(umap, aes(x=UMAP_1, y=UMAP_2,color= ARX))+
  geom_point(size=1) + 
  theme_classic()+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())+
  labs(x = "UMAP1", y = "UMAP2", title = "ARX",col='Expression')+
  theme(plot.title = element_text(hjust = 0.5,size = 25))+
  scale_colour_distiller(palette='Reds',direction = 1,limit=c(-0.1,8.3),breaks=c(0,2,4,6,8))#+
p4<-ggplot(umap, aes(x=UMAP_1, y=UMAP_2,color= STAR))+
  geom_point(size=1) + 
  theme_classic()+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+
  theme(axis.line.y = element_blank(),axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.title.y = element_blank())+
  theme(axis.title.x = element_text(size = 25, color = "black"))+
  labs(x = "UMAP1", y = "UMAP2", title = "STAR")+
  theme(plot.title = element_text(hjust = 0.5,size = 25))+
  scale_colour_distiller(palette='Reds',direction = 1,limit=c(-0.1,8.3),breaks=c(0,2,4,6,8))+
  theme(legend.position="none")
#S2C#
pdf('S2C, STR_ovary.marker.pdf',width=8.4,height=8)
patchwork::wrap_plots(p1,p3,p2,p4,ncol = 2)
dev.off()
remove(str)

#str EU#
str_EU<-subset(str,origin!='Ovary')

#S2D menstrual cycle#
str_EU$bigphase<-''
str_EU$bigphase[str_EU$phase%in%c('early_proliferative','late_proliferative')]<-'Proliferative'
str_EU$bigphase[str_EU$phase%in%c('early_secretory','mid_secretory')]<-'Secretory'
str_EU$bigphase<-factor(str_EU$bigphase,levels = c('Proliferative','Secretory'))
pdf('S2D, str_EU_menstrual_cycle.pdf',width=5.0,height=4)
DimPlot(str_EU, reduction = "umap", pt.size = 1,label = F,group.by = 'bigphase',
        cols = c("#B09C85FF","#E64B35FF"))+
  theme_bw()+
  theme(axis.text= element_blank(), axis.ticks=element_blank())+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  labs(x = "UMAP1", y = "UMAP2",title = '')
dev.off()

#S2E epi phase related genes#
quake<-read_xlsx('G:/endometriosis/slingshot/quake_supp_table.xlsx')
quake<-quake[quake$Cell_type=='Unciliated epithelia',]
quake<-quake[quake$Gene%in%row.names(epi@assays$RNA@data),]
rownames(quake)<-quake$Gene
quake_remove_late_sec<-quake[1:1279,]
quake_remove_late_sec %>% group_by(Peak_phase) %>% top_n(n = 100, wt = Peak_height)-> quake_top100
quake_top100_phase_data<-data.frame(epi@assays$RNA@data,check.names = F,check.rows = F)
quake_top100_phase_data<-quake_top100_phase_data[quake_top100$Gene,]
quake_top100_phase_data<-scale(quake_top100_phase_data)
quake_top100_phase_data<-as.data.frame(t(quake_top100_phase_data))
quake_top100_phase_data$origin<-meta$origin
quake_top100_phase_data$phase<-meta$phase
quake_top100_phase_data$patient_origin<-meta$patient_origin
quake_top100_phase_data<-arrange(quake_top100_phase_data,quake_top100_phase_data$phase)
quake_top100_phase_data<-arrange(quake_top100_phase_data,desc(quake_top100_phase_data$origin))

plot1_col<-data.frame(row.names=rownames(quake_top100_phase_data),origin=quake_top100_phase_data$origin,phase=quake_top100_phase_data$phase)
plot1_row <- data.frame(row.names=rownames(t(quake_top100_phase_data[,1:400])),
                        GeneClass = factor(rep(c("Early_proliferative","Late_proliferative",'Early_secretory','Mid_secretory'), c(100,100,100,100))))
plot1_color<-list(origin=c('Normal'='#197EC0FF','Eutopic'='#C80813FF'),
                  phase=c("early_proliferative"="#B09C85FF","late_proliferative"="#00A087FF",'early_secretory'="#E64B35FF",'mid_secretory'="#8491B4FF"),
                  GeneClass=c("Early_proliferative"="#B09C85FF","Late_proliferative"="#00A087FF",'Early_secretory'="#E64B35FF",'Mid_secretory'="#8491B4FF"))

pdf("S2E, Epi_pheatmap_phase_related_top100genes(N+EU).pdf",width=7,height=6)
pheatmap(t(quake_top100_phase_data[,1:400]),cluster_cols = F,cluster_rows=F,show_colnames = F,show_rownames =F,
         scale='row',
         color =colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(50),
         annotation_col = plot1_col,annotation_row = plot1_row,
         annotation_colors =  plot1_color,
         breaks = unique(c(seq(-2,2,length=50))),
         main = "Phase related genes in epithelium")
dev.off()

#S2F str phase related genes#
quake<-read_xlsx('G:/endometriosis/slingshot/quake_supp_table.xlsx')
quake<-quake[quake$Cell_type=='Stromal fibroblasts',]
quake<-quake[quake$Gene%in%row.names(str_EU@assays$RNA@data),]
rownames(quake)<-quake$Gene

quake_phase_data<-data.frame(str_EU@assays$RNA@data,check.names = F,check.rows = F)
quake_phase_data<-quake_phase_data[quake$Gene,]
quake_phase_data<-scale(quake_phase_data)
quake_phase_data<-as.data.frame(t(quake_phase_data))
quake_phase_data$origin<-meta$origin
quake_phase_data$phase<-meta$phase
quake_phase_data<-arrange(quake_phase_data,quake_phase_data$phase)
quake_phase_data<-quake_phase_data[c(1:783,991:1184,784:990,1185:1399),]
quake_phase_data<-arrange(quake_phase_data,desc(quake_phase_data$origin))

plot2_col<-data.frame(row.names=rownames(quake_phase_data),origin=quake_phase_data$origin,phase=quake_phase_data$phase)
plot2_row <- data.frame(row.names=rownames(t(quake_phase_data[,1:518])),
                        GeneClass = factor(rep(c("Early_proliferative","Late_proliferative",'Early_secretory','Mid_secretory'), c(274,75,18,151))))
plot2_color<-list(origin=c('Normal'='#197EC0FF','Eutopic'='#C80813FF'),
                  phase=c("early_proliferative"="#B09C85FF","late_proliferative"="#00A087FF",'early_secretory'="#E64B35FF",'mid_secretory'="#8491B4FF"),
                  GeneClass=c("Early_proliferative"="#B09C85FF","Late_proliferative"="#00A087FF",'Early_secretory'="#E64B35FF",'Mid_secretory'="#8491B4FF"))

pdf("S2F, Str_pheatmap_phase_related_genes(N+EU).pdf",width=7,height=6)
pheatmap(t(quake_phase_data[,1:518]),cluster_cols = F,cluster_rows= F,show_colnames = F,show_rownames = F,
         scale='row',
         color =colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(50),
         annotation_col = plot2_col,annotation_row = plot2_row,
         annotation_colors =  plot2_color,
         breaks = unique(c(seq(-2,2,length=50))),
         main = "Phase related genes in stroma")
dev.off()

#S2G#
epi_data<-data.frame(epi@assays$RNA@data)
epi_data<-data.frame(t(epi_data),check.names = F)
pseu1<-read.table('EPI_eu_pseudotime.txt',header = T)
epi_wnt<-c('LAMC2','MMP2','MMP7')
epi_data_wnt<-epi_data[,colnames(epi_data)%in%(epi_wnt)]
epi_data_wnt$Identity<-rownames(epi_data_wnt)
pseu_wnt<-merge(pseu1,epi_data_wnt,by = 'Identity')

p1<-ggplot(pseu_wnt)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = epi_wnt[1], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt,size=0.1,aes_string(x = 'Pseudotime', y = epi_wnt[1], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_wnt[1])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p2<-ggplot(pseu_wnt)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_wnt[2], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt,size=0.1,aes_string(x = 'Pseudotime', y = epi_wnt[2], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_wnt[2])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p3<-ggplot(pseu_wnt)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_wnt[3], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt,size=0.1,aes_string(x = 'Pseudotime', y = epi_wnt[3], color= 'origin',alpha=0.001))+
  labs(y='Expression',title = epi_wnt[3])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
pdf("S2G, pseu_epi_EU_PvsN_wnt.pdf", width = 9, height = 6)
patchwork::wrap_plots(p1,p2,p3,ncol = 2)
dev.off()

#S2H#
epi_cellcycle<-c('CCN2','CCNB1')
epi_data_cellcycle<-epi_data[,colnames(epi_data)%in%(epi_cellcycle)]
epi_data_cellcycle$Identity<-rownames(epi_data_cellcycle)
pseu_cellcycle<-merge(pseu1,epi_data_cellcycle,by = 'Identity')
p1<-ggplot(pseu_cellcycle)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = epi_cellcycle[1], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_cellcycle,size=0.1,aes_string(x = 'Pseudotime', y = epi_cellcycle[1], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_cellcycle[1])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p2<-ggplot(pseu_cellcycle)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_cellcycle[2], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_cellcycle,size=0.1,aes_string(x = 'Pseudotime', y = epi_cellcycle[2], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_cellcycle[2])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
pdf("S2H, pseu_epi_EU_PvsN_cellcycle.pdf", width = 5, height = 6)
patchwork::wrap_plots(p1,p2,ncol = 1)
dev.off()

#S2I#
epi_cilia<-c('FOXJ1','RFX2','CDHR3','PIFO')
epi_data_cilia<-epi_data[,colnames(epi_data)%in%(epi_cilia)]
epi_data_cilia$Identity<-rownames(epi_data_cilia)
pseu_cilia<-merge(pseu1,epi_data_cilia,by = 'Identity')

p1<-ggplot(pseu_cilia)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = epi_cilia[1], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_cilia,size=0.1,aes_string(x = 'Pseudotime', y = epi_cilia[1], color= 'origin',alpha=0.001))+
  labs(y='Expression',title = epi_cilia[1])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p2<-ggplot(pseu_cilia)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_cilia[2], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_cilia,size=0.1,aes_string(x = 'Pseudotime', y = epi_cilia[2], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_cilia[2])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p3<-ggplot(pseu_cilia)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_cilia[3], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_cilia,size=0.1,aes_string(x = 'Pseudotime', y = epi_cilia[3], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_cilia[3])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p4<-ggplot(pseu_cilia)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_cilia[4], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_cilia,size=0.1,aes_string(x = 'Pseudotime', y = epi_cilia[4], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_cilia[4])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
pdf("S2I, pseu_epi_EU_PvsN_cilia.pdf", width = 9, height = 6)
patchwork::wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()

#S2J#
str_data<-data.frame(str_EU@assays$RNA@data)
str_data<-data.frame(t(str_data),check.names = F)
pseu2<-read.table('STR_eu_pseudotime.txt',header = T)
str_TF<-c('GAB2','FCGR2B','FPR1')
str_data_TF<-str_data[,colnames(str_data)%in%(str_TF)]
str_data_TF$Identity<-rownames(str_data_TF)
pseu_TF<-merge(pseu2,str_data_TF,by = 'Identity')

p1<-ggplot(pseu_TF)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = str_TF[1], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_TF,size=0.1,aes_string(x = 'Pseudotime', y = str_TF[1], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = str_TF[1])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(10,14.17,14.82),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.text = element_text(color = 'black'),axis.title = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))+
  theme(legend.position="none")
p2<-ggplot(pseu_TF)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = str_TF[2], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_TF,size=0.1,aes_string(x = 'Pseudotime', y = str_TF[2], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = str_TF[2])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(10,14.17,14.82),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.text = element_text(color = 'black'),axis.title = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))+
  theme(legend.position="none")
p3<-ggplot(pseu_TF)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = str_TF[3], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_TF,size=0.1,aes_string(x = 'Pseudotime', y = str_TF[3], color= 'origin',alpha=0.001))+
  labs(y='Expression',title = str_TF[3])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(10,14.17,14.82),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.text = element_text(color = 'black'),axis.title = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))+
  theme(legend.position="none")
pdf("S2J, pseu_str_EU_PvsN_TF.pdf", width = 9, height = 6)
patchwork::wrap_plots(p1,p2,p3,ncol = 2)
dev.off()

pal<-c('#BC3C29','#0072B5','#E18727','#20854E','#7876B1','#6F99AD')
str_midsec<-subset(str_EU,phase%in%c('mid_secretory') & origin%in%c('Normal','Eutopic'))
str_midsec$loci_phase<-paste0(str_midsec$origin,'_',str_midsec$phase)
#S2K#
str_rep_posi<-read.table('G:/endometriosis/gene_list/endometirum_receptivity_posmarkers.txt',header=T)
str_rep_posi<-data.frame(gene=str_rep_posi[str_rep_posi$gene%in%c(str.genes),])
str_rep_posi$gene<-as.character(str_rep_posi$gene)
str=AddModuleScore(object = str_midsec,features = list(str_rep_posi$gene),ctrl = 5,name = 'str_rep_posi')
pdf('S2K, str_midsec_receptivity_posi_score.pdf',width=4,height=5)
str_midsec$loci_phase<-factor(str_midsec$loci_phase,levels=(c('Normal_mid_secretory','Eutopic_mid_secretory')))
VlnPlot(str_midsec,'str_rep_posi1',group.by = 'loci_phase',cols = pal[4:5],pt.size = 0)+
  labs(title = 'Receptivity positive signature')+theme(legend.position="none")+
  theme(axis.text.x = element_text(angle=45))+
  theme(axis.text = element_text(color = 'black'))
dev.off()

#S2L#
str_rep_nega<-read.table('G:/endometriosis/gene_list/endometrium_receptivity_negmarkers.txt',header=T)
str_rep_nega<-data.frame(gene=str_rep_nega[str_rep_nega$gene%in%c(str.genes),])
str_rep_nega$gene<-as.character(str_rep_nega$gene)
str=AddModuleScore(object = str_midsec,features = list(str_rep_nega$gene),ctrl = 5,name = 'str_rep_nega')
pdf('S2L, str_midsec_receptivity_nega_score.pdf',width=4,height=5)
str_midsec$loci_phase<-factor(str_midsec$loci_phase,levels=(c('Normal_mid_secretory','Eutopic_mid_secretory')))
VlnPlot(str_midsec,'str_rep_nega1',group.by = 'loci_phase',cols = pal[4:5],pt.size = 0)+
  labs(title = 'Receptivity negative signature')+theme(legend.position="none")+
  theme(axis.text.x = element_text(angle=45))+
  theme(axis.text = element_text(color = 'black'))
dev.off()
