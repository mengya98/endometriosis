setwd('G:/endometriosis/Analysis/EPI')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(ggsci)

epi<-readRDS('./epi.rds')

#S2A epi QC#
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

pdf('S2A, epi_QCplot.pdf',width=16,height=6)
t2+t1+t3+t4+t5+t6+t7+plot_layout(widths = c(1,1,1,3,1,1,1),ncol = 7)
dev.off()

#S2B NNMT+ cells proportions in N/EU/EC#
NNMT_exp<-meta2
NNMT_exp$exp<-data.frame(epi@assays$RNA@data['NNMT',])
NNMT_exp<-NNMT_exp[,c(8,16)]
colnames(NNMT_exp)[2]<-'expression'
NNMT_exp$NNMT<-'NNMT+'
NNMT_exp$NNMT[NNMT_exp$expression<=3]<-'NNMT-'
pdf("S2B, NNMT+.exp_3.origin.proportion.pdf",width = 4, height = 3)
ggplot(NNMT_exp,aes(x=origin,fill=NNMT))+
  geom_bar(stat = 'count',
           position = 'fill',
           width = 0.6)+
  labs(x = "",y = "Proportion", title = "",fill='NNMT expression')+
  theme_classic()+
  scale_fill_manual(values = c('#9FDAA4','#468141'))+
  #guides(fill=guide_legend(title=),size)+
  theme(axis.text = element_text(size = 10, color = "black"),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        legend.title = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        legend.text = element_text(color="black", size = 10))+
  theme(axis.text = element_text(colour = 'black'))
dev.off()

#S2C#
#UMAP seperated for N EU
epi_N<-subset(epi,origin%in%c('Normal'))
pdf('S2C, epi_N_celltype.pdf',width=6.6,height=5)
DimPlot(epi_N, reduction = "umap", group.by = 'celltype',pt.size = 1,label = F,
        cols = c('#7876B1FF','#CD9B1D',pal_simpsons()(5)[c(2,5)]))+
  theme_bw()+
  labs(x = "UMAP1", y = "UMAP2",title = 'Normal')+
  theme(axis.text= element_blank(), axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5))+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))
dev.off()
remove(epi_N)

epi_EU<-subset(epi,origin%in%c('Eutopic'))
pdf('S2C, epi_EU_celltype.pdf',width=6.6,height=5)
DimPlot(epi_EU, reduction = "umap", group.by = 'celltype',pt.size = 1,label = F,
        cols = c('#7876B1FF','#CD9B1D',pal_simpsons()(5)[c(2,5)]))+
  theme_bw()+
  labs(x = "UMAP1", y = "UMAP2",title = 'Eutopic')+
  theme(axis.text= element_blank(), axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5))+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))
dev.off()
remove(epi_EU)

#S2D epi Menstrual cycle#
pdf('S2D, epi_menstrual_cycle.pdf',width=5.0,height=4)
DimPlot(epi, reduction = "umap", pt.size = 1,label = F,group.by = 'bigphase',cols = c("#B09C85FF","#E64B35FF"))+
  labs(x = "UMAP1", y = "UMAP2",title = '')+
  theme_bw()+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text= element_blank(), axis.ticks=element_blank())
dev.off()

#S2E#
sec_EU_N<-rownames(meta2[meta2$bigphase%in%c('Secretory') & meta2$control=='N' & meta2$seurat_clusters%in%c(0,1,2,4,5,6,7,8),])
sec_EU_P<-rownames(meta2[meta2$bigphase%in%c('Secretory') & meta2$control=='P' & meta2$seurat_clusters%in%c(0,1,2,4,5,6,7,8),])
volcano_sec<-FindMarkers(epi,ident.1 = sec_EU_P, ident.2 = sec_EU_N,only.pos = F,min.pct = 0.25, logfc.threshold = 0)
volcano_sec<-volcano_sec[volcano_sec$p_val_adj<1,]
volcano_sec$gene<-rownames(volcano_sec)
volcano_sec<-volcano_sec[order(volcano_sec$avg_log2FC,decreasing = T),]

volcano_sec$label_genes<-''
volcano_sec$label_genes_sig<-''
volcano_sec$label_genes_sig<-'con'
volcano_sec$label_genes_sig[volcano_sec$p_val_adj <= 0.05 & (volcano_sec$avg_log2FC >= 0.25 | volcano_sec$avg_log2FC <= -0.25)]<-'mark'
volcano_sec$label_genes_sig[volcano_sec$p_val_adj <= 1e-10 & (volcano_sec$avg_log2FC >= 4 | volcano_sec$avg_log2FC <= -3)]<-'label'
volcano_sec$label_genes_sig[volcano_sec$p_val_adj <= 1e-18]<-'label'
for (i in 1:nrow(volcano_sec)) {
  if (volcano_sec[i,8]=='label') {
    volcano_sec[i,7]<-volcano_sec[i,6]
  }
}
pdf(file = "S2E, volcano_sec_PvsN_DEGs.pdf",width=4,height=4)
ggplot(volcano_sec,aes(x=avg_log2FC,y=-log10(p_val_adj),color=label_genes_sig,size=label_genes_sig))+geom_point()+
  theme_classic()+
  labs(x = 'log2(FoldChange)', y = '-log10(FDR)', title = "Differentially expressed genes in secretory phase")+
  geom_text_repel(data = volcano_sec, aes(x = avg_log2FC, 
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

#S2F#
epi_NEU<-subset(epi,origin%in%c('Normal','Eutopic'))
meta<-data.frame(epi_NEU@meta.data)
load('./monocle.Rdata')
component<-data.frame(t(cds_1@reducedDimS))
component$celltype<-meta$celltype
target<-rownames(component[component$X1 > -10 & component$X1 < -1.5,])
meta$target<-'Non-target'
meta$identity<-rownames(meta)
meta$target[meta$identity%in%target]<-'Target'
meta_pro<-meta[meta$celltype=='Glandular proliferative',]
pdf("S2F, EPI.target.origin.pro.proportion.pdf",width = 4.2, height = 4)
ggplot(meta_pro,aes(x=origin,fill=target))+
  geom_bar(stat = 'count',
           position = 'fill',
           width = 0.6)+ 
  labs(x = "",y = "Proportion", title = "",fill='Epi subsets')+
  theme_classic()+
  scale_fill_manual(values = c('#9FDAA4','#468141'))+
  #guides(fill=guide_legend(title=),size)+
  theme(axis.text = element_text(size = 10, color = "black"),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        legend.title = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        legend.text = element_text(color="black", size = 10))+
  theme(axis.text = element_text(colour = 'black'))
dev.off()
remove(epi_NEU,meta)

#S2G#
meta2_EU_pro<-meta2_EU[meta2_EU$bigphase=='Proliferative',]
pdf("S2G, EPI.cilia.origin.pro.proportion.pdf",width = 3.5, height = 3)
ggplot(meta2_EU_pro,aes(x=origin,fill=cilia))+
  geom_bar(stat = 'count',
           position = 'fill',
           width = 0.6)+
  labs(x = "",y = "Proportion", title = "",fill='Epi subsets')+
  theme_classic()+
  scale_fill_manual(values = c('#9FDAA4','#468141'))+
  #guides(fill=guide_legend(title=),size)+
  theme(axis.text = element_text(size = 10, color = "black"),
        #axis.text.x = element_text(angle = 45,hjust = 1),
        legend.title = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        legend.text = element_text(color="black", size = 10))+
  theme(axis.text = element_text(colour = 'black'))
dev.off()
remove(meta2_EU_pro)
