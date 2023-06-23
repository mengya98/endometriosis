setwd('G:/endometriosis/Analysis/EPI')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(ggsci)

epi<-readRDS('./epi.rds')
meta<-data.frame(epi@meta.data)

#S4A Volcano plot based on P15 P23 P36#
P15_P23_EC<-rownames(meta[meta$patient%in%c('P15','P23') & meta$origin%in%c('Ectopic'),])
P15_P23_EU<-rownames(meta[meta$patient%in%c('P15','P23') & meta$origin%in%c('Eutopic'),])
P36_EC<-rownames(meta[meta$patient=='P36' & meta$origin%in%c('Ectopic'),])
P36_EU<-rownames(meta[meta$patient=='P36' & meta$origin%in%c('Eutopic'),])
#P15_P23#
P15_P23_ECvsEU<-FindMarkers(epi,ident.1 = P15_P23_EC, ident.2 = P15_P23_EU,only.pos = F,min.pct = 0.25, logfc.threshold = 0)
P15_P23_ECvsEU$gene<-rownames(P15_P23_ECvsEU)
P36_ECvsEU<-FindMarkers(epi,ident.1 = P36_EC, ident.2 = P36_EU,only.pos = F,min.pct = 0.25, logfc.threshold = 0)
P36_ECvsEU$gene<-rownames(P36_ECvsEU)
P15P23_P36_ECvsEU<-merge(P15_P23_ECvsEU,P36_ECvsEU,by = 'gene')
P15P23_P36_ECvsEU<-P15P23_P36_ECvsEU[,c(1,3,8)]
colnames(P15P23_P36_ECvsEU)<-c('gene','Proliferative','Secretory')

label_genes<-unique(c('LAMB3','LAMA4','SPARC','NRP1','NNMT','RDH10','STRA6',
                      'MGP','RARRES1','RUNX3','CYP1B1','RBP1','PGR','FHL2'))
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

pdf("F4A, Volcano_P15P23_P36_ECvsEU.pdf",width=4.0,height=4)
ggplot(P15P23_P36_ECvsEU,aes(x=Proliferative,y=Secretory,color=label_genes_sig,size=label_genes_sig))+geom_point()+
  theme_classic()+
  scale_color_manual(name='',values =c("label"="#91331FFF","con"="grey"))+
  scale_size_manual(values = c(0.5,1))+
  labs(x = 'EC vs EU in proliferative phase (log2)', y = 'EC vs EU in secretory phase (log2)', title = "")+
  geom_text_repel(data = P15P23_P36_ECvsEU, aes(x = Proliferative, 
                                                y = Secretory, 
                                                label = label_genes),
                  size = 2,box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.3, "lines"), 
                  max.overlaps = 2000,
                  segment.color = "black", 
                  color="black",
                  show.legend = FALSE)+
  theme(legend.position = 'none')
dev.off() 
