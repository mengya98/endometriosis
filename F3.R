setwd('G:/endometriosis/Analysis/p40/EU')

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(ggsci)
library(readxl)

#F3A GO term in proliferative phase#
pro_term<-data.frame(term=c('Cilium organization','positive regulation of cell growth','WNT SIGNALING','regulation of apoptotic signaling pathway'),
                     LogP=c(14.8559300379,3.452741703,2.39344536,2.1239648141))
pro_term$term<-factor(pro_term$term,levels = rev(pro_term$term))
pdf("F3A, EPI_EU_PvsN_pro_GOterms_upregulation.pdf",width = 4.6,height = 3)
ggplot(pro_term,aes(x=term,y=LogP,fill=LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.6,fill="#91331FFF")+
  theme(axis.text = element_text(size = 8, color = "black"))+
  # theme(axis.title.x=element_text(size=16))+
  # theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=16,hjust =1.2))+
  labs(x = "",y = "-LogP", title = "Go terms of EU_pro vs N_pro in EPI")
dev.off()

#F3B#
epi_EU<-readRDS('./EPI_EU.rds')
pseu1<-read.table('EPI_eu_pseudotime.txt',header = T)
epi_data<-data.frame(epi_EU@assays$RNA@data)
epi_data<-data.frame(t(epi_data),check.names = F)
epi_wnt_cellcycle<-c('LRP5','FZD4','CCNA1','TOP2A')
epi_data_wnt_cellcycle<-epi_data[,colnames(epi_data)%in%(epi_wnt_cellcycle)]
epi_data_wnt_cellcycle$Identity<-rownames(epi_data_wnt_cellcycle)
pseu_wnt_cellcycle<-merge(pseu1,epi_data_wnt_cellcycle,by = 'Identity')
p1<-ggplot(pseu_wnt_cellcycle)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = epi_wnt_cellcycle[1], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt_cellcycle,size=0.1,aes_string(x = 'Pseudotime', y = epi_wnt_cellcycle[1], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_wnt_cellcycle[1])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p2<-ggplot(pseu_wnt_cellcycle)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_wnt_cellcycle[2], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt_cellcycle,size=0.1,aes_string(x = 'Pseudotime', y = epi_wnt_cellcycle[2], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_wnt_cellcycle[2])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p3<-ggplot(pseu_wnt_cellcycle)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_wnt_cellcycle[3], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt_cellcycle,size=0.1,aes_string(x = 'Pseudotime', y = epi_wnt_cellcycle[3], color= 'origin',alpha=0.001))+
  labs(y='Expression',title = epi_wnt_cellcycle[3])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
p4<-ggplot(pseu_wnt_cellcycle)+
  theme_classic()+
  geom_smooth(aes_string(x = 'Pseudotime', y = epi_wnt_cellcycle[4], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt_cellcycle,size=0.1,aes_string(x = 'Pseudotime', y = epi_wnt_cellcycle[4], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = epi_wnt_cellcycle[4])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(4,6.9,10),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position="none")
pdf("F3B, pseu_EU_PvsN_wnt_cellcycle.pdf", width = 9, height = 6)
patchwork::wrap_plots(p1,p3,p2,p4,ncol = 2)
dev.off()

#F3C#
selected_rows<-read.table('rows_use.txt',sep='\t')
selected_rows=selected_rows$V1
selected_columns<-read.table('col_use.txt',sep='\t')
selected_columns=selected_columns$V1

all_pval<-read.table('./pvalues.txt', header=T, stringsAsFactors=F, sep='\t', comment.char='', check.names=F)
all_means<-read.table('./means.txt', header=T, stringsAsFactors=F, sep='\t', comment.char='', check.names=F)

intr_pairs<-all_pval$interacting_pair
all_pval<-all_pval[,-c(1:11)]
all_means<-all_means[,-c(1:11)]

if(is.null(selected_rows)){
  selected_rows=intr_pairs
}

if(is.null(selected_columns)){
  selected_columns=colnames(all_pval)
}

sel_pval<-all_pval[match(selected_rows, intr_pairs), selected_columns]
sel_means<-all_means[match(selected_rows, intr_pairs), selected_columns]

df_names<-expand.grid(selected_rows, selected_columns)
pval<-unlist(sel_pval)
pval[pval==0] = 0.0009
plot.data<-cbind(df_names,pval)
pr<-unlist(as.data.frame(sel_means))
pr[pr==0] = 1
plot.data<-cbind(plot.data,log2(pr))
colnames(plot.data)<-c('pair', 'clusters', 'pvalue', 'mean')

r_n=3:6
rm_n=c(48+r_n,54+r_n,60+r_n,66+r_n,72+r_n,78+r_n,84+r_n,90+r_n)
l_n=1:2
lm_n=c(l_n,6+l_n,12+l_n,18+l_n,24+l_n,30+l_n,36+l_n,42+l_n)
plot.data_x=plot.data[c(rm_n,lm_n),]
plot.data_x_r=plot.data_x[1:32,]
plot.data_x_r$clusters=c(rep('epi_Eutopic_early_proliferative|str_Eutopic_early_proliferative',4),rep('epi_Eutopic_late_proliferative|str_Eutopic_late_proliferative',4),
                         rep('epi_Eutopic_early_secretory|str_Eutopic_early_secretory',4),rep('epi_Eutopic_mid_secretory|str_Eutopic_mid_secretory',4),
                         rep('epi_Normal_early_proliferative|str_Normal_early_proliferative',4),rep('epi_Normal_late_proliferative|str_Normal_late_proliferative',4),
                         rep('epi_Normal_early_secretory|str_Normal_early_secretory',4),rep('epi_Normal_mid_secretory|str_Normal_mid_secretory',4))

plot.data_x_r$pair=rep(c('FZD3_WNT5A','FZD4_WNT2','FGFR2_FGF7','VEGFB_NRP1'),8)
plot.data=rbind(plot.data_x_r,plot.data_x[33:48,])
plot.data$clusters=factor(plot.data$clusters,levels<-c('epi_Eutopic_early_proliferative|str_Eutopic_early_proliferative','epi_Eutopic_late_proliferative|str_Eutopic_late_proliferative',
                                                        'epi_Eutopic_early_secretory|str_Eutopic_early_secretory','epi_Eutopic_mid_secretory|str_Eutopic_mid_secretory',
                                                        'epi_Normal_early_proliferative|str_Normal_early_proliferative','epi_Normal_late_proliferative|str_Normal_late_proliferative',
                                                        'epi_Normal_early_secretory|str_Normal_early_secretory','epi_Normal_mid_secretory|str_Normal_mid_secretory'))
plot.data$pair=factor(plot.data$pair,levels = unique(plot.data$pair))

plot.data$clusters=gsub("proliferative","Pro",plot.data$clusters)
plot.data$clusters=gsub("secretory","Sec",plot.data$clusters)
unique(plot.data$clusters)
plot.data$clusters[plot.data$clusters=="epi_Eutopic_early_Pro|str_Eutopic_early_Pro"]="epi_str_Eutopic_early_Pro"
plot.data$clusters[plot.data$clusters=="epi_Eutopic_late_Pro|str_Eutopic_late_Pro"]="epi_str_Eutopic_late_Pro"
plot.data$clusters[plot.data$clusters=="epi_Eutopic_early_Sec|str_Eutopic_early_Sec"]="epi_str_Eutopic_early_Sec"
plot.data$clusters[plot.data$clusters=="epi_Eutopic_mid_Sec|str_Eutopic_mid_Sec"]="epi_str_Eutopic_mid_Sec"
plot.data$clusters[plot.data$clusters=="epi_Normal_early_Pro|str_Normal_early_Pro"]="epi_str_Normal_early_Pro"
plot.data$clusters[plot.data$clusters=="epi_Normal_late_Pro|str_Normal_late_Pro"]="epi_str_Normal_late_Pro"
plot.data$clusters[plot.data$clusters=="epi_Normal_early_Sec|str_Normal_early_Sec"]="epi_str_Normal_early_Sec"
plot.data$clusters[plot.data$clusters=="epi_Normal_mid_Sec|str_Normal_mid_Sec"]="epi_str_Normal_mid_Sec"
plot.data$clusters=factor(plot.data$clusters,levels = c("epi_str_Eutopic_early_Pro","epi_str_Eutopic_late_Pro","epi_str_Eutopic_early_Sec","epi_str_Eutopic_mid_Sec",
  "epi_str_Normal_early_Pro","epi_str_Normal_late_Pro","epi_str_Normal_early_Sec","epi_str_Normal_mid_Sec"))
my_palette <- colorRampPalette(rev(c("#9E0142","#D53E4F","#F46D43",'orange',"#FDAE61","#FEE08B","#E6F598","#BABABA","#878787","#4D4D4D","#1A1A1A")), alpha=TRUE)(n=399)

pdf('F3C, epi_str_interaction.pdf',height = 4.3,width = 7.2)  
ggplot(plot.data,aes(y=pair,x=clusters)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=16, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.8, linetype = "solid", colour = "black"))
dev.off()

#F3D#
epi_cilia<-c('TPPP3','SNTN')
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
pdf("F3D, pseu_epi_EU_PvsN_cilia.pdf", width = 5, height = 6)
patchwork::wrap_plots(p1,p2,ncol = 1)
dev.off()

#F3G#
str<-readRDS('./STR_EU.rds')
pseu2<-read.table('STR_eu_pseudotime.txt',header = T)
str_data<-data.frame(str@assays$RNA@data)
str_data<-data.frame(t(str_data),check.names = F)
str_wnt<-c('SFRP4','ID2','ID3','TP53')
str_data_wnt<-str_data[,colnames(str_data)%in%(str_wnt)]
str_data_wnt$Identity<-rownames(str_data_wnt)
pseu_wnt<-merge(pseu2,str_data_wnt,by = 'Identity')
p1<-ggplot(pseu_wnt)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = str_wnt[1], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt,size=0.1,aes_string(x = 'Pseudotime', y = str_wnt[1], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = str_wnt[1])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(10,14.17,14.82),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.text = element_text(color = 'black'),axis.title = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))+
  theme(legend.position="none")
p2<-ggplot(pseu_wnt)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = str_wnt[2], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt,size=0.1,aes_string(x = 'Pseudotime', y = str_wnt[2], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = str_wnt[2])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(10,14.17,14.82),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.text = element_text(color = 'black'),axis.title = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))+
  theme(legend.position="none")
p3<-ggplot(pseu_wnt)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = str_wnt[3], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt,size=0.1,aes_string(x = 'Pseudotime', y = str_wnt[3], color= 'origin',alpha=0.001))+
  labs(y='Expression',title = str_wnt[3])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(10,14.17,14.82),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.text = element_text(color = 'black'),axis.title = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))+
  theme(legend.position="none")
p4<-ggplot(pseu_wnt)+
  theme_classic()+
  geom_smooth(method = 'loess',aes_string(x = 'Pseudotime', y = str_wnt[4], color= 'origin'),se = FALSE)+
  geom_point(data=pseu_wnt,size=0.1,aes_string(x = 'Pseudotime', y = str_wnt[4], color= 'origin',alpha=0.001))+
  labs(x='',y='Expression',title = str_wnt[4])+
  scale_color_manual(values = c('Normal'='#197EC0FF','Eutopic'='#C80813FF'))+
  geom_vline(xintercept = c(10,14.17,14.82),size=0.1)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.text = element_text(color = 'black'),axis.title = element_text(color = 'black'),axis.ticks = element_line(color = 'black'))+
  theme(legend.position="none")
pdf("F3G, pseu_str_EU_PvsN_wnt.pdf", width = 9, height = 6)
patchwork::wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()

#F3H#
str_TF<-c('HSPA1A','HSPA1B','HSPA6')
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
pdf("F3H, pseu_str_EU_PvsN_TF.pdf", width = 14, height = 3)
patchwork::wrap_plots(p1,p2,p3,ncol = 3)
dev.off()


