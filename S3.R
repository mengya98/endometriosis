setwd('G:/endometriosis/Analysis/p40/EPI/remove_cilia/all/')

library(ggplot2) 
library(readxl)

#S3A#
P23_GO<-read_xlsx('./DEG/EC/GEO/P23_ec_up_meta.xlsx',sheet = 3)
P23_GO$Description<-factor(P23_GO$Description,levels = rev(P23_GO$Description))
pdf("S3A, P23_ECvsEU_upregulation_GOterms.pdf",width = 6.7,height = 4)
ggplot(P23_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#D2AF81FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-LogP", title = "Go terms of EC vs EU in P23")
dev.off()

#S3B#
P36_GO<-read_xlsx('./DEG/EC/GEO/P36_ec_up_meta.xlsx',sheet = 3)
P36_GO$Description<-factor(P36_GO$Description,levels = rev(P36_GO$Description))
pdf("S3B, P36_ECvsEU_upregulation_GOterms.pdf",width = 7.5,height = 4)
ggplot(P36_GO,aes(x=Description,y=-LogP,fill=-LogP))+
  coord_flip()+theme_classic()+
  geom_bar(stat="identity",width=0.8,fill="#D2AF81FF")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(plot.title=element_text(size=20,hjust =0.5))+
  labs(x = "",y = "-LogP", title = "Go terms of EC vs EU in P36")
dev.off()