setwd('./kinship/')
a <- read.table('../PCA/ggplot_pca_new.txt',header = T)
b <- read.csv('kinship_cluster_k3.csv')
library('tidyverse')
c <-left_join(a,b,by=c('id'='sample')) 

install.packages("do")
library(do)
c$Cluster <- Replace(c$Cluster,"1", "A") 
c$Cluster <- Replace(c$Cluster,"2", "B") 
c$Cluster <- Replace(c$Cluster,"3", "C") 


ggplot(data = c,aes(x=PC1,y=PC2,color=Cluster)

p2 <- ggplot(c,aes(x=PC1,y=PC2,color=factor(Cluster)))+geom_point()+ 
  # stat_ellipse(level = 0.95, size = 1) + 
  #stat_ellipse(aes(fill=Cluster),
               #type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  xlab("PC1(5.33%)")+
  ylab("PC2(4.06%)")+
  theme_classic()+theme(axis.title.x=element_text(face="bold", color="black",size=20,family="serif"),
                        axis.title.y = element_text(colour="black",size=20,face="bold",family="serif"),
                        axis.text.x = element_text(colour="black",size=16,face="bold",family="serif"),
                        axis.text.y = element_text(colour="black",size=16,face="bold",family="serif"))+theme(axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(size = 16,face="bold",family="serif"),legend.title = element_text(size = 16,,face="bold"))

p2+guides(color=guide_legend(title = "population"))

