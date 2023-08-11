library(ggplot2)
library(patchwork)
library(ggpubr)

#1. 绘制初级版箱线图

#绘图，以HT（单倍型）为x轴，表型AUC为y轴绘制初级版箱线图
gwas_BW <- read_csv("C:/Users/Admin/Desktop/gwas-xiiangxiantu/gwas_BW.csv")
ggboxplot(gwas_BW,x = "type",y = "data")

my_comparisons=list(c("smallest group", "biggest group"))

#2. 绘制改良版箱线图

#绘制改良版箱线图

#add = "jitter"-添加散点

#stat_compare_means(comparisons = my_comparisons-引入不同HT差异显著性T检验的P-value，


p1 <- ggboxplot(gwas_BW,x = "type",y = "data",color ='type', add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,method = "t.test",label = "p.signif")+
  
  labs(x = 'Groups', y = "Body Weight (g)")+theme(axis.title.x=element_text(face="bold", color="black",size=20,family="serif"),
                                                  axis.title.y = element_text(colour="black",size=20,face="bold",family="serif"),
                                                  axis.text.x = element_text(colour="black",size=16,face="bold",family="serif"),
                                                  axis.text.y = element_text(colour="black",size=16,face="bold",family="serif"))
p1

gwas_BL <- read_csv("C:/Users/Admin/Desktop/gwas-xiiangxiantu/gwas_BL.csv")
my_comparisons=list(c("smallest group", "biggest group"))

p2 <- ggboxplot(gwas_BL,x = "type",y = "data",color ='type', add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,method = "t.test",label = "p.signif")+
  
  labs(x = 'Groups', y = "Body Length (cm)")+theme(axis.title.x=element_text(face="bold", color="black",size=20,family="serif"),
                                                   axis.title.y = element_text(colour="black",size=20,face="bold",family="serif"),
                                                   axis.text.x = element_text(colour="black",size=16,face="bold",family="serif"),
                                                   axis.text.y = element_text(colour="black",size=16,face="bold",family="serif"))

p2
gwas_BH <- read_csv("C:/Users/Admin/Desktop/gwas-xiiangxiantu/gwas-BH.csv")
my_comparisons=list(c("smallest group", "biggest group"))

p3 <- ggboxplot(gwas_BH,x = "type",y = "data",color ='type', add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,method = "t.test",label = "p.signif")+
  
  labs(x = 'Groups', y = "Body Height (cm)")+theme(axis.title.x=element_text(face="bold", color="black",size=20,family="serif"),
                                                   axis.title.y = element_text(colour="black",size=20,face="bold",family="serif"),
                                                   axis.text.x = element_text(colour="black",size=16,face="bold",family="serif"),
                                                   axis.text.y = element_text(colour="black",size=16,face="bold",family="serif"))
P3
p3

p <- p1+p2+p3

p

ggsave('test.tiff', p, device = "tiff", dpi = 300, 
       width=16, height=8, unit = "in")
getwd()
