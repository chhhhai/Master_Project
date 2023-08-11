rm(list=ls())

library(clusterProfiler) # 富集分析R包
library(stringr) # 标签换行
library(AnnotationDbi)
library(DOSE)
library(ggplot2) # 绘图
library(ggrepel) # 标签相关
library(org.My.eg.db, lib.loc = "E:/EC/orf.protein.fa_eggnog/R_Library/")
setwd('C:/Users/Admin/Desktop/EC_GW_gwas_go_kegg/')
SNP_500kb_152genes <- read.csv("C:/Users/Admin/Desktop/EC_GW_gwas_go_kegg/SNP_500kb_152genes_2.csv")

gene_1 <- SNP_500kb_152genes$id
# GO富集分析
go <- enrichGO(gene = SNP_500kb_152genes$id,
               OrgDb = org.My.eg.db,
               keyType = 'GID',
               ont = 'ALL',
               qvalueCutoff = 0.6,
               pvalueCutoff = 0.6)
go.res <- data.frame(go) # 将GO结果转为数据框，方便后续分析（不转也可以，看个人习惯）
getwd()
write.csv(go.res,"Table_GO_result.csv",quote = F) # 输出GO富集分析结果
# 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:9,]
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:2,]
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:3,]
go.df <- rbind(goBP,goCC,goMF)
# 使画出的GO term的顺序与输入一致
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# 绘图
go_bar <- ggplot(data = go.df, # 绘图使用的数据
                 aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "Gene Number",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距




#KEGG 富集
emapper <- read_delim(
  file = 'E:/EC/orf.protein.fa_eggnog/out.emapper.annotations',
  "\t", escape_double = FALSE, col_names = FALSE,
  comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                Gene_Symbol = X9, 
                GO = X10, 
                KO = X12, 
                Pathway = X13, 
                OG = X7, 
                Gene_Name = X8)

pathway2gene <- dplyr::select(emapper, Pathway, GID) %>%
  separate_rows(Pathway, sep = ',', convert = F) %>%
  filter(str_detect(Pathway, 'ko')) %>%
  mutate(Pathway = str_remove(Pathway, 'ko'))
library(magrittr)
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
  colnames(keggpathid2name.df) <- c("path_id","path_name")
  return(keggpathid2name.df)
}
pathway2name <- get_path2name()
de_ekp <- enricher(gene_1,
                   TERM2GENE = pathway2gene,
                   TERM2NAME = pathway2name,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
kegg.res <- data.frame(de_ekp) 
write.csv(kegg.res,"Table_KEGG_result.csv",quote = F) # 输出GO富集分析结果
names(kegg.res)[names(kegg.res) == 'ID'] <- 'Pathways_ID'

# 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
kegg <- subset(kegg.res,)[1:2,]
kegg.df <- rbind(kegg)
# 使画出的GO term的顺序与输入一致
kegg.df$Description <- factor(kegg.df$Description,levels = rev(kegg.df$Description))
# 绘图
kegg_bar <- ggplot(data = kegg.df, # 绘图使用的数据
                   aes(x = Description, y = Count,fill = Pathways_ID))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = " KEGG pathways",y = "Gene Number",title = "Barplot of Enriched  KEGG pathways")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距

mytheme =theme(axis.title.x=element_text(face="bold", color="black",size=20,family="serif"),
               axis.title.y = element_text(colour="black",size=20,face="bold",family="serif"),
               axis.text.x = element_text(colour="black",size=16,face="bold",family="serif"),
               axis.text.y = element_text(colour="black",size=16,face="bold",family="serif"))

library(patchwork)

p1 <- go_bar+mytheme

p1
p2 <- kegg_bar+mytheme

p1/p2

ggsave(p1/p2,filename = "GO_Barplot.tiff",width = 16,height = 14,dpi=300)