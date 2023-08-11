rm(list = ls())
setwd("C:/Users/Admin/Desktop/12s_R/")
getwd()
library(MicrobiotaProcess)
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
otu <- "otu_table.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "group.txt"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeqza=tree)

p_rare <- ggrarecurve(obj=ps_dada2, 
                      factorNames="group",
                      indexNames=c("Observe"), 
                      chunks=300) +
  theme(legend.spacing.y=unit(0.02,"cm"),
        legend.text=element_text(size=4))+
  theme_bw()

p_rare

tiff(filename = "p_rare.tif",
     width = 3200, height = 3200, units = "px", pointsize = 12,
     compression = "lzw", res = 400)
p_rare
dev.off()



alphaobj <- get_alphaindex(ps_dada2)

head(as.data.frame(alphaobj))

p1 <- ggbox(alphaobj, geom="boxplot", factorNames="group",indexNames=c("Chao1"),compare = F) + 
  scale_fill_manual(values=c("#2874C5", "#EABF00","red","green","#00BFFF","yellow","pink","violet"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
p1

p2 <- ggbox(alphaobj, geom="boxplot", factorNames="group",indexNames=c("Shannon"),compare = F) + 
  scale_fill_manual(values=c("#2874C5", "#EABF00","red","green","#00BFFF","yellow","pink","violet"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
p2

p3 <- ggbox(alphaobj, geom="boxplot", factorNames="group",indexNames=c("Simpson"),compare = F) + 
  scale_fill_manual(values=c("#2874C5", "#EABF00","red","green","#00BFFF","yellow","pink","violet"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
p3

p4 <- ggbox(alphaobj, geom="boxplot", factorNames="group",indexNames=c("ACE"),compare = F) + 
  scale_fill_manual(values=c("#2874C5", "#EABF00","red","green","#00BFFF","yellow","pink","violet"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
p4

(p1 | p2) / (p3 | p4)

p <- (p1 | p2) / (p3 | p4)

tiff(filename = "alphaobj.tif",
     width = 3200, height = 3200, units = "px", pointsize = 12,
     compression = "lzw", res = 400)
p
dev.off()

genustax <- get_taxadf(obj=ps_dada2, taxlevel=6)
genusbar <- ggbartax(obj=genustax, facetNames="group", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
genusbar
pdf("genusbar.pdf",width=12,height=8)
genusbar
dev.off()


genustax <- get_taxadf(obj=ps_dada2, taxlevel=7)
genusbar_1 <- ggbartax(obj=genustax, facetNames="group", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
genusbar_1
pdf("genusbar_1.pdf",width=12,height=8)
genusbar_1
dev.off()




pcares <- get_pca(obj=ps_dada2, method="hellinger")
pcaplot <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
                      pc=c(1,2),factorNames=c("group"), ellipse=TRUE) + 
  scale_color_manual(values=c("#2874C5", "#EABF00","red","green","#00BFFF","yellow","pink","violet"))
pcaplot
tiff(filename = "pca.tif",
     width = 3200, height = 3200, units = "px", pointsize = 12,
     compression = "lzw", res = 400)
pcaplot
dev.off()

pcoares <- get_pcoa(obj=ps_dada2, 
                    distmethod="euclidean", method="hellinger")
pcoaplot <- ggordpoint(obj=pcoares, biplot=TRUE,
                       speciesannot=TRUE,pc = c(2,1),
                       factorNames=c("group"), ellipse=T)
pcoaplot
tiff(filename = "pcoa.tif",
     width = 3200, height = 3200, units = "px", pointsize = 12,
     compression = "lzw", res = 400)
pcoaplot
dev.off()