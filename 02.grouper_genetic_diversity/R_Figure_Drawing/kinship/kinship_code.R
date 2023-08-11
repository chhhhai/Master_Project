setwd('C:/Users/Admin/Desktop/kinship/')

library(reshape2)
tmp <- read.table(gzfile("root.gcta.grm.gz"), header = F, stringsAsFactors = F)
ids <- read.table("root.gcta.grm.id", header = F, stringsAsFactors = F)
tmp <- tmp[,c(1,2,4)]
result_matrix <- acast(tmp, V1~V2, value.var="V4", drop = F)
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
result_full <- makeSymm(result_matrix)
diag(result_full) <- 2
result_df <- as.data.frame(result_full)
row.names(result_df) <- ids$V2
colnames(result_df) <- ids$V2
write.table(result_df, file = "gcta.kinship.txt", row.names = T, col.names = NA, sep = "\t", quote = F)

rm(list = ls())
gcta<- read.table("gcta.kinship.txt")
library("pheatmap")
#pheatmap(gcta, fontsize_row = 0.3, fontsize_col = 0.3,filename = "test.pdf")保存热图文件
p=pheatmap(gcta, fontsize_row = 0.3, fontsize_col = 0.3)
row_cluster <- cutree(p$tree_row,k=3)#k=3 分为三个亚群
newOrder <- gcta[p$tree_row$order,]
newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder), names(row_cluster))]
colnames(newOrder)[ncol(newOrder)]="Cluster"
write.csv(newOrder, "FPKM_cluster.csv")