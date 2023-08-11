library(CMplot)
library(readr)



GWAS_result_4 <- read_delim("C:/Users/Admin/Desktop/2023_taishi_data_temp/temp/GWAS_ach_data/ach-gwas-12-22/ach-genotype_begin_s100_no_scaffold-R/ach-BL_BW_BH_merge.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)


nrow(GWAS_result_4)


CMplot(GWAS_result_4,type="p",plot.type="c",chr.labels=paste("Chr",c(1:24),sep=""),r=0.4,cir.legend=TRUE,outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg", memo="",dpi=300,file.output=T,verbose=TRUE,width=10,height=10)

CMplot(GWAS_result_4, plot.type=c('q'), multracks=T,threshold=c(1,0.05)/nrow(GWAS_result_4),threshold.lty=c(2,1), 
       threshold.lwd=c(1,1), threshold.col=c("blacK","red"), amplify=F,bin.size=1e6, 
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),signal.cex=c(0.5,0.5),
       file="jpg",memo="",dpi=300,file.output=T,verbose=TRUE,width=8,height=8)

CMplot(GWAS_result_4, plot.type=c('m'), multracks=T,threshold=c(1,0.05)/nrow(GWAS_result_4),threshold.lty=c(2,1), 
       threshold.lwd=c(1,1), threshold.col=c("blacK","red"), amplify=F,bin.size=1e6, 
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),signal.cex=c(0.5,0.5),
       file="jpg",memo="",dpi=300,file.output=T,verbose=TRUE,width=16,height=4)


CMplot(GWAS_result_4,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       file.output=T, verbose=TRUE)