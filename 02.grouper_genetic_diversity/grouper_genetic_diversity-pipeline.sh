# plink 质控ddrad-seq数据（斜带石斑鱼）

```{bash eval =F}
1.创建文件夹
mkdir ~/data/ach/EC-gwas-plink/2022_10_23_population_diversity
sed 's/.sort.bam//g' ach_TEST_var_flt_qual10_dp30000.vcf > ach_TEST_var_flt_qual10_dp30000_sample.vcf

vcftools --vcf ach_TEST_var_flt_qual10_dp30000_sample.vcf --plink --out ach_TEST_var_flt_qual10_dp30000_sample #有博主说，采用vcftools将vcf转化为plink格式可能会存在数据缺失
plink --file ach_TEST_var_flt_qual10_dp30000_sample --out c
plink --file ach_TEST_var_flt_qual10_dp30000_sample --missing
#采用plink转化vcf文件 （pink文件形式：map ped 属于文本形式；bed(二进制),bim,fam为一类）
plink --vcf ach_TEST_var_flt_qual10_dp30000.vcf --recode --out g --allow-extra-chr #正常plink文件
plink --vcf ach_TEST_var_flt_qual10_dp30000.vcf  --out g --allow-extra-chr --chr-set 50   #保存为二进制文件 默认操作23条染色体，需要修改，因为石斑鱼24条染色体
#将二进制转为正常plink文件形式 bed bim fam -----转化为----map ped 文件
plink --bfile g --recode --out test --allow-extra-chr --chr-set 50 
#查看 基因型个体和SNP数量
wc -l test.map test.ped
   1035043 test.map
       288 test.ped
   1035331 total
可以看出，共有 288 个基因型个体，共有 882048 个 SNP 数据
#对 SNP 缺失率进行筛选
plink --file test --missing --chr-set 50  --allow-extra-chr 

##先过滤 SNP在个体中缺失率高于 10%的 SNP
plink --bfile g --geno 0.10 --make-bed --out g_2 --allow-extra-chr --chr-set 50 

plink --bfile g_2 --recode --out test_2 --allow-extra-chr --chr-set 50 
wc -l test_2.map test_2.ped
   153692 test_2.map
      288 test_2.ped
   153980 total
可以看出，过滤了 1035043-153692 = 881351 个位点
#对样本缺失率进行筛选
##过滤缺失率高于10%的个体
plink --bfile g_2 --mind 0.10 --make-bed --out g_3 --allow-extra-chr --chr-set 50
plink --bfile g_3 --recode --out test_3 --allow-extra-chr --chr-set 50
wc -l test_3.map test_3.ped
   153692 test_3.map
      266 test_3.ped
   153958 total
可以看出，过滤了 288 - 266 = 12个体

#去掉 MAF 小于 0.05 的位点
plink --bfile g_3 --maf 0.05 --make-bed --out g_4 --allow-extra-chr --chr-set 50
plink --bfile g_4 --recode --out test_4 --allow-extra-chr --chr-set 50
wc -l test_4.map test_4.ped
   66901 test_4.map
     266 test_4.ped
   67167 totall
可以看出，过滤了 153692-66901 = 86791 个位点

#数据质控：哈温平衡
##计算所有位点的 HWE 的 P 值
plink --bfile g_4 --hardy --allow-extra-chr --chr-set 50
##提取哈温 p 值小于 0.0001 的位点
awk '{if($9 < 0.0001) print $0}' plink.hwe >plinkzoomhwe.hwe
##设定过滤标准 1e-4
plink --bfile g_4 --hwe 1e-4 --make-bed --out g_5 --allow-extra-chr --chr-set 50
plink --bfile g_5 --recode --out test_5 --allow-extra-chr --chr-set 50 
wc -l test_5.map test_5.ped
   61236 test_5.map
     266 test_5.ped
   61502 total
可以看出，过滤了 66901-61236 = 5665个位点


#Step6_转化成vcf
plink --bfile g_5 --allow-extra-chr --recode vcf-iid --out ach_TEST_var_flt_all_qual10_dp30000_plink_s6


#Step7_过滤indel/SNP
vcftools --vcf ach_TEST_var_flt_all_qual10_dp30000_plink_s6.vcf --remove-indels --recode --recode-INFO-all --out ach_TEST_var_flt_all_qual10_dp30000_plink_s7_SNP

vcftools --vcf ach_TEST_var_flt_all_qual10_dp30000_plink_s6.vcf --keep-only-indels --out ach_TEST_var_flt_all_qual10_dp30000_plink_s7_indel --recode --recode-INFO-all

#Step8_删除scaffold
sed 's/.sort.bam//g' ach_TEST_var_flt_all_qual10_dp30000_plink_s7_SNP.recode.vcf > ach_TEST_var_flt_all_qual10_dp30000_plink_s8_SNP.vcf

awk '{ if(/^#/ ||$1 >=1 && $1 <= 24) print $0}' ach_TEST_var_flt_all_qual10_dp30000_plink_s8_SNP.vcf >ach_TEST_var_flt_all_qual10_dp30000_plink_s8_2_SNP.vcf
#awk /^#/为正则匹配，将匹配文件中所有是#开头 https://blog.csdn.net/XiaodunLP/article/details/90262319


#Step9_给VCF文件加上ID
perl /home/xiajh2/data/ach/EC-gwas-plink/2022_10_23_population_diversity/VCF_add_id-master/VCF_add_id.pl ach_TEST_var_flt_all_qual10_dp30000_plink_s8_2_SNP.vcf ach_TEST_var_flt_all_qual10_dp30000_plink_s9_SNP.vcf

#step10_转换成bed文件
plink --vcf ach_TEST_var_flt_all_qual10_dp30000_plink_s9_SNP.vcf --make-bed --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP


plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --recode --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --allow-extra-chr --chr-set 50 
wc -l ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP.map ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP.ped


plink 转化vcf文件，最后生成bed fam bim文件，需要对fam文件进行修改，使用两表关联（R或excel）

#gemma 好像操作有点问题 待定使用 fam文件如果多余 也可能造成报错，比如本是266个个体，性状值选择了288个个体这样会导致出错！

../gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -gk 2 -o kin

../gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -k /home/xiajh2/data/ach/EC-gwas-plink/2022_10_23_population_diversity/output/kin.sXX.txt -lmm 1 -o yield




```



```{bash eval =F}
#admixture
 for i in {1..22};do admixture --cv ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP.bed $i |tee log${i}.out;done
 grep -h 'CV'  log*.out > CV.txt


#pca计算
##plink计算
plink --threads 16 --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --pca 10 --out pca --allow-extra-chr --chr-set 50

##GCTA计算PCA
gcta64 --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --make-grm  --out file_temp
gcta64 --grm file_temp --pca --out file_temp_pca

cut -d " " -f 2,3,4 file_temp_pca.eigenvec > ggplot_pca.txt

head -n 2 file_temp_pca.eigenval
14.5278
12.5329

#计算出SUM_TOTAL_EIGENVECS
awk '{sum += $1};END {print sum}' file_temp_pca.eigenval
272.498
PC1= (14.5278/272.498)*100%= 5.33%
PC2= (12.5329/272.498)*100%= 4.60%

table=read.table("ggplot_pca.txt",header = F,sep = " ")
library(ggplot2)
ggplot(table,aes(x=V2,y=V3))+geom_point()+
  xlab("PC1(5.33%)")+
  ylab("PC2(4.06%)")+
  theme_classic()

p2 <- ggplot(table,aes(x=V2,y=V3))+geom_point()+
  xlab("PC1(5.33%)")+
  ylab("PC2(4.06%)")+
  theme_classic()+theme(axis.title.x=element_text(face="bold", color="black",size=20,family="serif"),
              axis.title.y = element_text(colour="black",size=20,face="bold",family="serif"),
              axis.text.x = element_text(colour="black",size=16,face="bold",family="serif"),
              axis.text.y = element_text(colour="black",size=16,face="bold",family="serif"))+theme(axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(size = 16,face="bold",family="serif"),legend.title = element_text(size = 16,,face="bold"))


#计算亲缘关系

gcta64 --make-grm-gz --out root.gcta --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --autosome-num 24
####################################################R###########################################################################
setwd('C:/Users/Admin/Desktop/基于ddrad-seq对斜带石斑鱼放流群体的遗传多样性分析和基于简化基因组D-loop比对-pipeline/kinship/')

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

#计算遗传多样性指数
plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --hardy --allow-extra-chr --chr-set 50

###修改总群体的fam文件，原文件修改为后缀加bake--为了区分ABC
plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --recode  --allow-extra-chr --chr-set 50 --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_1

###修改总群体的fam文件，原文件修改为后缀加bake--全为ALL
plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --recode  --allow-extra-chr --chr-set 50 --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_ALL

##提取A群体计算
plink --file ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --keep population_A.txt --recode --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_A
plink --file ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_A --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_A
plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_A --hardy --allow-extra-chr --chr-set 50 --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_A
###将二进制文件种群改为A（修改fam文件--excel修改），然后再转为正常plink文件
plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_A --recode  --allow-extra-chr --chr-set 50 --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_A_1


##提取B群体计算
plink --file ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --keep population_B.txt --recode --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_B
plink --file ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_B --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_B

plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_B --hardy --allow-extra-chr --chr-set 50 --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_B
##提取C群体计算
plink --file ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --keep population_C.txt --recode --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_C
plink --file ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_C --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_C

plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_C --hardy --allow-extra-chr --chr-set 50 --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_population_C

#通过ROH计算近交系数
plink --allow-extra-chr --chr-set 50 --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --homozyg-window-snp 50 --homozyg-snp 50 --homozyg-kb 300 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --homozyg-window-het 03 --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP_ROH


```



