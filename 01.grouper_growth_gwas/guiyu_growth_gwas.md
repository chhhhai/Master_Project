# LDR_guiyu_growth_gwas

```{bash eval =F}
#1.拆分数据
for ((i=1; i<=8; i++))
do
cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/${i}
mkdir ${i}_radtags
/software/stacks-2.4/process_radtags -1 ${i}_1.clean.fq.gz -2 ${i}_2.clean.fq.gz -o ./${i}_radtags/ -b /home/xjh3/data/ach/guiyu_gwas/LDR/AD_seq_barcode.txt -r -c -q --renz_1 ecoRI --renz_2 mspI -i gzfastq > ${i}.radtags.log 2>&1
done



#2.mapping(参考基因组：GCF_020085105.1)
##2.1建立索引
bowtie2-build /home/xjh3/data/ach/guiyu_gwas/ncbi_ref/GCF_020085105.1_ASM2008510v1_genomic.fna /home/xjh3/data/ach/guiyu_gwas/ncbi_ref/gw.fa_index


for ((j=1; j<=8; j++))
do
echo cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/${j}
	for i in *.rem.1.fq.gz
	do
	i=${i%.rem.1.fq.gz*}
	echo bowtie2 -x /home/xjh3/data/ach-EC-rna-seq/ref/Eco_Chr.fa_index -1 ${i}.1.fq.gz -2 ${i}.2.fq.gz -S ${i}.sam -p 12
	done	
done


for ((i=1; i<=8; i++))
do
echo cd ${i}
	for i in *.rem.1.fq.gz
	do
	i=${i%.rem.1.fq.gz*}
	echo bowtie2 -x /home/xjh3/data/ach-EC-rna-seq/ref/Eco_Chr.fa_index -1 ${i}.1.fq.gz -2 ${i}.2.fq.gz -S ${i}.sam -p 1 >${i}.log 2>&1
	done	
done


```



# 使用梁旭方课题组的参考基因组进行gwas分析

```{bash eval =F}
#2.mapping(liangxufang_ref)
bowtie2-build /home/xjh3/data/ach/guiyu_gwas/hzau_lxf_ref/sinChu7.fa /home/xjh3/data/ach/guiyu_gwas/hzau_lxf_ref/hzau_gw.fa_index

for ((j=1; j<=8; j++))
do
cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/${j}/${j}_radtags
	for i in *.rem.1.fq.gz
	do
	i=${i%.rem.1.fq.gz*}
	bowtie2 -x /home/xjh3/data/ach/guiyu_gwas/hzau_lxf_ref/hzau_gw.fa_index -1 ${i}.1.fq.gz -2 ${i}.2.fq.gz -S ${i}.sam -p 12 >${i}.log 2>&1
	samtools view -bS -t /home/xjh3/data/ach/guiyu_gwas/hzau_lxf_ref/hzau_gw.fa_index ${i}.sam > ${i}.bam
	samtools sort ${i}.bam -o ${i}.sort.bam
	rm ${i}.sam
	rm ${i}.bam
	done	
done

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/1/1_radtags

i=1
cat list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/2/2_radtags

i=25
cat /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/1/1_radtags/list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/3/3_radtags
i=49
cat /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/1/1_radtags/list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/4/4_radtags
i=73
cat /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/1/1_radtags/list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/5/5_radtags
i=97
cat /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/1/1_radtags/list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/6/6_radtags
i=121
cat /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/1/1_radtags/list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/7/7_radtags
i=145
cat /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/1/1_radtags/list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/8/8_radtags
i=169
cat /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/1/1_radtags/list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done



mv *.sort.bam ../02.sort.bam/

for ((j=1; j<=8; j++))
do
cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/00.CleanData/${j}/${j}_radtags
cp *.sort.bam ~/data/ach/guiyu_gwas/20230526_ldr_gwas/01.sort.bam
done

bcftools mpileup -Ou -f /home/xjh3/data/ach/guiyu_gwas/hzau_lxf_ref/sinChu7.fa  *.sort.bam | \
    bcftools call -Ou -mv | \
bcftools filter -s LowQual -e '%QUAL<10 || DP>2000' > ach_TEST_var_flt_qual10_dp30000.vcf


###传统艺能
awk '{if (/^#/ || $6>10) print }' ach_TEST_var_flt_qual10_dp30000.vcf > ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf

awk 'END{print NR}' ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf
277050

sed 's/;/ /g' ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf > ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf2
sed 's/DP=/DP= /g' ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf2 > ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf3

####每个样品平均4个READS x 181 =724#########
awk '{if ($9>724) print }' ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf3 > ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000.vcf3
sed 's/[0-9]*\,[0-9]*\,[0-9]*\:[0-9]*:[0-9]\t/replace2\t/g'  ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000.vcf3 > tt2
sed 's/:0:/:?:/g' tt2 > tt3
sed 's/:1:/:?:/g' tt3 > tt4
sed 's/:2:/:?:/g' tt4 > tt5
sed 's/:3:/:?:/g' tt5 > tt6
sed 's/:4:/:?:/g' tt6 > tt7
cp tt7 ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4.vcf3
awk 'END{print NR}' ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4.vcf3
87255
####基因型PL 得分>9#####
awk '{print $0,"xxx"}' ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4.vcf3 >tt8
sed 's/ /\t/g' tt8 >tt8a
sed 's/:[0-9]*\,[0-9]*\,[1-9]\t/:replace2\t/g'  tt8a > tt8b
sed 's/:[1-9]\,[0-9]*\,[0-9]*\t/:replace2\t/g' tt8b >tt9
sed 's/:[0-9]*\,[1-9]\,[0-9]*\t/:replace2\t/g' tt9 >tt10
sed 's/:0,0,0\t/:replace2\t/g' tt10 >tt11
sed 's/0\/0:replace2/?\t/g' tt11 > tt12
sed 's/0\/1:replace2/?\t/g' tt12 > tt13
sed 's/1\/1:replace2/?\t/g' tt13 > tt14
sed 's/.\/.:replace2/?\t/g' tt14 > tt15
sed 's/:[0-9]*\,[0-9]*\,[0-9]*/\t/g'  tt15 > tt16
cp tt16 ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4.vcf4

###########统计缺失数据的比例##############
awk '{print gsub("?","?"), $0}' tt16 >tt17
grep -v 'INDEL'  tt17 > tt18
awk 'END{print NR}' tt18
19757
cp tt18 ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_noindel.vcf5
awk '{print $1/181,$0}' tt18 > tt19
awk '{if($1<=0.25) print $0}' tt19 > tt20
awk 'END{print NR}' tt20
7469

awk -F\0/0 '{print "0/0",NF-1,$0}' tt20 > tt21
awk -F\0/1 '{print "0/1",NF-1,$0}' tt21 > tt22
awk -F\1/1 '{print "1/1",NF-1,$0}' tt22 > tt23
awk '{print 181-$8,$0}' tt23 >tt24
awk '{if($3/$1<0.9 && $5/$1<0.9 && $7/$1<0.9) print $0}' tt24 > ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss25_af09.vcf4
awk 'END{print NR}' ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss25_af09.vcf4
3828
cp ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss25_af09.vcf4 tt24a
grep -v ",A" tt24a >tt25
grep -v ",G" tt25>tt26
grep -v ",C" tt26>tt27
grep -v ",T" tt27>tt28
cp tt28 ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss25_af09_final.vcf4 
awk 'END{print NR}' tt28
3598


awk '{if($1<=0.1) print $0}' tt19 > bb20
awk 'END{print NR}' bb20
21935
awk -F\0/0 '{print "0/0",NF-1,$0}' bb20 > bb21
awk -F\0/1 '{print "0/1",NF-1,$0}' bb21 > bb22
awk -F\1/1 '{print "1/1",NF-1,$0}' bb22 > bb23
awk '{print 181-$8,$0}' bb23 >bb24
awk '{if($3/$1<0.9 && $5/$1<0.9 && $7/$1<0.9) print $0}' bb24 > ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss10_af09.vcf4

awk 'END{print NR}' ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss10_af09.vcf4
7431
cp ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss10_af09.vcf4 bb24a
grep -v ",A" bb24a >bb25
grep -v ",G" bb25>bb26
grep -v ",C" bb26>bb27
grep -v ",T" bb27>bb28

cp bb28 ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss10_af09_final.vcf4 
awk 'END{print NR}' bb28
6106
```



```{bash eval =F}
##主成分成分选择
###twstats 方法（使用tassel将hapmap文件转化为map和ped文件（plink格式）(先导入基因型文件-----save as plink文件)）
plink --file ldr.plk --pca 50 --out ach_pca --allow-extra-chr#因为存在scaffold等所以需要加--allow-extra-chr
conda activate ach
conda install eigensoft
twstats -t twtable -i ach_pca.eigenval -o ach_pca_50
less ach_pca_50 #选择 P < 0.05 的前 K个主成分
###基于可解释方差###
#利用ped map 文件生成bed bim fam 文件
plink --file ldr_trans.plk --out ldr_trans.plk --allow-extra-chr

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("gdsfmt")

library(gdsfmt)
library(SNPRelate)
bed.fn <- "ach_change.plk.bed"
fam.fn <- "ach_change.plk.fam"
bim.fn <- "ach_change.plk.bim"

# 将 PLINK 文件转为 GDS 文件
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "test.gds")

# 读取 GDS 文件
genofile <- snpgdsOpen("test.gds")

# 根据 LD 过滤 SNPs，阈值根据需要设定
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)

# 选择 SNP pruning 后要保留的 SNP
snpset.id <- unlist(unname(snpset))

# 计算 PCA，num.thread 是并行的线程数
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=10)

# 以百分比形式输出 variance proportion
print(pca$varprop*100)

library(ggplot2)
K= 30
qplot(x = 1:K, y = (pca$varprop[1:K]), col = "red", xlab = "PC", ylab = "Proportion of explained variance") + 
  geom_line() + guides(colour = FALSE) +
  ggtitle(paste("Scree Plot - K =", K))

```



```{bash eval =F}
#传统艺能结果很烂，选择使用plink+gemma分析 尝试
1.创建文件夹
mkdir /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink
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
   477382 test.map
      181 test.ped
   477563 total
可以看出，共有 181 个基因型个体，共有 477563 个 SNP 数据

#对 SNP 缺失率进行筛选
plink --file test --missing --chr-set 50  --allow-extra-chr 

##先过滤 SNP在个体中缺失率高于 10%的 SNP
plink --bfile g --geno 0.10 --make-bed --out g_2 --allow-extra-chr --chr-set 50 

plink --bfile g_2 --recode --out test_2 --allow-extra-chr --chr-set 50 
wc -l test_2.map test_2.ped
   18169 test_2.map
     181 test_2.ped
   18350 total
可以看出，过滤了 477563-18169 = 459394 个位点

#对样本缺失率进行筛选
##过滤缺失率高于10%的个体
plink --bfile g_2 --mind 0.10 --make-bed --out g_3 --allow-extra-chr --chr-set 50
plink --bfile g_3 --recode --out test_3 --allow-extra-chr --chr-set 50
wc -l test_3.map test_3.ped
   18169 test_3.map
     159 test_3.ped
   18328 total
可以看出，过滤了 181 - 159 = 22个体
#去掉 MAF 小于 0.05 的位点
plink --bfile g_3 --maf 0.05 --make-bed --out g_4 --allow-extra-chr --chr-set 50
plink --bfile g_4 --recode --out test_4 --allow-extra-chr --chr-set 50
wc -l test_4.map test_4.ped
   7753 test_4.map
    159 test_4.ped
   7912 total
可以看出，过滤了 18169-7753 = 10416个位点

#数据质控：哈温平衡
##计算所有位点的 HWE 的 P 值
plink --bfile g_4 --hardy --allow-extra-chr --chr-set 50
##提取哈温 p 值小于 0.0001 的位点
awk '{if($9 < 0.0001) print $0}' plink.hwe >plinkzoomhwe.hwe
##设定过滤标准 1e-4
plink --bfile g_4 --hwe 1e-4 --make-bed --out g_5 --allow-extra-chr --chr-set 50
plink --bfile g_5 --recode --out test_5 --allow-extra-chr --chr-set 50 
wc -l test_5.map test_5.ped
   7274 test_5.map
    159 test_5.ped
   7433 total
可以看出，过滤了 7753-7274 = 479个位点

#Step6_转化成vcf
plink --bfile g_5 --allow-extra-chr --recode vcf-iid --out ach_TEST_var_flt_all_qual10_dp30000_plink_s6


#Step7_过滤indel/SNP
vcftools --vcf ach_TEST_var_flt_all_qual10_dp30000_plink_s6.vcf --remove-indels --recode --recode-INFO-all --out ach_TEST_var_flt_all_qual10_dp30000_plink_s7_SNP

vcftools --vcf ach_TEST_var_flt_all_qual10_dp30000_plink_s6.vcf --keep-only-indels --out ach_TEST_var_flt_all_qual10_dp30000_plink_s7_indel --recode --recode-INFO-all

#Step8_删除scaffold
sed 's/.sort.bam//g' ach_TEST_var_flt_all_qual10_dp30000_plink_s7_SNP.recode.vcf > ach_TEST_var_flt_all_qual10_dp30000_plink_s8_SNP.vcf

sed 's/sinChu7-LG0//g' ach_TEST_var_flt_all_qual10_dp30000_plink_s8_SNP.vcf > ach_TEST_var_flt_all_qual10_dp30000_plink_s8_0_SNP.vcf

sed 's/sinChu7-LG//g' ach_TEST_var_flt_all_qual10_dp30000_plink_s8_0_SNP.vcf > ach_TEST_var_flt_all_qual10_dp30000_plink_s8_1_SNP.vcf

awk '{ if(/^#/ ||$1 >=1 && $1 <= 24) print $0}' ach_TEST_var_flt_all_qual10_dp30000_plink_s8_1_SNP.vcf >ach_TEST_var_flt_all_qual10_dp30000_plink_s8_2_SNP.vcf
#awk /^#/为正则匹配，将匹配文件中所有是#开头 https://blog.csdn.net/XiaodunLP/article/details/90262319


#Step9_给VCF文件加上ID
perl /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/VCF_add_id-master/VCF_add_id.pl ach_TEST_var_flt_all_qual10_dp30000_plink_s8_2_SNP.vcf ach_TEST_var_flt_all_qual10_dp30000_plink_s9_SNP.vcf

#step10_转换成bed文件
plink --vcf ach_TEST_var_flt_all_qual10_dp30000_plink_s9_SNP.vcf --make-bed --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP


plink --bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --recode --out ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP --allow-extra-chr --chr-set 50 
wc -l ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP.map ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP.ped

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight_transform

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -gk 2 -o kin

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -k /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight_transform/output/kin.sXX.txt -lmm 4 -o yield

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -k /home/xiajh2/data/ach/EC-gwas-plink/output/kin.sXX.txt -lmm 1 -o yield_s10

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -gk 1 -o kin_1

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight_transform/output/kin_1.cXX.txt


/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -k /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight_transform/output/kin_1.cXX.txt -lmm 4 -o yield_1



cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -gk 1 -o kin_1

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -k /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight_transform/output/kin_1.cXX.txt -lmm 4 -o yield_1


cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodylength

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -gk 1 -o kin_1

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -k /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight_transform/output/kin_1.cXX.txt -lmm 4 -o yield_1

cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/totallength

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -gk 1 -o kin_1

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -k /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight_transform/output/kin_1.cXX.txt -lmm 4 -o yield_1


cd /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyheight

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -gk 1 -o kin_1

/home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/gemma-0.98.5-linux-static-AMD64 -bfile ach_TEST_var_flt_all_qual10_dp30000_plink_s10_SNP -k /home/xjh3/data/ach/guiyu_gwas/20230526_ldr_gwas/02.plink/traits/bodyweight_transform/output/kin_1.cXX.txt -lmm 4 -o yield_1
```

