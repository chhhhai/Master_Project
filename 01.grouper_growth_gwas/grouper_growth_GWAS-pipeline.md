# 斜带石斑鱼苗期生长GWAS数据处理

## 1.12G数据预处理

```{bash eval=FALSE}
bowtie2-build /home/xiajh2/data/ach/ncbi_data/ref/Eco_Chr_anchored20170428.fa /home/xiajh2/data/ach/ncbi_data/ref/Eco_Chr.fa_index

for i in *_1.clean.fq.gz
do
i=${i%_1.clean.fq.gz*}
bowtie2 -x /home/xjh3/data/ach-EC-rna-seq/ref/Eco_Chr.fa_index -1 ${i}_1.clean.fq.gz -2 ${i}_2.clean.fq.gz -S ${i}.sam -p 12 >${i}.log 2>&1
done

#拆分数据
for i in *_1.clean.fq.gz
do
i=${i%_1.clean.fq.gz*}
mkdir ${i}
/software/stacks-2.4/process_radtags -1 ${i}_1.clean.fq.gz -2 ${i}_2.clean.fq.gz -o ./${i}/ -b AD_seq_barcode.txt -r -c -q --renz_1 ecoRI --renz_2 mspI -i gzfastq > ${i}.radtags.log 2>&1
done
```

# 2.120G数据处理

```{bash eval=FALSE}
操作同上
#生成sam文件
cd /home/xjh3/data/ach-EC-gwas/01.clean_data/ach-EC-growth-ach-EC-growth1_BKDL210048487-1a-1/
for i in *.rem.1.fq.gz
do
i=${i%.rem.1.fq.gz*}
bowtie2 -x /home/xjh3/data/ach-EC-rna-seq/ref/Eco_Chr.fa_index -1 ${i}.1.fq.gz -2 ${i}.2.fq.gz -S ${i}.sam -p 1 >${i}.log 2>&1
done

for i in *.rem.1.fq.gz
do
i=${i%.rem.1.fq.gz*}
samtools view -bS -t /home/xjh3/data/ach-EC-rna-seq/ref/Eco_Chr.fa_index ${i}.sam > ${i}.bam
done

for i in *.rem.1.fq.gz
do
i=${i%.rem.1.fq.gz*}
samtools sort ${i}.bam -o ${i}.sort.bam
done

for i in *.rem.1.fq.gz
do
i=${i%.rem.1.fq.gz*}
samtools view -bS -t /home/xjh3/data/ach-EC-rna-seq/ref/Eco_Chr.fa_index ${i}.sam > ${i}.bam
samtools sort ${i}.bam -o ${i}.sort.bam
done

#批量修改练习数据
awk '{print "sample_"$1".sam"}' barcode.txt > sample.txt
awk '{print "touch "$1}' sample.txt >1.sh
sh 1.sh
i=1
cat sample.txt | while read id
do
mv $id s$i.sam
let i=i+1
done

#批量修改每个子库的barcode文件（sample_xxxxx.sam为例子）
sed -i "s/\s*$//g" barcode.txt ##去除末尾特殊字符^M
awk '{print "sample_"$1".sam"}' barcode.txt > sample.txt
sed  's/.sam/.sort.bam/g' sample.txt>list.txt

i=265
cat list.txt | while read id
do
mv $id s$i.sort.bam
let i=i+1
done

mv *.sort.bam ../02.sort.bam/

#https://blog.csdn.net/chengsi9809/article/details/100701566
bcftools mpileup -Ou -f /home/xjh3/data/ach-EC-rna-seq/ref/Eco_Chr_anchored20170428.fa *.sort.bam | \
    bcftools call -Ou -mv | \
bcftools filter -s LowQual -e '%QUAL<10 || DP>2000' > ach_TEST_var_flt_qual10_dp30000.vcf
参数注释：
bcftools mpileup -Ou -f
#mpileup 用pileup方法监测变异
#-Ou 参数可以在bcftools管道符命令中使用。添加后，可以减少VCF、BCF格式间的转换，加速处理过程
#-f 参数后添加fa格式的参考基因组文件
bcftools call -Ou -mv 
#-m bcftools优化后的变异检测参数
#-v 只输出有变异的位点
bcftools filter -s LowQual -e '%QUAL<10 || DP>2000' > ach_TEST_var_flt_qual10_dp30000.vcf
#-s 满足筛选条件注释为特定字符
#-e 添加表达式命令 添加筛选条件
#这个命令意思是将QUAL值小于10，深度高于2000位点的FILTER列标记为LowQual
##https://www.jianshu.com/p/5781e7d74c40  去重复

##snpEff使用----https://www.jianshu.com/p/e03095642ef0----https://www.jianshu.com/p/77c3a2fae4ab
xjh3@xjh3:/software/snpEff$ java -jar snpEff.jar build -gtf22 -v Eco_Chr_anchored20170428

java -jar snpEf/snpEff.jar -v Eco_Chr_anchored20170428 Ecoli.vcf.gz > Ecoli.anno.vcf.gz
java -Xmx16g -jar /software/snpEff/snpEff.jar eff -v Eco_Chr_anchored20170428 -i bed  ach_reduce_gwas.bed >snpeff.log 2>&1
```

![1213213](C:\Users\83617\Desktop\扩增子测序材料\1213213.png)

```{bash eval=FALSE}
####筛选genotyping quality >10#############
awk '{if (/^#/ || $6>10) print }' ach_TEST_var_flt_qual10_dp30000.vcf > ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf

awk 'END{print NR}' ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf
826525
sed 's/;/ /g' ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf > ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf2
sed 's/DP=/DP= /g' ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf2 > ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf3

####每个样品平均4个READS x 250 =1000#########
awk '{if ($9>1000) print }' ach_TEST_var_flt_qual10_dp30000_filter_q10.vcf3 > ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000.vcf3
sed 's/[0-9]*\,[0-9]*\,[0-9]*\:[0-9]*:[0-9]\t/replace2\t/g'  ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000.vcf3 > tt2
sed 's/:0:/:?:/g' tt2 > tt3
sed 's/:1:/:?:/g' tt3 > tt4
sed 's/:2:/:?:/g' tt4 > tt5
sed 's/:3:/:?:/g' tt5 > tt6
sed 's/:4:/:?:/g' tt6 > tt7
cp tt7 ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4.vcf3
awk 'END{print NR}' ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4.vcf3
237395
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
154781
cp tt18 ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_noindel.vcf5
awk '{print $1/250,$0}' tt18 > tt19
awk '{if($1<=0.25) print $0}' tt19 > tt20
awk 'END{print NR}' tt20
45519

awk '{print $1/288,$0}' tt18 > tt19_1
awk '{if($1<=0.25) print $0}' tt19_1 > tt20_1
awk 'END{print NR}' tt20_1
51918

awk -F\0/0 '{print "0/0",NF-1,$0}' tt20_1 > tt21
awk -F\0/1 '{print "0/1",NF-1,$0}' tt21 > tt22
awk -F\1/1 '{print "1/1",NF-1,$0}' tt22 > tt23
awk '{print 250-$8,$0}' tt23 >tt24
awk '{if($3/$1<0.9 && $5/$1<0.9 && $7/$1<0.9) print $0}' tt24 > ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss25_af09.vcf4
awk 'END{print NR}' ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss25_af09.vcf4
15071
cp ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss25_af09.vcf4 tt24a
grep -v ",A" tt24a >tt25
grep -v ",G" tt25>tt26
grep -v ",C" tt26>tt27
grep -v ",T" tt27>tt28
cp tt28 ach_TEST_var_flt_qual10_dp30000_filter_q10_dp1000_dp4_gp4_miss25_af09_final.vcf4 
awk 'END{print NR}' tt28
13746


awk '{if($1<=0.1) print $0}' tt19 > bb20
awk 'END{print NR}' bb20
21935
awk -F\0/0 '{print "0/0",NF-1,$0}' bb20 > bb21
awk -F\0/1 '{print "0/1",NF-1,$0}' bb21 > bb22
awk -F\1/1 '{print "1/1",NF-1,$0}' bb22 > bb23
awk '{print 250-$8,$0}' bb23 >bb24
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




i=1
cat list.mapping.txt | while read id
do
cp $id S$i.log
let i=i+1
done


##主成分成分选择
###twstats 方法（使用tassel将hapmap文件转化为map和ped文件（plink格式）(先导入基因型文件-----save as plink文件)）
plink --file ach_change.plk --pca 50 --out ach_pca --allow-extra-chr#因为存在scaffold等所以需要加--allow-extra-chr
conda activate ach
conda install eigensoft
twstats -t twtable -i ach_pca.eigenval -o ach_pca_50
less ach_pca_50 #选择 P < 0.05 的前 K个主成分
###基于可解释方差###
#利用ped map 文件生成bed bim fam 文件
plink --file ach-genotype_begin_s100.plk --out ach-genotype_begin_s100.plk --allow-extra-chr

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



```{bash eval = F}
#采用徐洲更方法

cd ~/data/ach/ach-EC-gwas/01.clean_data/02.sort.bam
cd xuzhougen_pipeline/

/home/xjh3/data/hdd/software/seqkit/seqkit seq -n genome.fa|grep 'Chr' |head -n 10| \
while read region
do
bcftools mpileup -f /home/xjh3/data/ach/guiyu_gwas/ref/genome.fa  \
     --regions ${region} \
     -Ou /home/xjh3/data/ach/guiyu_gwas/DNA/*.sort.sambamba_markdup.bam | \
     bcftools call -mv -Ob -o /home/xjh3/data/ach/guiyu_gwas/DNA/${region}.bcf &
done

/home/xjh3/data/hdd/software/seqkit/seqkit seq -n genome.fa|grep 'Chr' |sed -n '11,25p'| \
while read region
do
bcftools mpileup -f /home/xjh3/data/ach/guiyu_gwas/ref/genome.fa  \
     --regions ${region} \
     -Ou /home/xjh3/data/ach/guiyu_gwas/DNA/*.sort.sambamba_markdup.bam | \
     bcftools call -mv -Ob -o /home/xjh3/data/ach/guiyu_gwas/DNA/${region}.bcf &
done


bcftools concat --naive -o merged.bcf *.bcf

bcftools filter  -g3 -G10 -e'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || MQ < 30 || MQSB <=0.1'  merged.bcf > filter.vcf

```

