# lncRNA鉴定

## RNA-seq数据下载

1.登陆ncbi网站，进入SRA数据库，输入对应物种的拉丁文，如斜带石斑鱼（*Epinephelus coioides*）

![](C:\Users\83617\Desktop\lncRNA\QQ截图20210411221418.png)

2.选择RNA，导出SRA数据库所有的与物种相关的RNA-seq，Run Selector ！

![](C:\Users\83617\Desktop\lncRNA\2.png)



3.通过得到的数据的表格选择自己需要的数据，通过Accession List生成的文件统计自己rna-seq的run号,便于后期批量处理

![](C:\Users\83617\Desktop\lncRNA\3.png)

![](C:\Users\83617\Desktop\lncRNA\4.png)

4.利用得到的Run（例如：SRR10058572）或Bioproject（例如：PRJNA563518），得到aspera格式的下载链接

4.1 通过https://sra-explorer.info/?#输入run，生成aspera格式的下载链接（网站不稳定，一般需要翻墙）

![](C:\Users\83617\Desktop\lncRNA\5.png)

4.2通过https://www.ebi.ac.uk/ena/browser/home进入ENA数据库，利用Bioproject，生成生成aspera格式的下载链接

![](C:\Users\83617\Desktop\lncRNA\7.png)

![](C:\Users\83617\Desktop\lncRNA\8.png)

5.将得到的aspera格式的下载链接制作成sra_download.txt

![](C:\Users\83617\Desktop\lncRNA\9.png)

6.批量生成ascp下载命令

注：有时候将下载链接复制到txt或excel行尾会出现特殊符号，导致运行命令报错，需要将sra_download.txt存在的特殊符号删除（cat -A 可检查文件是否存在特殊符号）

```{bash eval=FALSE}
sed -i "s/\s*$//g" sra_download.txt
```

6.1利用awk生产成批量处理的脚本

```{bash eval=FALSE}
awk '{print "ascp -QT -k 1 -l 300m -P33001 -i /home/data/gmb45/miniconda3/etc/asperaweb_id_dsa.openssh era-fasp@" $1 " /home/data/gmb45/pracice_data}' sra.download.txt >sra.download.sh

bash sra_download.sh
```

6.2利用循环语句进行批量处理

```{bash eval=FALSE}
cat sra_download.txt |while read id
do
ascp -QT -l 300m -P33001 \
-i /home/xiajh2/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@$id /home/xiajh2/data/ach/ncbi_data/1.raw_data/
done>sra_download.sh

nohup bash sra.download.sh >sra.download.log 2>&1 &
```

7.得到fastq.gz文件

```{bash eval=FALSE}
ls `pwd`/rawdata/*fastq.gz > 2.fastp_data/fastq.list
```

#

## 数据过滤，比对和组装（fastp+hisat2+stringtie）

1.建立hisat2的索引

```{bash eval=FALSE}
# make exon 
hisat2_extract_exons.py Eco_chr20170428.gtf > Eco_chr20170428.exon

# make splice site
hisat2_extract_splice_sites.py Eco_chr20170428.gtf > Eco_chr20170428.ss

# make snp and haplotype(鱼基本无这个，snp.txt需要从数据库下载)
hisat2_extract_snps_haplotypes_UCSC.py Eco_Chr_anchored20170428.fa snp151Common.txt snp151Common 

# build index
hisat2-build -p 6 --snp snp151Common.snp --haplotype snp151Common.haplotype --exon Eco_chr20170428.exon  --ss Eco_chr20170428.ss Eco_Chr_anchored20170428.fa Eco_Chr_anchored20170428.fa

```



2.批量运行fastp+hisat2+stringtie

```{bash eval=FALSE}
cat fastq.list| while read id
do
/home/data/fastp -i ${id}_1.fastq.gz -o ../2.fastp_data/${id}_1.fastp.fq.gz \
-I ${id}_2.fastq.gz \
-O ../2.fastp_data/${id}_2.fastp.fq.gz \
-q 20 -w 12 -n 6 --compression=6 \
-h ../2.fastp_data/${id}.html -j ../2.fastp_data/${id}.json
hisat2 -p 8 \
-x ../ref/Eco_Chr_anchored20170428.fa \
-1 ../2.fastp_data/${id}_1.fastp.fq.gz \
-2 ../2.fastp_data/${id}_2.fastp.fq.gz \
-S ../3.bam/${id}.sam > ../3.bam/${id}.log 2>&1
samtools sort -O BAM -o ../3.bam/${id}.sort.bam -@6 -T ../3.bam/${id}.sort.bam.temp ../3.bam/${id}.sam
rm ../3.bam/${id}.sam
stringtie ../3.bam/${id}.sort.bam -p 4 -G /home/data/ref/Eco_chr20170428.gtf -o ../4.stringtie_gtfs/${id}.gtf -l ${id} -B
mkdir ../4.stringtie_gtfs/${id}_ballgown
mv ../4.stringtie_gtfs/*.ctab ../4.stringtie_gtfs/${id}_ballgown
done
```

## 合并gtf后鉴定和筛选新的lncRNA

1.合并gtf

```{bash eval=FALSE}
ls `pwd`/4.stringtie_gtfs/*.gtf >mergelist_117samples.txt
stringtie --merge -p 10 -G /home/xjh3/data/ach-EC-rna-seq/ref/Eco_chr20170428.gtf -o ../5.merge_gtf/merge_gtf ../4.stringtie_gtfs/mergelist_117samples.txt > ../4.stringtie_gtfs/merge_gtf.log 2>&1
```

1.1Gffcompare 获取转录本组装情况---------(2021-12-14)

```{bash eval=FALSE}
/software/gffcompare/gffcompare-0.12.1.Linux_x86_64/gffcompare -r /home/xiajh2/data/ach/ncbi_data/ref/Eco_Chr20170428.gff -o merged merge_gtf ##生成了merged.annnoted.gtf
```

![1639467342(1)](C:\Users\83617\Desktop\lncRNA\1639467342(1).png)

```{bash eval=FALSE}
# 统计class code 类型
awk '$3!~/class/ {print $3}' merged.merge_gtf.tmap | sort -V | uniq -c
```

![1639469416(1)](C:\Users\83617\Desktop\lncRNA\1639469416(1).png)





2.鉴定和筛选新的lncRNA

```{bash eval=FALSE}
#初步鉴定新的转录本
cd 5.merge_gtf
/home/xjh3/data/hdd_salinity_ncbi_sra/software/gffcompare/gffcompare-0.12.1.Linux_x86_64/gffcompare -r /home/xjh3/data/ach-EC-rna-seq/ref/Eco_Chr20170428.gff -o merged merge_gtf ##生成了merged.annnoted.gtf
awk '$3=="u"{print $0}' merged.merge_gtf.tmap >filter_by_u.tmap
awk '$3=="u"{print $0}' merged.merge_gtf.tmap|wc #30683个转录本
#过滤，只保留exon>1且长度>200bp的转录本
awk '($6>1 && $10>=200){print $0}' filter_by_u.tmap > filter2_by_exon_length.tmap
wc -l filter2_by_exon_length.tmap #15127个转录本
```

3.提取出新的lncRNA的fasta文件

```{bash eval=FALSE}
#获取过滤出的转录本的ID
awk '{print $5}' filter2_by_exon_length.tmap >filter2_transcript_ID
#得到gtf
grep -Ff filter2_transcript_ID -w merged_annotated.gtf >filter2_transcript.gtf
#利用exon组成gtf
awk '$3=="exon"{print $0}' filter2_transcript.gtf >filter2_transcript_exon.gtf
#根据exon位置信息提取基因组序列，组装成转录本序列
gffread -w filter2_transcript_exon.fa -g ../ref/Eco_Chr_anchored20170428.fa filter2_transcript_exon.gtf
grep ">"filter2_transcript_exon.fa |wc
15127
```

问：可以直接先把整个gtf文件转成fa文件，利用seqtk直接利用ID提取出fa？(个人觉得seqkit比blast好用)

```{bash eval=FALSE}
seqtk subseq merge.fasta filter2_transcript_ID >filter2_transcript_exon.fasta
```

##鉴定出所有的lncRNA的转录本

1.将merged_gtf转化成merge.fasta

```{bash eval=FALSE}
gffread -w merge.fasta -g ../ref/Eco_Chr_anchored20170428.fa merge_gtf
grep ">" merge.fasta  > merge_fasta.name.list
157452

#过滤大于200bp的转录本
samtools faidx merge_rm_NN.fasta
awk '{if($2>=200){print $1}}' merge_rm_NN.fasta.fai >merge.fasta.name.list
awk 'END{print NR}' merge.fasta.name.list
157363
/home/xiajh2/data/ach/lncrna-EC/merge_data/seqtk   subseq merge.fasta merge.fasta.name.list >merge_rm.fasta
#去除含有N端序列
sed 'N;s/\n/_/' merge_rm.fasta | grep -v N | tr "_" "\n" > merge_rm_NN.fasta
cat merge_rm.fasta |paste - - |grep -v N | tr "\t" "\n" | more
136724
#cpc
/software/CPC2_standalone-1.0.1/bin/CPC2.py -i merge_rm_NN.fasta -o cpc2output
grep "noncoding" cpc2output.txt |awk '{print $1}'>cpc.list
wc cpc.list
47782
#cnci
sudo /usr/bin/python2.7 /home/xiajh2/data/hdd/RNA-seq/software/CNCI-master/CNCI.py -f ../merge_rm_NN.fasta -o CNCI_fasta -m ve -p 20 > CNCI_fasta.log 2>&1
mv CNCI.index CNCI_index.list
grep "noncoding" CNCI_index.list | awk '{print $1}' > CNCI_noncoding.list
wc CNCI_noncoding.list
45622
#plek
python /software/PLEK.1.2/PLEK.py -fasta merge_rm_NN.fasta -out merge_rm_NN_PLEK.fasta -thread 2
grep "Non-coding" merge_rm_NN_PLEK.fasta |awk '{print $3}'>PLEK.list
sed 's/>//g' PLEK.list >plek.list
wc plek.list
49485
```

2.取cpc，cnci以及plek结果的交集

![](C:\Users\83617\Desktop\lncRNA\cpc_plek_cnci.png)

```{bash eval=FALSE}
sort CNCI_noncoding.list plek.list | uniq -d > unite_cnci_plek_lncRNA.list
wc -l unite_cnci_plek_lncRNA.list
37357
sort unite_cnci_plek_lncRNA.list cpc.list | uniq -d > unite_cnci_plek_cpc2_lncRNA.list
wc -l unite_cnci_plek_cpc2_lncRNA.list
34944
/home/xiajh2/data/ach/lncrna-EC/merge_data/seqtk subseq merge_rm_NN.fasta unite_cnci_plek_cpc2_lncRNA.list >plek_CNCI_cpc2.fasta
```

3.与nr，uniref90 ，pfam数据库比对

```{bash eval=FALSE}
/home/xiajh2/data/hdd/RNA-seq/tilapia_uniprot_Nr_proteome/diamond blastx -d /home/xiajh2/data/hdd/RNA-seq/tilapia_uniprot_Nr_proteome/nr_db.fasta -q ./plek_CNCI_cpc2.fasta -p 8 -o nr_fmt6.txt > diamond.log 2>&1


/home/xiajh2/data/hdd/RNA-seq/tilapia_uniprot_Nr_proteome/diamond blastx -d /home/xiajh2/data/hdd/RNA-seq/tilapia_uniprot_Nr_proteome/uniref90_db.fasta -q./plek_CNCI_cpc2.fasta -p 8 -o uniref90_fmt6.txt > uniref90.log 2>&1


/software/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t plek_CNCI_cpc2.fasta -m 20
pfam_scan.pl -fasta /home/xiajh2/data/ach/ncbi_data/pfam_nr_uniref90/plek_CNCI_cpc2.fasta.transdecoder_dir/longest_orfs.pep -dir /home/xiajh2/data/hdd/RNA-seq/Pfam/dir_PfamA.hmm  -outfile Pfam_hits_transdecoder_m20_lncRNA.fasta            -cpu 16 > Pfam_hits.log 2>&1

#过滤掉大于1e10-05的转录本
awk '$11<=1e-05{print $1}' nr_fmt6.txt |sort|uniq > nr_out.list
wc -l nr_out.list
13812
awk '$11<=1e-05{print $1}' uniref90_fmt6.txt |sort|uniq > uniref90_out.list
wc -l uniref90_out.list
14360
sort nr_out.list uniref90_out.list | uniq -d > unite_nr_uniref_lncRNA.txt           awk 'END{print NR}' unite_nr_uniref_lncRNA.txt
13592
grep -v '#' Pfam_hits_transdecoder_m20_lncRNA.fasta |awk '$13<=1e-05{print $1 }'|sed 's/...$//g'|uniq >pfam_treat_coding.txt #存在空白行已删除
sed 's/\.$//' pfam_treat_coding.txt >pfam_treat_coding_1.txt
sort pfam_treat_coding_1.txt|uniq >pfam_treat_coding_2.txt
wc -l pfam_treat_coding_2.txt
4139
sort unite_nr_uniref_lncRNA.txt pfam_treat_coding_2.txt| uniq -d > unite_nr_uniref_pfam_lncRNA.txt
wc -l unite_nr_uniref_pfam_lncRNA.txt
3997
```

![](C:\Users\83617\Desktop\lncRNA\pfam_nr_uniref.png)

4.将低于1e10-5的转录本（具有编码蛋白能力的转录本），从cpc，cnci以及plek结果的交集结果去除

```{bash eval=FALSE}
#取差集
sort unite_cnci_plek_cpc2_lncRNA.list  unite_nr_uniref_pfam_lncRNA.txt  unite_nr_uniref_pfam_lncRNA.txt |uniq -u> filter_database_noncoding.txt
wc -l filter_database_noncoding.txt
30997
```

5.去除ORF>100aa的转录本

```{bash eval=F}
grep -Ff filter_database_100aa_noncoding.txt -w /home/xiajh2/data/ach/ncbi_data/merge_gtf/merge_gtf > final_predicted_ORF_nocoding.gtf



gffcompare -r /home/xiajh2/data/ach/ncbi_data/ref/Eco_Chr20170428.gff -o final_predicted final_predicted_ORF_nocoding.gtf

awk '$3!~/class/ {print $3}' final_predicted.final_predicted_ORF_nocoding.gtf.tmap | sort -V | uniq -c

   2294 i
    407 j
     82 k
     10 m
     18 n
    211 o
    389 p
     44 s
  16385 u
   3644 x
    494 y
    588 =


#过滤，只保留class_code="u","x","i","j","o"的 transcripts
awk '{if ($3=="u" || $3=="x" || $3=="i" || $3=="j" || $3=="o"){print $0}}' final_predicted.final_predicted_ORF_nocoding.gtf.tmap  > filter1_by_uxijo.tmap
$ wc -l filter1_by_uxijo.tmap
22941 filter1_by_uxijo.tmap
#获取剩余transcripts的exons位置信息，提取序列并组装成转录本序列
# 获取剩余的transcripts的ID
awk '{print $5}' filter1_by_uxijo.tmap > filter1_transcript_ID
$ wc -l filter1_transcript_ID
22941 filter1_transcript_ID

# 剩余的transcripts得到gtf
grep -w -Ff filter1_transcript_ID -w final_predicted.annotated.gtf > filter1_transcript.gtf

# 把filter2_transcript.gtf中的class_code "=" 替换为L

# 去除剩余为去除的class_code "=" 
awk -F ';' '{if ($8!=" L"){print $0}}' filter1_transcript.gtf > filter1
mv filter1 filter1_transcript.gtf
wc -l filter1_transcript.gtf
62223 filter1_transcript.gtf
# 过滤，只保留exon>1并且长度>200bp的transcripts
awk '($10>=200){print$0}' filter1_by_uxijo.tmap > filter2_by_exon_length.tmap
wc -l filter2_by_exon_length.tmap

awk '{print $5}' filter2_by_exon_length.tmap > filter2_transcript_ID

grep -Ff filter2_transcript_ID -w filter1_transcript.gtf > filter.gtf

########过滤掉低表达量的lncRNA#########
ls /home/xiajh2/data/ach/ncbi_data/hisat_data|grep "85821"|grep "bam" >VNN.sample.txt

awk '{print "mv /home/xiajh2/data/ach/ncbi_data/hisat_data/" $1" ./"}' VNN.sample.txt >VNN.sample.sh

sh  VNN.sample.sh

featureCounts -T 8 -a \
/home/xiajh2/data/ach/ncbi_data/2021-12-19_from_ORF/filter.gtf  \
-o ./raw_count.txt -p -B -C -f -t transcript -g transcript_id  \
./*.bam > transcript_featureCount.log 2>&1

# R 语言计算FPKM，筛选：FPKM > 0 in at least one sample ，得到lncRNA_id.txt
rm(list=ls())
# make count table 
raw_df <- read.table(file = "~/lncRNA_project/test/08.featurecounts/raw_count.txt",header = T,skip = 1,sep = "\t")

count_df <- raw_df[ ,c(7:ncol(raw_df))]
metadata <- raw_df[ ,1:6] # 提取基因信息count数据前的几列
rownames(count_df) <- as.character(raw_df[,1])
colnames(count_df) <- paste0("SRR85821",75:83)

# calculate FPKM
countToFpkm <- function(counts, effLen)
{
  N <- colSums(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

options(scipen = 200) # 表示在200个数字以内都不使用科学计数法
fpkm = countToFpkm(count_df, metadata$Length)
# View(fpkm)

# FPKM > 0 in at least one sample
count_df.filter <- count_df[rowSums(fpkm)>0,]

write.table(rownames(count_df.filter),file="filter6_by_fpkm_id", sep="\t",quote=F) ##这个人写的有问题，已经修改
vi filter6_by_fpkm_id
awk '{print $2}' filter6_by_fpkm_id > filter_by_fpkm_id

sed -i "s/\s*$//g"   filter_by_fpkm_id

# linux 里提取最终lncRNA的gtf文件
grep -Ff filter_by_fpkm_id -w /home/xiajh2/data/ach/ncbi_data/merge_gtf/2021-12-14-final-transcripts/filter.gtf > lncRNA.gtf

# 根据exon位置信息提取基因组序列，组装成转录本序列
gffread -w lncRNA.fa -g /home/xiajh2/data/ach/ncbi_data/ref/Eco_Chr_anchored20170428.fa  lncRNA.gtf
$ grep -c '^>' lncRNA.fa
10408

####linux featurecounts定量####
#featureCounts-by-gff----mrna
featureCounts -t mRNA -g ID  \
-Q 10 --primary -s 0 -p -f -T 8 \
-a /home/xiajh2/data/ach/ncbi_data/ref/Eco_Chr20170428.gff \
-o ./mRNA_by_gff_count.txt \
/home/xiajh2/data/ach/ncbi_data/merge_gtf/2021-12-14-final-transcripts/featurecounts/*.bam \
> mRNA_by_gff_featureCounts.log 2>&1 
#featureCounts-by-gff----lncRNA
featureCounts -t transcript -g transcript_id  \
-Q 10 --primary -s 0 -p -f -T 8 \
-a ./lncRNA.gtf \
-o ./lncRNA_raw_count.txt \
./*.bam \
> lncRNA_featureCounts.log 2>&1
```



 ```{bash eval=F}
 #windows Rsudio进行DEseq2差异分析（生信宝典脚本）
 ##靶基因预测
 1.cis
 #lncRNA_low_control_id.txt由deseq2分析得到，用excel提取转录本id，vi生成
 grep -Ff lncRNA_low_control_id.txt -w lncRNA_raw_count.txt |awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > VNN_lncRNA_low_control.bed
 
 #mrna.bed获得
 grep -w 'mRNA' /home/xiajh2/data/ach/ncbi_data/ref/Eco_Chr20170428.gff|grep -v "C"|grep -v "S"|awk '{print $1 "\t" $4 "\t" $5 "\t" $9}'|sed 's/ID=//g'|sed 's/;//g'|less > mRNA_id.bed
 
 bedtools intersect -a mRNA_id.bed -b VNN_lncRNA_low_control.bed -wo > mRNA_overlapping_lncRNA.txt
 bedtools window -a mRNA_id.bed -b  VNN_lncRNA_low_control.bed -l 100000 -r 100000 > mRNA_up_down_stream_lncRNA.txt
 
 #DElncRNA和DEmRNA做相关性分析
 ###提取DElncRNA的表达量文件（FPKM）
 mv diff_merge.table diff_merge_lncRNA_low_control.table
 awk '{print $2 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14}' diff_merge_lncRNA_low_control.table>diff_merge_sort_lncRNA_low_control.table
 grep  -Ff lncRNA_low_control_id.txt -w diff_merge_sort_lncRNA_low_control.table |sort -k 1|uniq >lncRNA_low_control_FPKM.txt ##出现重复的原因不清楚，遍历出现了问题？
 
 ###提取cis_mRNA的表达量文件（FPKM）
 awk '{print $4}' mRNA_up_down_stream_lncRNA.txt >cis_mRNA_id.txt
 
 awk '{print $2 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14}' diff_merge_mRNA_low_control.table > diff_merge_sort_mRNA_low_control.table
 
 grep  -Ff cis_mRNA_id.txt -w diff_merge_sort_mRNA_low_control.table |sort -k 1|uniq >mRNA_low_control_FPKM.txt
 
 
 ###windows R 关联分析
 library(Hmisc)
 library(reshape2)
 library(fdrtool)
 
 options(scipen = 200)
 getwd()
 setwd('E:/中山大学研究生（2021）/实验/lncRNA-EC/VNN(SRR8582175-SRR8582183)/lncRNA-co-mRNA(cis)')
 mRNA<-read.table("mRNA_low_control_FPKM.txt",header =FALSE,check.names =FALSE,stringsAsFactors =FALSE,sep="\t",row.names=1)
 lncRNA<-read.table("lncRNA_low_control_FPKM.txt",header =FALSE,check.names =FALSE,stringsAsFactors =FALSE,sep="\t",row.names=1)
 n_lncRNA<-nrow(lncRNA)
 n_mRNA<-nrow(mRNA)
 alldata<-rbind(lncRNA,mRNA)
 mydata<-as.matrix(t(alldata))
 p<-rcorr(mydata)
 cor<-p[["r"]][1:n_lncRNA,(n_lncRNA+1):nrow(alldata)]
 pvalue<-p[["P"]][1:n_lncRNA,(n_lncRNA+1):nrow(alldata)]
 cor<-cbind(rownames(cor),cor)
 
 colnames(cor)[1]<-"LncRNA"
 pvalue<-cbind(rownames(pvalue),pvalue)
 colnames(pvalue)[1]<-"LncRNA"
 new_cor<-melt(as.data.frame(cor),id.vars="LncRNA",variable.name ="mRNA",value.name ="correlation")
 new_p<-melt(as.data.frame(pvalue),id.vars="LncRNA",variable.name ="mRNA",value.name ="p.value")
 
 
 all_cor<-merge(new_cor, new_p, by.x=c("LncRNA","mRNA"), by.y=c("LncRNA","mRNA"))
 all_cor<-na.omit(all_cor)
 all_cor$p.value<-as.numeric(all_cor$p.value)
 all_cor$correlation<-as.numeric(all_cor$correlation)
 
 
 write.table(all_cor,paste("sig_DElncRNAs_coexpression_DEmRNA_result.txt",sep=""),quote=FALSE,row.names=FALSE)
 
 all_fli<-subset(all_cor,abs(all_cor$correlation)>0.8&all_cor$p.value <0.001)
 
 
 2.trans
 
 得到所有mrna的表达量矩阵（生成all_mRNA_low_control_FPKM.txt文件）
 ###Linux R 关联分析
 library(Hmisc)
 library(reshape2)
 library(fdrtool)
 
 options(scipen = 200)
 getwd()
 mRNA<-read.table("mRNA_low_control_FPKM.txt",header =FALSE,check.names =FALSE,stringsAsFactors =FALSE,sep="\t",row.names=1)
 lncRNA<-read.table("lncRNA_low_control_FPKM.txt",header =FALSE,check.names =FALSE,stringsAsFactors =FALSE,sep="\t",row.names=1)
 n_lncRNA<-nrow(lncRNA)
 n_mRNA<-nrow(mRNA)
 alldata<-rbind(lncRNA,mRNA)
 mydata<-as.matrix(t(alldata))
 p<-rcorr(mydata)
 cor<-p[["r"]][1:n_lncRNA,(n_lncRNA+1):nrow(alldata)]
 pvalue<-p[["P"]][1:n_lncRNA,(n_lncRNA+1):nrow(alldata)]
 cor<-cbind(rownames(cor),cor)
 
 colnames(cor)[1]<-"LncRNA"
 pvalue<-cbind(rownames(pvalue),pvalue)
 colnames(pvalue)[1]<-"LncRNA"
 new_cor<-melt(as.data.frame(cor),id.vars="LncRNA",variable.name ="mRNA",value.name ="correlation")
 new_p<-melt(as.data.frame(pvalue),id.vars="LncRNA",variable.name ="mRNA",value.name ="p.value")
 
 
 all_cor<-merge(new_cor, new_p, by.x=c("LncRNA","mRNA"), by.y=c("LncRNA","mRNA"))
 all_cor<-na.omit(all_cor)
 all_cor$p.value<-as.numeric(all_cor$p.value)
 all_cor$correlation<-as.numeric(all_cor$correlation)
 
 
 write.table(all_cor,paste("trans_DElncRNAs_coexpression_DEmRNA_result.txt",sep=""),quote=FALSE,row.names=FALSE)
 
 all_fli<-subset(all_cor,abs(all_cor$correlation)>0.99&all_cor$p.value <0.001)
 write.table(all_fli,paste("filt_0.99_trans_DElncRNAs_coexpression_DEmRNA_result.txt",sep=""),quote=FALSE,row.names=FALSE)
 ```

```{bash eval =F}
#Cytoscape作图 互作网络（mRNA-lncRNA）
##思路：找到与cis调控 tris调控 交集lncRNA 与 mRNA，然后选择与该lncRNA相关
```

