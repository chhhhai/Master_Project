# 鱼类多样性pipeline

```bash eval=F
#训练数据库 数据库来源：https://zenodo.org/record/6522134#.YoHoRXVBybg
qiime feature-classifier extract-reads \
  --i-sequences 12S-seqs-derep-uniq.qza \
  --p-f-primer GTCGGTAAAACTCGTGCCAGC \
  --p-r-primer CATAGTGGGGTATCTAATCCCAGTTTG \
  --o-reads ref-seqs.qza
  
time qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy 12S-tax-derep-uniq.qza \
  --o-classifier Mifish_version_12S.qza

mkdir qiime_analyze
cd qiime_analyze
qiime tools import \
   --type SampleData[PairedEndSequencesWithQuality] \
   --input-path ../12S_manifest \
   --output-path reads.qza \
   --input-format PairedEndFastqManifestPhred33V2
qiime demux summarize --i-data reads.qza --o-visualization read_summary.qzv   


mkdir /home/xjh3/data/ach/fish_diversity/figaro_data/

#github下载figaro https://www.jianshu.com/p/66ebd1d4558a
cd /home/xjh3/data/ach/fish_diversity/figaro_data

gunzip *.gz

for name in `ls *_1.clean.fq`;do echo $name ${name%_1.clean.fq}_R1.fastq;done

for name in `ls *_1.clean.fq`;do mv $name ${name%_1.clean.fq}_R1.fastq;done

for name in `ls *_2.clean.fq`;do mv $name ${name%_2.clean.fq}_R2.fastq;done

for name in `ls *_R1.fastq`;do echo $name ${name%-*}_R1.fastq;done



python /home/xjh3/data/ach/fish_diversity/figaro/figaro/figaro/figaro.py -i figaro_data -o figaro_data -f 21 -r 27 -a 330 -F zymo

time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs reads.qza \
  --p-n-threads 12 \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-trunc-len-f 0 --p-trunc-len-r 0 \
  --o-table dada2-table.qza \
  --o-representative-sequences dada2-rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
  
cp dada2-table.qza table.qza
cp dada2-rep-seqs.qza rep-seqs.qza

mkdir phyloseq
qiime tools export \
--input-path table.qza \
--output-path phyloseq

biom convert \
-i phyloseq/feature-table.biom \
-o phyloseq/otu_table.tsv \
--to-tsv

cd phyloseq
sed -i '1d' otu_table.tsv
sed -i 's/#OTU ID/ASV/' otu_table.tsv
cd ../

##########################R语言操作###############################################
library(pacman)
pacman::p_load(tidyverse,magrittr,stringr)
otu <- "otu_table.tsv" %>%
  read.delim(check.names = FALSE,header = T,sep="\t")

rown <- paste0("ASV",seq_len(nrow(otu)))
otu[,1] <- rown
colnames(otu)[1] <- paste0("asv",colnames(data)[1])
write.table (otu,file ="otu_table.tsv", sep ="\t", row.names = F) 
################################################################################
sed 's/"//g' otu_table.tsv > otu_table_1.tsv

mkdir temp
cd temp
biom convert -i ../otu_table_1.tsv -o feature-table.biom --to-hdf5 --table-type="OTU table"
qiime tools import \
  --input-path feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path otu_table.qza
cd ..
qiime tools export \
--input-path /home/xjh3/data/ach/fish_diversity/qiime_analyze/rep-seqs.qza \
--output-path phyloseq
cd phyloseq


less dna-sequences.fasta |paste - -|sed '1i ASVID,seq' > rep.fa
###################################R############################################
library(pacman)
pacman::p_load(tidyverse,magrittr,stringr)
rep <- "rep.fa" %>%
  read.delim(check.names = FALSE, row.names = 1) %>%
  set_rownames(paste0(">ASV", seq_len(nrow(.))))
write.table (rep,file ="rep.xls", sep ="\t", row.names = T)  
################################################################################
less rep.xls|sed '1d'|sed 's/"//g'|\
sed 's/\r//g'|tr "\t" "\n" > rep-seqs.fasta

time qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path rep-seqs.fasta \
--output-path rep-seqs.qza

mv ./temp/otu_table.qza ./
time qiime feature-table summarize \
--i-table /home/xjh3/data/ach/fish_diversity/qiime_analyze/table.qza \
--o-visualization table.qzv \

time qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

qiime feature-classifier classify-sklearn \
  --i-classifier /home/xjh3/data/ach/fish_diversity/database/Mifish_version_12S.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
  
qiime tools export \
--input-path taxonomy.qza \
--output-path phyloseq
```



```{bash eval=F}
#换了五月份的数据库重做2022_5（采用这个，上面建议不用看了）
time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /home/xjh3/data/ach/fish_diversity/qiime_analyze/reads.qza \
  --p-n-threads 12 \
  --p-trim-left-f 21 --p-trim-left-r 27 \
  --p-trunc-len-f 0 --p-trunc-len-r 0 \
  --o-table dada2-table.qza \
  --o-representative-sequences dada2-rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
  
qiime feature-classifier extract-reads \
  --i-sequences 12S-seqs-derep-uniq.qza \
  --p-f-primer GTCGGTAAAACTCGTGCCAGC \
  --p-r-primer CATAGTGGGGTATCTAATCCCAGTTTG \
  --o-reads ref-seqs.qza
  
time qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy 12S-tax-derep-uniq.qza \
  --o-classifier Mifish_version_12S.qza
  
  
cp dada2-table.qza table.qza
cp dada2-rep-seqs.qza rep-seqs.qza

mkdir phyloseq
qiime tools export \
--input-path table.qza \
--output-path phyloseq

biom convert \
-i phyloseq/feature-table.biom \
-o phyloseq/otu_table.tsv \
--to-tsv
cd phyloseq; sed -i '1d' otu_table.tsv
sed -i 's/#OTU ID/ASV/' otu_table.tsv
cd ../

library(pacman)
pacman::p_load(tidyverse,magrittr,stringr)
otu <- "otu_table.tsv" %>%
  read.delim(check.names = FALSE,header = T,sep="\t")

rown <- paste0("ASV",seq_len(nrow(otu)))
otu[,1] <- rown
colnames(otu)[1] <- paste0("asv",colnames(data)[1])
write.table (otu,file ="otu_table_1.tsv", sep ="\t", row.names = F)   

biom convert -i otu_table_1.tsv -o feature-table.biom --to-hdf5 --table-type="OTU table"
qiime tools import \
  --input-path feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path otu_table.qza
  
qiime tools export \
--input-path rep-seqs.qza \
--output-path phyloseq

less dna-sequences.fasta |paste - -|sed '1i ASVID,seq' > rep.fa

library(pacman)
pacman::p_load(tidyverse,magrittr,stringr)
rep <- "rep.fa" %>%
  read.delim(check.names = FALSE, row.names = 1) %>%
  set_rownames(paste0(">ASV", seq_len(nrow(.))))
write.table (rep,file ="rep.xls", sep ="\t", row.names = T)  

less rep.xls|sed '1d'|sed 's/"//g'|\
sed 's/\r//g'|tr "\t" "\n" > rep-seqs.fasta

time qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path rep-seqs.fasta \
--output-path rep-seqs.qza

time qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

qiime feature-classifier classify-sklearn \
  --i-classifier /home/xjh3/data/ach/fish_diversity/2022_5_16_analyze/database_2022_5/Mifish_version_12S.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
  
#剩下的R作图——————https://www.jianshu.com/p/d38b2f5aec43
导入数据出了问题发现是otu_table_1.tsv格式不正确，需要将引号删除!!!!!!

sed 's/"//g' otu_table_1.tsv >otu_table_2.tsv

biom convert -i otu_table_2.tsv -o feature-table.biom --to-hdf5 --table-type="OTU table"
qiime tools import \
  --input-path feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path otu_table.qza
  
重新生成otu_table.qza替换
开始R作图



#由于顺序不是自己想要的 使用了excel排列了列，也可以使用R tidyverse select函数 但是懒得搞了，excel 快，懒得思考了 ，sort完 文件为otu_table_3.tsv
#将三个生物学重复合并成一组，只看组间相对丰度差异，文件为otu_table_4.tsv（excel加合）
mkdir /home/xjh3/data/ach/fish_diversity/2022_5_16_analyze/phyloseq/sort_and_only-group

sed -i "s/\s*$//g" otu_table_3.tsv

sed -i "s/\s*$//g" otu_table_4.tsv

biom convert -i otu_table_3.tsv -o feature-table_3.biom --to-hdf5 --table-type="OTU table"
qiime tools import \
  --input-path feature-table_3.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path otu_table_3.qza
  
biom convert -i otu_table_4.tsv -o feature-table_4.biom --to-hdf5 --table-type="OTU table"
qiime tools import \
  --input-path feature-table_4.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path otu_table_4.qza  
```

