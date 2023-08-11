## 2022.10_鱼类多样性

```{bash eval = F}
qiime tools import \
   --type SampleData[PairedEndSequencesWithQuality] \
   --input-path ../new_12s_manifest \
   --output-path reads.qza \
   --input-format PairedEndFastqManifestPhred33V2
qiime demux summarize --i-data reads.qza --o-visualization read_summary.qzv   

time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /home/xjh3/data/ach/2022_10_fish_diversity/2022_11_08_analysis/reads.qza \
  --p-n-threads 12 \
  --p-trim-left-f 21 --p-trim-left-r 27 \
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

library(pacman)
pacman::p_load(tidyverse,magrittr,stringr)
otu <- "otu_table.tsv" %>%
  read.delim(check.names = FALSE,header = T,sep="\t")

rown <- paste0("ASV",seq_len(nrow(otu)))
otu[,1] <- rown
colnames(otu)[1] <- paste0("asv",colnames(data)[1])
write.table (otu,file ="otu_table_1.tsv", sep ="\t", row.names = F) 

sed 's/"//g' otu_table_1.tsv >otu_table_2.tsv
sed -i "s/\s*$//g" otu_table_2.tsv
biom convert -i otu_table_2.tsv -o feature-table.biom --to-hdf5 --table-type="OTU table"

otu_table_3.tsv

biom convert -i otu_table_3.tsv -o feature-table.biom --to-hdf5 --table-type="OTU table"
qiime tools import \
  --input-path feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path otu_table.qza
  
cd ..  
qiime tools export \
--input-path rep-seqs.qza \
--output-path phyloseq
cd phyloseq
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
```

