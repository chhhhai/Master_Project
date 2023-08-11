# 基于D-loop区域高通量测序分析——pipeline

```{bash eval = F}
#1.拆分barcode软件的选择
##由于使用stack软件拆分数据后为1kb，感觉可能不适用于这种只是拆分单一区域的序列（参考扩增子测序拆分）
###选择使用 fqkit seqtk_demultiplex barcodeSpliter
####最后选择barcodeSpliter（https://github.com/atlasbioinfo/barcodeSpliter）支持双端barcode拆分
pip install barcodeSpliter

barcodeSpliter HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1_1.clean.fq.gz HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1_2.clean.fq.gz ../barcode.tsv

for (( i=2;i <= 12;i++ )) 
do
echo cd HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}
done

for (( i=2;i <= 12;i++ )) 
do
echo cd /home/xiajh2/data/ach/HTS-D-loop/HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}
echo barcodeSpliter HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}_1.clean.fq.gz HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}_2.clean.fq.gz ../barcode.tsv
done

for (( i=2;i <= 12;i++ )) 
do
cd /home/xiajh2/data/ach/HTS-D-loop/HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}
barcodeSpliter HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}_1.clean.fq.gz HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}_2.clean.fq.gz ../barcode.tsv
done


#2.构建索引（选择基于引物开发选择的D-loop区域）
bowtie2-build MT_D-loop.fa MT_D-loop.fa_index

bowtie2 -x  ../MT_D-loop.fa_index -1 HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1_1.clean.fq.gz -2 HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1_2.clean.fq.gz -S HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1.sam -p 6

samtools view -bS -t MT_D-loop.fa_index HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1.sam > HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1.bam

samtools sort HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1.bam -o HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_1.sort.bam

bcftools mpileup -Ou -f MT_D-loop.fa *.sort.bam |     bcftools call -Ou -mv | bcftools filter -s LowQual -e '%QUAL<10 || DP>2000' > ach_TEST_var_flt_qual10_dp30000.vcf


bowtie2 -x  ../MT_D-loop.fa_index -1 Barcode01.R1.fastq.gz -2  Barcode01.R2.fastq.gz -S Barcode01.sam -p 6

samtools view -bS -t ../MT_D-loop.fa_index Barcode01.sam > Barcode01.bam

samtools sort Barcode01.bam -o Barcode01.sort.bam

bowtie2 -x  ../MT_D-loop.fa_index -1 Barcode02.R1.fastq.gz -2  Barcode02.R2.fastq.gz -S Barcode02.sam -p 6 > Barcode02.log 2>&1

samtools view -bS -t ../MT_D-loop.fa_index Barcode02.sam > Barcode02.bam

samtools sort Barcode02.bam -o Barcode02.sort.bam
```



```{bash eval =F}
for i in 1 
do
echo cd HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}
i=i+1
done

bowtie2 -x /home/xjh3/data/ach/guiyu_gwas/ref/genome.fa_index -1 ${i}.R1.clean.fq.gz -2 ${i}.R2.clean.fq.gz -S ${i}.sam -p 12 >${i}.log 2>&1
samtools view -bS -t /home/xjh3/data/ach/guiyu_gwas/ref/genome.fa_index ${i}.sam > ${i}.bam
samtools sort ${i}.bam -o ${i}.sort.bam
rm ${i}.sam
rm ${i}.bam
done
```



```{bash eval =F}


cd-hit -i s1.sort.fa -o s1.sort.db100.fa -c 1 -n 5 -M 4000 -d 0 -T 4

flash Barcode01.R1.fastq.gz Barcode01.R2.fastq.gz -p 33 -r 250 -s 100 -o flash_out

seqtk seq -A flash_out.extendedFrags.fastq > flash_out.extendedFrags.fasta

cd-hit -i flash_out.extendedFrags.fasta -o s1.sort.db100.fa -c 1 -n 5 -M 4000 -d 0 -T 4

seqkit rmdup -s flash_out.extendedFrags.fasta -o flash_out.extendedFrags.rmdup.fasta


####https://blog.csdn.net/woodcorpse/article/details/106554033

cutadapt -g GCATGTTAGATATATAATGTAATTGTAATGGTTTG -G ACGGTTCTGGAATAAGTGCTCGGCA -o Barcode01_trim.R1.fastq.gz -p Barcode01_trim.R2.fastq.gz Barcode01.R1.fastq.gz Barcode01.R2.fastq.gz

flash Barcode01_trim.R1.fastq.gz Barcode01_trim.R2.fastq.gz -p 33 -r 200 -x 0  -o flash_trim_out

seqtk seq -A flash_trim_out.extendedFrags.fastq > flash_trim_out.extendedFrags.fasta

seqkit rmdup -s flash_trim_out.extendedFrags.fasta -o flash_trim_out.extendedFrags.rmdup.fasta

cd-hit -i flash_trim_out.extendedFrags.rmdup.fasta -o s1.sort_gai.db100.fa -c 1 -n 5 -M 4000 -d 0.9 -T 4

/home/xjh3/data/ach/20230110_HTS_D-loop/gDNA/HTS-Dloop-gDNA-1_HTS-Dloop-gDNA-1_1/Output/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit --input-bed gwas.bed gwas.bim gwas.fam \
        --output-max gwas.phased.haps gwas.phased.sample
        
        
        
/home/xjh3/data/ach/20230110_HTS_D-loop/gDNA/HTS-Dloop-gDNA-1_HTS-Dloop-gDNA-1_1/Output/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit \
--input-vcf /home/xjh3/data/ach/20230110_HTS_D-loop/gDNA/HTS-Dloop-gDNA-1_HTS-Dloop-gDNA-1_1/Output/ach_gDNA_D-loop_var_flt_qual10_dp30000.vcf \
-O gwas.phased


```



```{bash eval = F}
#测试 pipeline
for i in `tail -n+2 doc/design.txt | cut -f 1`;do
  vsearch --fastq_mergepairs Barcode01_trim.R1.fastq.gz --reverse Barcode01_trim.R2.fastq.gz \
  --fastqout Barcode01_trim.merged.fq --relabel Barcode01_trim.
  
  
vsearch --fastx_filter Barcode01_trim.merged.fq \
  --fastq_maxee_rate 0.01 \
  --fastaout filtered.fa 
  
  vsearch --derep_fulllength filtered.fa \
  --sizeout --minuniquesize 8 \
  --output uniques.fa
  
 vsearch --cluster_fast uniques.fa \
  --id 1 --centroids asvs.fa \
  --relabel ASV_
  
vsearch --usearch_global filtered.fa \
  --db asvs.fa \
  --id 1 \
  --otutabout asvtab.txt
  
  vsearch --cluster_unoise uniques.fa  --centroids asvs_unoise3.fa --relabel ASV_ #OTU序列重命名
  
vsearch --usearch_global filtered.fa \
  --db asvs_unoise3.fa \
  --id 1 \
  --otutabout result/asv_unoise_otutab.txt
  

##pipeline v1.0
cutadapt -g GCATGTTAGATATATAATGTAATTGTAATGGTTTG -G ACGGTTCTGGAATAAGTGCTCGGCA -o Barcode01_trim.R1.fastq.gz -p Barcode01_trim.R2.fastq.gz Barcode01.R1.fastq.gz Barcode01.R2.fastq.gz

vsearch --fastq_mergepairs Barcode01_trim.R1.fastq.gz --reverse Barcode01_trim.R2.fastq.gz \
  --fastqout Barcode01_trim.merged.fq --relabel Barcode01_trim.
 
vsearch --fastx_filter Barcode01_trim.merged.fq \
  --fastq_maxee_rate 0.01 \
  --fastaout filtered.fa
  
vsearch --derep_fulllength filtered.fa \
  --sizeout --minuniquesize 8 \
  --output uniques.fa
  
vsearch --cluster_unoise uniques.fa  --centroids asvs_unoise3.fa --relabel ASV_ #OTU序列重命名
  
vsearch --usearch_global filtered.fa \
  --db asvs_unoise3.fa \
  --id 1 \
  --otutabout asv_unoise_otutab.txt

cutadapt -g AACCATTAGATATATAATGTAATTGTAATGGTTTG -G ATACGTCTGGAATAAGTGCTCGGCA -o Barcode02_trim.R1.fastq.gz -p Barcode02_trim.R2.fastq.gz Barcode02.R1.fastq.gz Barcode02.R2.fastq.gz


vsearch --fastq_mergepairs Barcode02_trim.R1.fastq.gz --reverse Barcode02_trim.R2.fastq.gz \
  --fastqout Barcode02_trim.merged.fq --relabel Barcode02_trim.
 
vsearch --fastx_filter Barcode02_trim.merged.fq \
  --fastq_maxee_rate 0.01 \
  --fastaout Barcode02_filtered.fa
  
vsearch --derep_fulllength Barcode02_filtered.fa \
  --sizeout --minuniquesize 8 \
  --output Barcode02_uniques.fa
  
vsearch --cluster_unoise Barcode02_uniques.fa  --centroids Barcode02_asvs_unoise3.fa --relabel ASV_ #OTU序列重命名
  
vsearch --usearch_global Barcode02_filtered.fa \
  --db Barcode02_asvs_unoise3.fa \
  --id 1 \
  --otutabout Barcode02_asv_unoise_otutab.txt
  
  
  
for (( i=1;i <= 12;i++ )) 
do
cd /home/xiajh2/data/ach/HTS-D-loop/HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-1_${i}/Output/
cp *.fastq.gz /home/xiajh2/data/ach/HTS-D-loop/02.all_fastq.gz
done

for (( i=1;i <= 12;i++ )) 
do
cd /home/xiajh2/data/ach/HTS-D-loop/HTS-Dloop-eDNA-1_HTS-Dloop-eDNA-2_${i}/Output/
cp *.fastq.gz /home/xiajh2/data/ach/HTS-D-loop/02.all_fastq.gz
done
for i in 
vsearch --fastq_mergepairs Barcode01_trim.R1.fastq.gz --reverse Barcode01_trim.R2.fastq.gz \
  --fastqout Barcode01_trim.merged.fq --relabel Barcode01_trim.
  
  
  




for i in *.R1.fastq.gz
do
i=${i%.R1.fastq.gz*}
cutadapt -g GCATGTTAGATATATAATGTAATTGTAATGGTTTG -G ACGGTTCTGGAATAAGTGCTCGGCA -o ${i}_trim.R1.fastq.gz -p ${i}_trim.R2.fastq.gz ${i}.R1.fastq.gz ${i}.R2.fastq.gz
done

for i in *.R1.fastq.gz
do
i=${i%.R1.fastq.gz*}
cutadapt -g AACCATTAGATATATAATGTAATTGTAATGGTTTG -G ATACGTCTGGAATAAGTGCTCGGCA -o ${i}_trim.R1.fastq.gz -p ${i}_trim.R2.fastq.gz ${i}.R1.fastq.gz ${i}.R2.fastq.gz
done

for i in *_trim.R1.fastq.gz
do
i=${i%_trim.R1.fastq.gz*}
vsearch --fastq_mergepairs ${i}_trim.R1.fastq.gz --reverse ${i}_trim.R2.fastq.gz --fastqout ${i}.merged.fq --relabel ${i}.
done

mkdir temp

cat *.merged.fq > temp/all.fq

vsearch --fastx_filter temp/all.fq \
  --fastq_maxee_rate 0.01 \
  --fastaout temp/filtered.fa

cd temp/

vsearch --derep_fulllength filtered.fa \
  --sizeout --minuniquesize 8 \
  --output uniques.fa
  
vsearch --cluster_unoise uniques.fa  --centroids asvs_unoise3.fa --relabel ASV_ #OTU序列重命名
  
vsearch --usearch_global filtered.fa \
  --db asvs_unoise3.fa \
  --id 1 \
  --otutabout asv_unoise_otutab.txt
```



```{bash eval =F}
vcftools --vcf ach_gDNA_D-loop_var_flt_qual10_dp30000.vcf --out vcf_pi --site-pi

vcftools --vcf ach_gDNA_D-loop_var_flt_qual10_dp30000.vcf 
	--keep sample.list 
	--window-pi 200 
	--out gDNA_pi





```

