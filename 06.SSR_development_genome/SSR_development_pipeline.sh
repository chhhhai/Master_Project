######MISA identify SSR #####
Syntax: p3_in.pl <FASTAfile.misa> 
<FASTAfile.misa> is one of the files which was created by MISA containing information about the type and the localization of each individual microsatellite. Since the DNA sequence is also submitted to Primer3, the original <FASTAfile> needs to be accessible as well. 
The created file <FASTAfile.p3in> is the input file for Primer3, which can then be invoked as "primer3_core <FASTAfile.p3in >FASTAfile.p3out". 
  
Syntax: p3_out.pl <FASTAfile.p3out> <FASTAfile.misa> 
This tool parses the Primer3 output file (FASTAfile.p3out) for the calculated primer sequences, their position and melting temperatures as well as the expected PCR product size and merges this information together with the data from the file <FASTAfile.misa> creating a new file (FASTAfile.results). 



####vim检查misa.ini 的GFF状态######
You forgot to say you are using MISA 2.0 and you had GFF: true in your misa.ini file. If you don't want the gff files, edit the misa.ini to GFF: false.

####下载基因组数据######
xiajh2@xiajh2:~/data/Epinephelus_coioides_grouper_data$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/536/245/GCA_900536245.1_Grouper/GCA_900536245.1_Grouper_genomic.fna.gz

xiajh2@xiajh2:~/data/Epinephelus_coioides_grouper_data$ gzip -d GCA_900536245.1_Grouper_genomic.fna.gz

####把/software/misa_ssr_primer_analysis 所有文件 复制到当前目录#####
cp /software/misa_ssr_primer_analysis/* ./

####去除文件名字中的空格#####
xiajh2@xiajh2:~/data/Epinephelus_coioides_grouper_data$
sed 's/ Epi.*//g'  GCA_900536245.1_Grouper_genomic.fna > GCA_900536245.1_Grouper_genomic_3.fa

####misa 软件分析######
xiajh2@xiajh2:/software/misa_ssr_primer_analysis$ 

perl ./misa.pl ./GCA_900536245.1_Grouper_genomic_3.fa

cp GCA_900536245.1_Grouper_genomic_3.fa.misa test4.fa.misa

cp GCA_900536245.1_Grouper_genomic_3.fa test4.fa


####因为如果直接使用p3_in.pl进行转换生成的文件会比较大，所以下面多了几步#提取misa文件中的染色体编号和开始，结束的位置，两边各延伸150bp，生成一个bed文件。####
cat test4.fa.misa |awk 'NR>1 {print $1"\t"$6-150"\t"$7+150}' > test4.fa.misa.bed
#####起始位置小于1会导致bedtools提取错误#####
awk '{if($2>0) print $0}' test4.fa.misa.bed > test4.fa.misa.bed2

#使用bedtools工具提取重复序列####
bedtools getfasta -fi test4.fa -bed test4.fa.misa.bed2 -fo test4.fa.misa.bed.fa

####去除不好的序列末端######
perl ./est_trimmer.pl test4.fa.misa.bed.fa -amb=2,50 -tr5=T,5,50 -tr3=A,5,50 -cut=100,2000

####重新分析SSR#####
perl ./misa.pl ./test4.fa.misa.bed.fa.results

cp test4.fa.misa.bed.fa.results test5.fa

cp test4.fa.misa.bed.fa.results.misa test5.fa.misa

perl ./p3_in.pl ./test5.fa.misa

523779 records created 449010 records created.

#####然后使用primer3(v 1.1.1;版本号很重要，不然需要修改命令格式)进行设计引物#####

primer3_core <test5.fa.p3in > test5.fa.p3out 


perl ./p3_out.pl  test5.fa.p3out  test5.fa.misa

Primer modelling was successful for 345897 sequences.
Primer modelling failed for 15916 sequences.


cp test5.fa.results GCA_900536245.1_Grouper_genomic_3.txt

cp test5.fa GCA_900536245.1_Grouper_genomic_3_new_SSR.fa


###结果过滤###
grep -v "(A)" GCA_900536245.1_Grouper_genomic_3.txt > T1
grep -v "(G)" T1 > T2
grep -v "(C)" T2 > T3
grep -v "(T)" T3 > T4
cp T4 GCA_900536245.1_Grouper_genomic_3_NO_SINGLE_NT_REPEATS.txt

awk 'END{print NR}' GCA_900536245.1_Grouper_genomic_3_NO_SINGLE_NT_REPEATS.txt

174414 261067

awk 'END{print NR}' GCA_900536245.1_Grouper_genomic_3.txt
361814 449011


####调出fasta 序列(序列名字中的：号会影响blast,所以要先取代)##########

awk '{print $1}' GCA_900536245.1_Grouper_genomic_3_NO_SINGLE_NT_REPEATS.txt >t1

sed 's/:/xxx/g' t1 >t2
sed 's/:/xxx/g' GCA_900536245.1_Grouper_genomic_3_new_SSR.fa > GCA_900536245.1_Grouper_genomic_3_new_SSR_2.fa


makeblastdb -in GCA_900536245.1_Grouper_genomic_3_new_SSR_2.fa -dbtype nucl -title GCA_900536245.1_Grouper_genomic_3_new_SSR_2.fa -parse_seqids -hash_index -out  GCA_900536245.1_Grouper_genomic_3_new_SSR_2.fa





blastdbcmd -db GCA_900536245.1_Grouper_genomic_3_new_SSR_2.fa -entry_batch t2 -outfmt "%f" -out  GCA_900536245.1_Grouper_genomic_3_NO_SINGLE_NT_REPEATS_2.fa

sed 's/xxx/:/g' GCA_900536245.1_Grouper_genomic_3_NO_SINGLE_NT_REPEATS_2.fa > GCA_900536245.1_Grouper_genomic_3_NO_SINGLE_NT_REPEATS_3.fa
