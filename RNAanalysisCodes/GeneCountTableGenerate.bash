### 统计每一个样本中mapping到常染色和性染色的reads数量。


## 1）将gtf文件中的基因位置信息提取出来。
awk '$3 == "gene"' /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.98.gtf | awk '{print $10, "\t", $1, "\t", $4, "\t", $5}' | sed 's/[",;]//g' | sort -k1,1n> /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigGenomeGenePosition.txt


# 统计每个样本比对到X染色体上的reads的总和
cd /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/
for sample in `ll *.gene.count| awk '{print $NF}'`
do
	join -o '0,1.2,2.2' "$sample" <(sort -k1,1n /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigGenomeGenePosition.txt) | sed 's/ /\t/g' | awk '$3 == "X"' |awk '{sum += $2} END {print sum}'
done > Xcount.txt


# 统计每个样本比对到常染色体上的reads的总和
cd /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/
for sample in `ll *.gene.count| awk '{print $NF}'`
do
	join -o '0,1.2,2.2,2.3,2.4' "$sample" <(sort -k1,1 /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigGenomeGenePosition.txt) | sed 's/ /\t/g' | awk '$3 != "X"' |awk '{sum += $2} END {print sum}'
done > AutosomeCount.txt


# 将样本名、比对到常染色体reads总和、比对到X染色体reads总和，一起构建一个表格。
paste <(ll *.gene.count| awk '{print $NF}') AutosomeCount.txt Xcount.txt > AllAutoXCount_table.txt


#####################################################################################（不用）
## 提取比对到母本基因组的样本数据（这里通过改变分隔符，来实现选定 样本名的第二个字符）
awk '$1 ~ "Ref[D,L]S"' AllAutoXCount_table.txt | awk -F "" '{if ($2 == "B") print $0, "Maternal", "Boar"; else if ($2 == "S") print $0, "Maternal", "Sow" } ' | awk -F '_Ref' '{print $1, $2}' | awk '{print $1, $3, $4, $5, $6}' > AllAutoXCountFemale_table.txt
## 提取比对到父本基因组的样本数据
awk '$1 ~ "Ref[D,L]B"' AllAutoXCount_table.txt | awk -F "" '{if ($2 == "B") print $0, "Paternal", "Boar"; else if ($2 == "S") print $0, "Paternal", "Sow" } ' | awk -F '_Ref' '{print $1, $2}' | awk '{print $1, $3, $4, $5, $6}' > AllAutoXCountMale_table.txt
#####################################################################################


###################################################################################
#### 统计mapped到X染色体基因上的reads数量以及位置(利用基因名进行搜索)。
cd /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment_insert/GeneCount
for sample in *.gene.count
do
join -o '0,1.2,2.2,2.3,2.4' "$sample" <(awk 'BEGIN{OFS="\t"}{print $4,$2,$3,$1}' /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigXChrGenePosition.txt|sort -k1,1) | sed 's/ /\t/g' | sort -k3,3n|cut -f1,2 > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment_insert/GeneCount/XChrGene/`echo $sample|cut -d "." -f1`.XchrGene.position.txt
done
#### 合并所有样本中比对到X染色体基因上的reads数
module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604
ext='XchrGene.position.txt' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"AllSamplesXChrGeneCountTable_Intersect.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.XchrGene.position.txt//g' AllSamplesXChrGeneCountTable_Intersect.txt|sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > AllSamplesXChrGeneCountTable_Intersect_final.txt


#### 移除公猪中的PAR区域（包括在X染色体与Y染色体中）
for sample in `ls ?B*Ref*gene.count | sed 's/.gene.count//g'`
do
join -o '0,1.2,2.2,2.3,2.4' "$sample".gene.count <(awk 'BEGIN{OFS="\t"}{print $4,$2,$3,$1}' /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigXChrGenePosition.txt|sort -k1,1) | sed 's/ /\t/g' | awk '!($5=="X" && $4<=6448126)' |cut -f1,2 > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment_insert/GeneCount/XChrGene/ParRemoved/"$sample".MaleParRemoved.gene.count.txt
done

#### 合并所有样本中比对到X染色体基因上的reads数
module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604
ext='.txt' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"AllSamplesXChrGeneCountTable_IntersectParRemoved.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.MaleParRemoved.gene.count.txt//g' AllSamplesXChrGeneCountTable_IntersectParRemoved.txt|sed 's/.XchrGene.position.txt//g'|sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > AllSamplesXChrGeneCountTable_Intersect_ParRemoved_final.txt
####################################################################################



#### 移除公猪的X染色体与Y染色体上的基因
for sample in `ls ?B*Ref*gene.count | sed 's/.gene.count//g'`
do
join -o '0,1.2,2.2,2.3,2.4' "$sample".gene.count <(sort -k1,1 /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigGenomeGenePosition.txt) | sed 's/ /\t/g' | awk '!($3=="X" || $3=="Y")' |cut -f1,2 > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/XYremoved/"$sample".XYRemoved.gene.count
done


#### 移除母猪的X染色体上的基因
for sample in `ls ?S*Ref*gene.count | sed 's/.gene.count//g'`
do
join -o '0,1.2,2.2,2.3,2.4' "$sample".gene.count <(sort -k1,1 /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigGenomeGenePosition.txt) | sed 's/ /\t/g' | awk '!($3=="X")' |cut -f1,2 > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/XXremoved/"$sample".XXRemoved.gene.count
done


#### 合并所有样本，所有染色体（包括常染色体）的gene.count文件（R里面操作）
# Load R v3.5.1
module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604

ext='gene.count' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"AllSamplesAllChrGeneCountTable.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.gene.count//g' AllSamplesAllChrGeneCountTable.txt | sed 's/.parRemoved//g'|sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > AllSamplesAllChrGeneCountTable_final.txt


###合并所有猪的常染色体
cp /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/XYRemoved/*XYRemoved.gene.count /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/XXRemoved/
cd XXRemoved

module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604

ext='.count' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){ read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"AllSamplesAutosomeGeneCountTable.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.gene.count//g' AllSamplesAutosomeGeneCountTable.txt |sed 's/.XYRemoved//g'|sed 's/.XXRemoved//g'|sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > AllSamplesAutosomeGeneCountTable_final.txt


####合并所有公猪
module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604

ext='.count' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){ read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"BoarSamplesXYremovedGeneCountTable.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.gene.count//g' BoarSamplesXYremovedGeneCountTable.txt | sed 's/XYRemoved//g'|sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > BoarSamplesXYremovedGeneCountTable_final.txt


### 合并所有母猪
module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604

ext='.count' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){ read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"SowSamplesGeneCountTable.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.gene.count//g' SowSamplesGeneCountTable.txt |sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > SowSamplesGeneCountTable_final.txt









##########################################################################################
#求同一时期的两个组合交集
cd /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified
sort <(cut -f6 DB1_informative.bed) <(cut -f6 LB1_informative.bed) |uniq -d > Stage1.txt
sort <(cut -f6 DB2_informative.bed) <(cut -f6 LB2_informative.bed) |uniq -d > Stage2.txt
sort <(cut -f6 DB3_informative.bed) <(cut -f6 LB3_informative.bed) |uniq -d > Stage3.txt
sort <(cut -f6 DB4_informative.bed) <(cut -f6 LB4_informative.bed) |uniq -d > Stage4.txt

#提取出同一时期不同组合的亲本的informative variants																					
cd /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified
for sample in DB1 DB2 DB3 DB4 DS1 DS2 DS3 DS4 LB1 LB2 LB3 LB4 LS1 LS2 LS3 LS4
do
awk '{print $6"\t"$0}' "$sample"_informative.bed | sort -k1,1 | join Stage`echo "$sample" |awk 'BEGIN{FS="_"}{$1=substr($1,3);print $1}'`.txt - | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > "$sample"_informative_insert.bed
done


#############################################################
##### 求每个时期不同亲本informative variants 覆盖的 exon 个数
cd /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified
for sample in `ls *informative_insert.bed`
do
sort -k1,1 "$sample" | /mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools intersect -a - -b <(awk '$3 == "exon"' /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.98.gtf) -wa | uniq | cut -f 6 | sort | uniq | wc -l
done

#### 基因个数
for sample in `ls *informative_insert.bed`
do
sort -k1,1 "$sample" |/mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools intersect -a - -b <(awk '$3 == "exon"' /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.98.gtf) -wa -wb | cut -f 6,16 | awk -F "\"" '{print $1"\t"$2}' | awk '{print $3"\t"$1}'| sort -k1,1 | /mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools groupby -g 1 -c 2 -o count_distinct |wc -l
done
##############################################################


#### 5.2 intersect BAM with the informative BED to find the informative reads
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
bam_path=$(echo $line |cut -d " " -f 2)
bed_path=$(echo $line |cut -d " " -f 3)
sbatch --export=sample_name="$sample_name",bam_path="$bam_path",bed_path="$bed_path" --output=log/ModifiedGenomeAnalysis_V2/Intersect/"$sample_name".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/Intersect/"$sample_name".sbatch.err /mnt/home/quanjian/quan/RNA/script/BedtoolsIntersect_V3.batch
done < IntersectSample_bam_V3.txt


#### 5.3.1 process the overlap to assign reads (also can process AssignRead.batch)
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
maleRef_var_read=$(echo $line |cut -d " " -f 2)
femaleRef_var_read=$(echo $line |cut -d " " -f 3)
join -a 1 -a 2 -t $'\t' -o '0,1.2,1.3,2.2,2.3' -e '-' "$maleRef_var_read" "$femaleRef_var_read" | perl /mnt/home/quanjian/quan/WGS/scripts/aseUtils/assignRead.pl > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned_insert/"$sample_name"read.genome.assign.txt
done < Variant_reads_V2.txt


#### 5.3.2 get the reads that are assigned to each genome
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
ref_1=$(echo $line |cut -d " " -f 2|cut -d "/" -f 10 | cut -d "." -f 1| grep -o ".\{3\}$")
ref_2=$(echo $line |cut -d " " -f 3|cut -d "/" -f 10 | cut -d "." -f 1| grep -o ".\{3\}$")
awk '$2 == 1 {print $1}' /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned_insert/"$sample_name"read.genome.assign.txt | awk -F "/" '{print $1}' | uniq | sort | uniq > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned_insert/Separate/"$sample_name"Ref"$ref_1".read
awk '$2 == 2 {print $1}' /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned_insert/"$sample_name"read.genome.assign.txt | awk -F "/" '{print $1}' | uniq | sort | uniq > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned_insert/Separate/"$sample_name"Ref"$ref_2".read
done < Variant_reads_V2.txt



#### 5.4 prepare constitutive exons (non-overlapping) for each genome


### 5.5.1 get read length
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
bam_path=$(echo $line |cut -d " " -f 2)
bed_path=$(echo $line |cut -d " " -f 3)
read_path=$(echo $line |cut -d " " -f 4)
sbatch --export=sample_name="$sample_name",bam_path="$bam_path",bed_path="$bed_path",read_path="$read_path" /mnt/home/quanjian/quan/RNA/script/GetReadLength_V2.batch
done < GetReadLength_V2.txt


### 5.5.2 get the gene assignment for the first read
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
bam_path=$(echo $line |cut -d " " -f 2)
bed_path=$(echo $line |cut -d " " -f 3)
overlap_bed_path=$(echo $line |cut -d " " -f 4)
read_len_path=$(echo $line |cut -d " " -f 5)
sbatch --export=sample_name="$sample_name",bam_path="$bam_path",bed_path="$bed_path",overlap_bed_path="$overlap_bed_path",read_len_path="$read_len_path" /mnt/home/quanjian/quan/RNA/script/GeneAssignFirstRead_V2.batch
done < GeneAssignment_V2.txt


### 5.5.3 get the gene assignment for the second read
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
bam_path=$(echo $line |cut -d " " -f 2)
bed_path=$(echo $line |cut -d " " -f 3)
overlap_bed_path=$(echo $line |cut -d " " -f 4)
read_len_path=$(echo $line |cut -d " " -f 5)
sbatch --export=sample_name="$sample_name",bam_path="$bam_path",bed_path="$bed_path",overlap_bed_path="$overlap_bed_path",read_len_path="$read_len_path" /mnt/home/quanjian/quan/RNA/script/GeneAssignSecondRead_V2.batch
done < GeneAssignment_V2.txt


### 5.5.4 combine both reads, only consider one fragment if both reads map to the same gene
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
first_path=$(echo $line |cut -d " " -f 2)
second_path=$(echo $line |cut -d " " -f 3)
sbatch --export=sample_name="$sample_name",first_path="$first_path",second_path="$second_path" /mnt/home/quanjian/quan/RNA/script/GeneCount_V2.batch
done < GeneCount_V2.txt



#### 移除公猪的X染色体与Y染色体上的基因
for sample in `ls ?B*Ref*gene.count | sed 's/.gene.count//g'`
do
join -o '0,1.2,2.2,2.3,2.4' "$sample".gene.count <(sort -k1,1 /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigGenomeGenePosition.txt) | sed 's/ /\t/g' | awk '!($3=="X" || $3=="Y")' |cut -f1,2 > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment_insert/GeneCount/XYremoved/"$sample".XYRemoved.gene.count
done

#### 移除母猪的X染色体上的基因
for sample in `ls ?S*Ref*gene.count | sed 's/.gene.count//g'`
do
join -o '0,1.2,2.2,2.3,2.4' "$sample".gene.count <(sort -k1,1 /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment/GeneCount/StageTissueRef/PigGenomeGenePosition.txt) | sed 's/ /\t/g' | awk '!($3=="X")' |cut -f1,2 > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment_insert/GeneCount/XXremoved/"$sample".XXRemoved.gene.count
done


#### 合并所有样本的gene.count文件（R里面操作）
# Load R v3.5.1
####合并所有公猪
module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604

ext='.count' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){ read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"BoarSamplesXYremovedGeneCountTableInsert.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.gene.count//g' BoarSamplesXYremovedGeneCountTableInsert.txt | sed 's/XYRemoved//g'|sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > BoarSamplesXYremovedGeneCountTableInsert_final.txt


### 合并所有母猪
module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604

ext='.count' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){ read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"SowSamplesGeneCountTableInsert.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.gene.count//g' SowSamplesGeneCountTableInsert.txt |sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > SowSamplesGeneCountTableInsert_final.txt



###合并所有猪的常染色体
cp /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment_insert/GeneCount/XYRemoved/*XYRemoved.gene.count /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/GeneAssignment_insert/GeneCount/XXRemoved/
cd XXRemoved

module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604
R

ext='.count' #can alter this to desired extension
files <- list.files(pattern=ext) #get name of files in a directory
listOfFiles <- lapply(files, function(x){ read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"AllSamplesAutosomeGeneCountTableInsert.txt", quote = F, row.names = T)
##在shell中对样本名进行修改
sed 's/.gene.count//g' AllSamplesAutosomeGeneCountTableInsert.txt |sed 's/.XYRemoved//g'|sed 's/.XXRemoved//g'|sed 's/ /\t/g'|awk '{if(FNR == 1){print "\t"$0}else{print $0}}'  > AllSamplesAutosomeGeneCountTableInsert_final.txt



#######用original gtf 位置进行统计
cut -f 6 DB3_informative_insert.bed |awk 'BEGIN{FS="_";OFS="\t"} {print $1,$2,$3}'|sort -k1,1n| /mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools intersect -a - -b <(awk '$3 == "exon"' /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.98.gtf) -wa | uniq |tr '\t' '_'|sort | uniq | wc -l

#######用lift gtf 位置进行统计
sort -k1,1 LB1_informative_insert.bed |/mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools intersect -a - -b <(awk '$3 == "exon"' /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/LB1/LB1_lift.gtf) -wa | uniq | cut -f 6 | sort | uniq | wc -l


####
for sample in DB1 DS1 LB1 LS1 DB2 DS2 LB2 LS2 DB3 DS3 LB3 LS3 DB4 DS4 LB4 LS4; do sort -k1,1 "$sample"_informative_insert.bed |/mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools intersect -a - -b <(awk '$3 == "exon"' /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/"$sample"/"$sample"_lift.gtf) -wa -wb|cut -f 6,16 | awk -F "\"" '{print $1"\t"$2}' | awk '{print $3"\t"$1}'| sort -k1,1 | /mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools groupby -g 1 -c 2 -o count_distinct |wc -l; done
