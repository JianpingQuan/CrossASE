# 5.1 find informative sites and make bed files
for path in /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified
do
for pairs in DB1_LS1 DB2_LS2 DB3_LS3 DB4_LS4 LB1_DS1 LB2_DS2 LB3_DS3 LB4_DS4
do
paste "$path"/`echo $pairs | cut -d "_" -f1`/`echo $pairs | cut -d "_" -f1`_lift.bed "$path"/`echo $pairs | cut -d "_" -f2`/`echo $pairs | cut -d "_" -f2`_lift.bed| awk '$4 ~ /[ATCGN]/ && $10 ~ /[ATCGN]/ && $4 != $10' | perl -wne 'chomp $_;
  @line = split /\t/, $_;
  if (length($line[3]) == length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1], "_SNP_", $line[3], "_", $line[9], "\n";
  } elsif (length($line[3]) > length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[9]), "_DEL_";
    print substr($line[3], length($line[9])), "_", "-" x (length($line[3]) - length($line[9])), "\n";
  } else {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[3]), "_INS_";
    print "-" x (length($line[9]) - length($line[3])), "_", substr($line[9], length($line[3])), "\n";
	}' > "$path"/`echo $pairs | cut -d "_" -f1`/`echo $pairs | cut -d "_" -f1`_informative.bed
done
done


for path in /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified
do
for pairs in DB1_LS1 DB2_LS2 DB3_LS3 DB4_LS4 LB1_DS1 LB2_DS2 LB3_DS3 LB4_DS4
do
paste "$path"/`echo $pairs | cut -d "_" -f2`/`echo $pairs | cut -d "_" -f2`_lift.bed "$path"/`echo $pairs | cut -d "_" -f1`/`echo $pairs | cut -d "_" -f1`_lift.bed| awk '$4 ~ /[ATCGN]/ && $10 ~ /[ATCGN]/ && $4 != $10' | perl -wne 'chomp $_;
  @line = split /\t/, $_;
  if (length($line[3]) == length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1], "_SNP_", $line[3], "_", $line[9], "\n";
  } elsif (length($line[3]) > length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[9]), "_DEL_";
    print substr($line[3], length($line[9])), "_", "-" x (length($line[3]) - length($line[9])), "\n";
  } else {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[3]), "_INS_";
    print "-" x (length($line[9]) - length($line[3])), "_", substr($line[9], length($line[3])), "\n";
	}' > "$path"/`echo $pairs | cut -d "_" -f2`/`echo $pairs | cut -d "_" -f2`_informative.bed
done
done



# 5.2 intersect BAM with the informative BED to find the informative reads
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
bam_path=$(echo $line |cut -d " " -f 2)
bed_path=$(echo $line |cut -d " " -f 3)

sbatch --export=sample_name="$sample_name",bam_path="$bam_path",bed_path="$bed_path" --output=log/ModifiedGenomeAnalysis_V2/Intersect/"$sample_name".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/Intersect/"$sample_name".sbatch.err /mnt/home/quanjian/quan/RNA/script/BedtoolsIntersect_V2.batch

done < IntersectSample_bam_V2.txt



/mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools intersect -a /mnt/gs18/scratch/users/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/DB115_Br1_RefLB3.align.bam -b /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/LB3/LB3_informative.bed -split | wc -l 



# 5.3
# process the overlap to assign reads (also can process AssignRead.batch)
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
maleRef_var_read=$(echo $line |cut -d " " -f 2)
femaleRef_var_read=$(echo $line |cut -d " " -f 3)
join -a 1 -a 2 -t $'\t' -o '0,1.2,1.3,2.2,2.3' -e '-' "$maleRef_var_read" "$femaleRef_var_read" | perl /mnt/home/quanjian/quan/WGS/scripts/aseUtils/assignRead.pl 2> /mnt/home/quanjian/quan/RNA/log/assignRead/"$sample_name".assRead.err > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned/"$sample_name"read.genome.assign.txt
done < Variant_reads_test.txt

# get the reads that are assigned to each genome (also can process )
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
ref_1=$(echo $line |cut -d " " -f 2|cut -d "/" -f 10 | cut -d "." -f 1| grep -o ".\{3\}$")
ref_2=$(echo $line |cut -d " " -f 3|cut -d "/" -f 10 | cut -d "." -f 1| grep -o ".\{3\}$")

awk '$2 == 1 {print $1}' /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned/"$sample_name"read.genome.assign.txt | awk -F "/" '{print $1}' | uniq | sort | uniq > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned/Separate/"$sample_name"Ref"$ref_1".read
awk '$2 == 2 {print $1}' /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned/"$sample_name"read.genome.assign.txt | awk -F "/" '{print $1}' | uniq | sort | uniq > /mnt/scratch/quanjian/RNA_seq/ModifiedGenomeAlnSamSort/Assigned/Separate/"$sample_name"Ref"$ref_2".read

done < Variant_reads.txt



# 5.4 prepare constitutive exons (non-overlapping) for each genome
# 5.4.1 generate fai index file through samtools
for sample in `ls /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/`
do
cd /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/"$sample"/
/mnt/research/qgg/software/samtools-1.6/samtools faidx "$sample"_lift.fa
done


# 5.4.2 generate the exon bed file of each parent modified gtf file.
for sample in `ls /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/`
do
cd /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/"$sample"/
awk '$3 == "exon"' "$sample"_lift.gtf | perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[8] =~ m/gene_id \"(.*?)\".*transcript_id \"(.*?)\";/) { print $line[0], "\t", $line[3] - 1, "\t", $line[4], "\t", $1, "\t", $2, "\t", $line[6], "\n";}' > "$sample"_lift.exon.bed
done


# 5.4.3 generate the exon const bed file of each parent
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample=$(echo $line |cut -f1)
sbatch --export=sample="$sample" --output=log/ModifiedGenomeAnalysis_V2/ConstExon/"$sample".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/ConstExon/"$sample".sbatch.err /mnt/home/quanjian/quan/RNA/script/ConstExon_2.batch
done < RefName.txt

# 5.4.4 generate the genome coverge bed file of each parent.

for sample in `ls /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/`
do
cd /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/"$sample"/
cat <(awk '$6 == "+"' "$sample"_lift.const.exon.bed | sort -k1,1 -k2,2 | /mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools genomecov -bg -i stdin -g "$sample"_lift.fa.fai | awk '$4 > 1 {print $1"\t"$2"\t"$3"\t.\t.\t+"}') <(awk '$6 == "-"' "$sample"_lift.const.exon.bed | sort -k1,1 -k2,2 | /mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools genomecov -bg -i stdin -g "$sample"_lift.fa.fai | awk '$4 > 1 {print $1"\t"$2"\t"$3"\t.\t.\t-"}')  > "$sample"_lift.multiCov.bed
done

# 5.4.5 generate the 
for sample in `ls /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/`
do
cd /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/"$sample"/
/mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools subtract -a "$sample"_lift.const.exon.bed -b "$sample"_lift.multiCov.bed -s | sort -k1,1 -k2,2n > "$sample"_lift.const.exon.no.overlap.bed
done


# 5.5
# 5.5.1
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
bam_path=$(echo $line |cut -d " " -f 2)
bed_path=$(echo $line |cut -d " " -f 3)
read_path=$(echo $line |cut -d " " -f 4)
sbatch --export=sample_name="$sample_name",bam_path="$bam_path",bed_path="$bed_path",read_path="$read_path" --output=log/ModifiedGenomeAnalysis_V2/GetReadLength/"$sample_name".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/GetReadLength/"$sample_name".sbatch.err /mnt/home/quanjian/quan/RNA/script/GetReadLength.batch
done < GetReadLength.txt


# 5.5.2 get the gene assignment for the first read
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
bam_path=$(echo $line |cut -d " " -f 2)
bed_path=$(echo $line |cut -d " " -f 3)
overlap_bed_path=$(echo $line |cut -d " " -f 4)
read_len_path=$(echo $line |cut -d " " -f 5)

sbatch --export=sample_name="$sample_name",bam_path="$bam_path",bed_path="$bed_path",overlap_bed_path="$overlap_bed_path",read_len_path="$read_len_path" --output=log/ModifiedGenomeAnalysis_V2/GeneAssignment/FirstRead/"$sample_name".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/GeneAssignment/FirstRead/"$sample_name".sbatch.err /mnt/home/quanjian/quan/RNA/script/GeneAssignFirstRead.batch
done < GeneAssignment.txt


# 5.5.3 get the gene assignment for the second read
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
bam_path=$(echo $line |cut -d " " -f 2)
bed_path=$(echo $line |cut -d " " -f 3)
overlap_bed_path=$(echo $line |cut -d " " -f 4)
read_len_path=$(echo $line |cut -d " " -f 5)

sbatch --export=sample_name="$sample_name",bam_path="$bam_path",bed_path="$bed_path",overlap_bed_path="$overlap_bed_path",read_len_path="$read_len_path" --output=log/ModifiedGenomeAnalysis_V2/GeneAssignment/SecondRead/"$sample_name".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/GeneAssignment/SecondRead/"$sample_name".sbatch.err /mnt/home/quanjian/quan/RNA/script/GeneAssignSecondRead.batch
done < GeneAssignment.txt


# combine both reads, only consider one fragment if both reads map to the same gene
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample_name=$(echo $line |cut -d " " -f 1)
first_path=$(echo $line |cut -d " " -f 2)
second_path=$(echo $line |cut -d " " -f 3)

sbatch --export=sample_name="$sample_name",first_path="$first_path",second_path="$second_path" --output=log/ModifiedGenomeAnalysis_V2/GeneCount/"$sample_name".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/GeneCount/"$sample_name".sbatch.err /mnt/home/quanjian/quan/RNA/script/GeneCount.batch

done < GeneCount.txt
