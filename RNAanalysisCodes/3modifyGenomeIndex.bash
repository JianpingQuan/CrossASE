# 1.1 Use the hisat2_extract_splice_sites.py to extract the splice sites from modified genome gtf file.
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample=$(echo $line | cut -d "/" -f 10|cut -d "_" -f1)
path=$(echo $line)
sbatch --export=sample="$sample",gtf_path="$path" --output=log/ModifiedGenomeAnalysis_V2/ExtractSpliceSite/"$sample".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/ExtractSpliceSite/"$sample".sbatch.err /mnt/home/quanjian/quan/RNA/script/Hisat2_MG_ExtractSpliceSite.batch
done < MG_sample_gtf_V2.txt

# 1.2 Use the hisat2_extract_splice_sites.py to extract the exons from modified genome gtf file.
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample=$(echo $line | cut -d "/" -f 10|cut -d "_" -f1)
path=$(echo $line)
sbatch --export=sample="$sample",gtf_path="$path" --output=log/ModifiedGenomeAnalysis_V2/ExtractExon/"$sample".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/ExtractExon/"$sample".sbatch.err /mnt/home/quanjian/quan/RNA/script/Hisat2_MG_ExtractExon.batch
done < MG_sample_gtf_V2.txt


# 1.3 Generate the index file combine the splice and exons files.
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample=$(echo $line |cut -f1|cut -d "/" -f10|cut -d "_" -f1)
fasta=$(echo $line|cut -d " " -f1)
ss=$(echo $line|cut -d " " -f2)
exon=$(echo $line|cut -d " " -f3)
sbatch --export=sample="$sample",fasta_path="$fasta",ss_path="$ss",exon_path="$exon" --output=log/ModifiedGenomeAnalysis_V2/Hisat2BuildIndex/"$sample".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/Hisat2BuildIndex/"$sample".sbatch.err /mnt/home/quanjian/quan/RNA/script/Hisat2_MG_indexBuild_V2.batch
done < MG_sample_fasta_ss_exon_V2.txt
