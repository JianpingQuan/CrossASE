# 2、比对
# 2.1 数据无补测，只有一对paired reads
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample=$(echo $line | cut -d " " -f 1)
r1=$(echo $line | cut -d " " -f 2)
r2=$(echo $line | cut -d " " -f 3)
ref_path=$(echo $line | cut -d " " -f 4)
ref_name=$(echo $line | cut -d " " -f 5)

sbatch --export=sample="$sample",r1="$r1",r2="$r2",ref_path="$ref_path",ref_name="$ref_name" --output=log/ModifiedGenomeAnalysis_V2/Alignment/"$sample"_"$ref_name".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/Alignment/"$sample"_"$ref_name".sbatch.err /mnt/home/quanjian/quan/RNA/script/Hisat2_MG_Align_OnePair_V2.batch

done < MG_AlgSample_OnePair.txt


# 2.2 数据有补测，有两对paired reads
cd /mnt/home/quanjian/quan/RNA
while read line
do
sample=$(echo $line | cut -d " " -f 1)
r1_1=$(echo $line | cut -d " " -f 2)
r2_1=$(echo $line | cut -d " " -f 3)
r1_2=$(echo $line | cut -d " " -f 4)
r2_2=$(echo $line | cut -d " " -f 5)
ref_path=$(echo $line | cut -d " " -f 6)
ref_name=$(echo $line | cut -d " " -f 7)

sbatch --export=sample="$sample",ref_path="$ref_path",ref_name="$ref_name",r1_1="$r1_1",r2_1="$r2_1",r1_2="$r1_2",r2_2="$r2_2" --output=log/ModifiedGenomeAnalysis_V2/Alignment/"$sample"_"$ref_name".sbatch.log --error=log/ModifiedGenomeAnalysis_V2/Alignment/"$sample"_"$ref_name".sbatch.err /mnt/home/quanjian/quan/RNA/script/Hisat2_MG_Align_TwoPairs_V2.batch

done < MG_AlgSample_TwoPairs.txt