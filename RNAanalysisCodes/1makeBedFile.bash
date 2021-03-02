# 1. Generate the Bed file of Autosome and Allosome
## Autosome and X Chromosome

for vcffile in AllSamplesAuto AllSamplesChrX
do
perl /mnt/home/quanjian/quan/WGS/scripts/aseUtils/prepSeqBED.pl --vcf /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/"$vcffile".vcf --group1 DB1,DB2,DB3,DB4,DS1,DS2,DS3,DS4 --group2 LB1,LB2,LB3,LB4,LS1,LS2,LS3,LS4 --freq 0.99 --geno 0.874 2> /mnt/home/quanjian/quan/WGS/log/makeBedFile/"$vcffile".prepBED.err > /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2/"$vcffile".bed &
done

## Y Chromosome

perl /mnt/home/quanjian/quan/WGS/scripts/aseUtils/prepSeqBED.pl --vcf /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/AllSamplesChrY.vcf --group1 DB1,DB2,DB3,DB4 --group2 LB1,LB2,LB3,LB4 --freq 0.99 --geno 0.74 2> /mnt/home/quanjian/quan/WGS/log/makeBedFile/AllSamplesChrY.prepBED.err > /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2/AllSamplesChrY.bed &


# 2. Generate the Bed file of pair samples use for informative judgment(the Uppercase variants in Y chromosome are all informative. So the Bed file of Y chromosome can be added into after informative judgment)

## Merge the autosome bed and X chromosome bed
cat /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2/AllSamplesAuto.bed <(sed '1d' /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2/AllSamplesChrX.bed) > /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2/AllSamplesAutoX.bed

## Extract the bed of pair samples from Autosome and X chromosome. And it's safe to remove those where both alleles are the same as reference alleles
for path in /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2
do
cut -f 1-5,6,18 "$path"/AllSamplesAutoX.bed |awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/DB1_LS1_AutoX.bed
cut -f 1-5,7,19 "$path"/AllSamplesAutoX.bed |awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/DB2_LS2_AutoX.bed
cut -f 1-5,8,20 "$path"/AllSamplesAutoX.bed |awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/DB3_LS3_AutoX.bed
cut -f 1-5,9,21 "$path"/AllSamplesAutoX.bed |awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/DB4_LS4_AutoX.bed
cut -f 1-5,14,10 "$path"/AllSamplesAutoX.bed |awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/LB1_DS1_AutoX.bed
cut -f 1-5,15,11 "$path"/AllSamplesAutoX.bed |awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/LB2_DS2_AutoX.bed
cut -f 1-5,16,12 "$path"/AllSamplesAutoX.bed |awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/LB3_DS3_AutoX.bed
cut -f 1-5,17,13 "$path"/AllSamplesAutoX.bed |awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/LB4_DS4_AutoX.bed
done

## Extract the BED of chromosome Y in male animal.

#for path in /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2
#do
#for i in 6 7 8 9 10 11 12 13
#do
#cut -f 1-5,"$i" "$path"/AllSamplesChrY.bed | awk '!(toupper($6) == toupper($4))' > "$path"/"$i"_ChrY.bed
#done
#done
# Then change file name manually.

# 3. generate a gene list that were covered the informative variants. Only when both columns are uppercase, the variants are informative.
## For Autosome and X chromosome

for pair in DB1_LS1 DB2_LS2 DB3_LS3 DB4_LS4 LB1_DS1 LB2_DS2 LB3_DS3 LB4_DS4
do
awk '$6 ~ /[ATCG]/ && $7 ~ /[ATCG]/ && $6 != $7' /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2/"$pair"_AutoX.bed |/mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools intersect -a - -b <(awk '$3 == "exon"' /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.98.gtf) -wa -wb | cut -f 5,16 | awk -F "\"" '{print $1"\t"$2}' | awk '{print $3"\t"$1}' | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.27.1/bin/bedtools groupby -g 1 -c 2 -o count_distinct > /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2/"$pair"_AutoX_InfoCoverGene.txt
done


# 4 Generate the bed file of each sample for lift gtf. The bed files have no header and include all variants, not just informative ones. The fourth and fifth columns are the reference and alternative.
# For female

for path in /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2
do
for file in `echo DB*_AutoX.bed`
do
cut -f1-5,7 "$file" | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > "$path"/`echo "$file" | cut -d "_" -f 2`_prelift.bed
done
done


for path in /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2
do
for file in `echo LB*_AutoX.bed`
do
cut -f1-5,6 "$file" | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > "$path"/`echo "$file" | cut -d "_" -f 2`_prelift.bed
done
done



# For male(not consider the chromosome Y)

for path in /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2
do
for file in `echo DB*_AutoX.bed`
do
cut -f1-5,6 "$file" | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > "$path"/`echo "$file" | cut -d "_" -f1`_prelift.bed
done
done


for path in /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2
do
for file in `echo LB*_AutoX.bed`
do
cut -f1-5,7 "$file" | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > "$path"/`echo "$file" | cut -d "_" -f1`_prelift.bed
done
done


# For male(Consider the Chromosome Y)
#for path in /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2
#do
#for file in `echo *_AutoX.bed`
#do
#cat "$path"/`echo "$file" | cut -d "_" -f1`_ChrY.bed | tail -n+2 | cat <(cut -f1-6 "$path"/"$file") - | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > "$path"/`echo "$file" | cut -d "_" -f1`_prelift.bed
#done
#done


