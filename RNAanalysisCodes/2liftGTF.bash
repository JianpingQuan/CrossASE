# 3 lift reference gtf and produce new reference genome/Sus_scrofa
# 3.1 copy the *prelift.bed file from /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2 into /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/"$sample"/

cd /mnt/scratch/quanjian/WGS/orginial_rawdata/Chr_allSample/sampleBed_V2/
for sample in `ls *_prelift.bed | cut -d "_" -f1`
do
cp "$sample"_prelift.bed /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/"$sample"/"$sample"_prelift.bed
done


# 3.2 run the liftGTF.pl to lift gtf and produce new reference based on above bed file.
cd /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/
for path in /mnt/scratch/quanjian/RNA_seq/F0_GenomeModified/
do
for sample in `realpath */*_prelift.bed | cut -d "/" -f10|cut -d "_" -f1`
do
perl /mnt/home/quanjian/quan/WGS/scripts/aseUtils/liftGTF.pl --gtf /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.98.gtf --bed "$path"/"$sample"/"$sample"_prelift.bed --fasta /mnt/research/qgg/resource/sscrofa11.1/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa --gtfout "$path"/"$sample"/"$sample"_lift.gtf --varout "$path"/"$sample"/"$sample"_lift.bed --fastaout "$path"/"$sample"/"$sample"_lift.fa 2> /mnt/home/quanjian/quan/RNA/log/liftGTF/"$sample"_lift.log &
done
done
