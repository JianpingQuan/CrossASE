##合并所有样本分布比对至父本和母本基因组上的gene count数据
module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604
files <- list.files()
listOfFiles <- lapply(files, function(x){ read.table(x, row.names=1)})
#The big reduction of all the files into a table
tbl <- Reduce(function(...) data.frame(merge(..., all = T, by = 0), row.names=1), listOfFiles)
tbl[is.na(tbl)] <- 0 #set all NA vals to 0
colnames(tbl) <- files #set the columns to the corresponding filenames (optional)
write.table(tbl,"GeneCount_intersect210223.txt", quote = F, row.names = T)