
###移除全样本全基因文件中位于性染色大的基因，另存为全样本常染色基因文件。
setwd("D:\\Documents\\geneCount\\20210223\\Intersect")
library(dplyr)
library(stringr)

###找出每个基因所属的染色体编号
gene_position <- read.csv("PigGenomeGenePosition.txt", header = F, sep = "\t")

colnames(gene_position) <- c("GeneName", "Chromosome", "Start", "End")

gene_position$Chromosome_noWhiteSpace <- as.numeric(gsub(" ","",gene_position$Chromosome)) ##将染色体号强行改成数字

gene_autosome <- gene_position %>% filter(Chromosome_noWhiteSpace >=1 & Chromosome_noWhiteSpace<=18)
gene_list <- gsub(" ","",gene_autosome$GeneName) ## 去掉基因名称后的空格

###提取常染色上的基因
gene_count <- read.csv("GeneCount_intersect210223_final.txt", header = T, sep = "\t")

col <- colnames(gene_count)
col[1] <- c("GeneName")
colnames(gene_count) <- col

gene_count_autosome <- filter(gene_count, GeneName %in% gene_list)

write.table(gene_count_autosome, "GeneCountIntersect210223_autosome.txt", quote = F, row.names = F)
