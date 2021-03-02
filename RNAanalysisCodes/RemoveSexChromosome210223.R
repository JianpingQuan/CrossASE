
###�Ƴ�ȫ����ȫ�����ļ���λ����Ⱦɫ��Ļ�������Ϊȫ������Ⱦɫ�����ļ���
setwd("D:\\Documents\\geneCount\\20210223\\Intersect")
library(dplyr)
library(stringr)

###�ҳ�ÿ������������Ⱦɫ����
gene_position <- read.csv("PigGenomeGenePosition.txt", header = F, sep = "\t")

colnames(gene_position) <- c("GeneName", "Chromosome", "Start", "End")

gene_position$Chromosome_noWhiteSpace <- as.numeric(gsub(" ","",gene_position$Chromosome)) ##��Ⱦɫ���ǿ�иĳ�����

gene_autosome <- gene_position %>% filter(Chromosome_noWhiteSpace >=1 & Chromosome_noWhiteSpace<=18)
gene_list <- gsub(" ","",gene_autosome$GeneName) ## ȥ���������ƺ�Ŀո�

###��ȡ��Ⱦɫ�ϵĻ���
gene_count <- read.csv("GeneCount_intersect210223_final.txt", header = T, sep = "\t")

col <- colnames(gene_count)
col[1] <- c("GeneName")
colnames(gene_count) <- col

gene_count_autosome <- filter(gene_count, GeneName %in% gene_list)

write.table(gene_count_autosome, "GeneCountIntersect210223_autosome.txt", quote = F, row.names = F)