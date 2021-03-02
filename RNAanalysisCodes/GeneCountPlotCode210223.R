
#### ��ÿһ�������ȶԵ������ױ��õ��� Informative reads count չʾ�ɶѵ�����ͼ��
## ���ݴ�������
setwd("D:\\Documents\\geneCount\\20210223\\Intersect")
library(ggplot2)
library(ggrepel)

rawdata <- read.table("GeneCount_intersect210223_final.txt", header = T, row.names = 1) # ��������Ⱦɫ��

colnames(rawdata) <- gsub("_[BCEFGHIJLMNO]_", "_", colnames(rawdata)) # ȥ��168��ʱ�ڵ��������
countSum <- colSums(rawdata) # �����
dat <- cbind(Sample = colnames(rawdata), SumCount = countSum)

dat <- as.data.frame(dat) # ������ת��Ϊ���ݿ�

## ����sample�а�����Ref��Ϣ�ж���ĸ����maternal�����Ǹ�����paternal����Ȼ����Ϣд��һ������Ref_origin��
Ref_Origin <- c() 
for(RN in 1:nrow(dat)){

	if(grepl(pattern = "RefDS", x = as.character(dat$Sample)[RN])) {
		Ref_Origin[RN] = "Maternal"
	}else if(grepl(pattern = "RefLS", x = as.character(dat$Sample)[RN])){
		Ref_Origin[RN] = "Maternal"
	} else{
		Ref_Origin[RN] = "Paternal"
	}
}

dat$Ref_Origin <- as.factor(Ref_Origin) # ��Ref_Originת��Ϊ���Ӻ����ӽ����ݿ�
rownames(dat) <- NULL 
#dat$Sample <- sub('..................$',"", dat$Sample) # ȥ��ÿһ��������������16���ַ�
dat$Sample <- sub('.......$',"", dat$Sample)
write.csv(dat, "AllSamplesIntersectGeneCount_plot.csv", quote = F, row.names = FALSE) 

#### �������������г�Ⱦɫ��һ��
dat <- read.csv("AllSamplesIntersectGeneCount_plot.csv", header = T)

p <- ggplot(dat, aes(Sample, SumCount, fill= Ref_Origin)) 

p1 <- p + geom_bar(stat='identity',position='fill',width = 0.8) + scale_fill_brewer(palette = 'Accent') + labs(title = "The count of mapped reads on All Chromosome between two personalize genomes")

p2 <- p1 + theme_bw() + theme(plot.title = element_text(size = 15, colour = "red"), axis.text.x = element_text(angle = 90, vjust=0.3, size = 5), panel.grid.major = element_line(color = "red", size = 0.2))

p3 <-  p2 + scale_y_continuous(expand = c(0,0))

p4 <- p3 + geom_hline(aes(yintercept=0.55), color = "blue", alpha = 0.6) + geom_hline(aes(yintercept=0.45), color = "blue", alpha = 0.6)

pdf("AllSamplesIntersectGeneCount210223.pdf", width = 20, height = 5)
p4
dev.off()



##################################################
####�������г�Ⱦɫ��

setwd("D:\\Documents\\geneCount\\20210223\\Intersect")
library(ggplot2)
library(ggrepel)

rawdata <- read.table("GeneCountIntersect210223_autosome.txt", header = T, row.names = 1) # ��������Ⱦɫ��

colnames(rawdata) <- gsub("_[BCEFGHIJLMNO]_", "_", colnames(rawdata)) # ȥ��168��ʱ�ڵ��������
countSum <- colSums(rawdata) # �����
dat <- cbind(Sample = colnames(rawdata), SumCount = countSum)

dat <- as.data.frame(dat) # ������ת��Ϊ���ݿ�

## ����sample�а�����Ref��Ϣ�ж���ĸ����maternal�����Ǹ�����paternal����Ȼ����Ϣд��һ������Ref_origin��
Ref_Origin <- c() 
for(RN in 1:nrow(dat)){
  
  if(grepl(pattern = "RefDS", x = as.character(dat$Sample)[RN])) {
    Ref_Origin[RN] = "Maternal"
  }else if(grepl(pattern = "RefLS", x = as.character(dat$Sample)[RN])){
    Ref_Origin[RN] = "Maternal"
  } else{
    Ref_Origin[RN] = "Paternal"
  }
}

dat$Ref_Origin <- as.factor(Ref_Origin) # ��Ref_Originת��Ϊ���Ӻ����ӽ����ݿ�
rownames(dat) <- NULL 
#dat$Sample <- sub('..................$',"", dat$Sample) # ȥ��ÿһ��������������16���ַ�
dat$Sample <- sub('.......$',"", dat$Sample)
write.csv(dat, "AllSamplesIntersectGeneCount_Autosome_plot.csv", quote = F, row.names = FALSE) 

#### �������������г�Ⱦɫ��һ��
dat <- read.csv("AllSamplesIntersectGeneCount_Autosome_plot.csv", header = T)

p <- ggplot(dat, aes(Sample, SumCount, fill= Ref_Origin)) 

p1 <- p + geom_bar(stat='identity',position='fill',width = 0.8) + scale_fill_brewer(palette = 'Accent') + labs(title = "The count of mapped reads on Autosome between two personalize genomes")

p2 <- p1 + theme_bw() + theme(plot.title = element_text(size = 15, colour = "red"), axis.text.x = element_text(angle = 90, vjust=0.3, size = 5), panel.grid.major = element_line(color = "red", size = 0.2))

p3 <-  p2 + scale_y_continuous(expand = c(0,0))

p4 <- p3 + geom_hline(aes(yintercept=0.55), color = "blue", alpha = 0.6) + geom_hline(aes(yintercept=0.45), color = "blue", alpha = 0.6)

pdf("AllSamplesIntersectGeneCount210223_autosome.pdf", width = 20, height = 5)
p4
dev.off()

