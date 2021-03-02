
#### 将每一个样本比对到两个亲本得到的 Informative reads count 展示成堆叠条形图。
## 数据处理部分
setwd("D:\\Documents\\geneCount\\20210223\\Intersect")
library(ggplot2)
library(ggrepel)

rawdata <- read.table("GeneCount_intersect210223_final.txt", header = T, row.names = 1) # 对于所有染色体

colnames(rawdata) <- gsub("_[BCEFGHIJLMNO]_", "_", colnames(rawdata)) # 去掉168天时期的样本序号
countSum <- colSums(rawdata) # 列求和
dat <- cbind(Sample = colnames(rawdata), SumCount = countSum)

dat <- as.data.frame(dat) # 将矩阵转化为数据框

## 根据sample中包含的Ref信息判断是母本（maternal）还是父本（paternal），然后将信息写入一个向量Ref_origin中
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

dat$Ref_Origin <- as.factor(Ref_Origin) # 将Ref_Origin转化为因子后添加进数据框
rownames(dat) <- NULL 
#dat$Sample <- sub('..................$',"", dat$Sample) # 去掉每一个样本名最后面的16个字符
dat$Sample <- sub('.......$',"", dat$Sample)
write.csv(dat, "AllSamplesIntersectGeneCount_plot.csv", quote = F, row.names = FALSE) 

#### 所有样本的所有常染色体一起
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
####对于所有常染色体

setwd("D:\\Documents\\geneCount\\20210223\\Intersect")
library(ggplot2)
library(ggrepel)

rawdata <- read.table("GeneCountIntersect210223_autosome.txt", header = T, row.names = 1) # 对于所有染色体

colnames(rawdata) <- gsub("_[BCEFGHIJLMNO]_", "_", colnames(rawdata)) # 去掉168天时期的样本序号
countSum <- colSums(rawdata) # 列求和
dat <- cbind(Sample = colnames(rawdata), SumCount = countSum)

dat <- as.data.frame(dat) # 将矩阵转化为数据框

## 根据sample中包含的Ref信息判断是母本（maternal）还是父本（paternal），然后将信息写入一个向量Ref_origin中
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

dat$Ref_Origin <- as.factor(Ref_Origin) # 将Ref_Origin转化为因子后添加进数据框
rownames(dat) <- NULL 
#dat$Sample <- sub('..................$',"", dat$Sample) # 去掉每一个样本名最后面的16个字符
dat$Sample <- sub('.......$',"", dat$Sample)
write.csv(dat, "AllSamplesIntersectGeneCount_Autosome_plot.csv", quote = F, row.names = FALSE) 

#### 所有样本的所有常染色体一起
dat <- read.csv("AllSamplesIntersectGeneCount_Autosome_plot.csv", header = T)

p <- ggplot(dat, aes(Sample, SumCount, fill= Ref_Origin)) 

p1 <- p + geom_bar(stat='identity',position='fill',width = 0.8) + scale_fill_brewer(palette = 'Accent') + labs(title = "The count of mapped reads on Autosome between two personalize genomes")

p2 <- p1 + theme_bw() + theme(plot.title = element_text(size = 15, colour = "red"), axis.text.x = element_text(angle = 90, vjust=0.3, size = 5), panel.grid.major = element_line(color = "red", size = 0.2))

p3 <-  p2 + scale_y_continuous(expand = c(0,0))

p4 <- p3 + geom_hline(aes(yintercept=0.55), color = "blue", alpha = 0.6) + geom_hline(aes(yintercept=0.45), color = "blue", alpha = 0.6)

pdf("AllSamplesIntersectGeneCount210223_autosome.pdf", width = 20, height = 5)
p4
dev.off()


