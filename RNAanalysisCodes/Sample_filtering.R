### 样本筛选
setwd("D:\\Documents\\geneCount\\20210223\\Intersect")
library(dplyr)
library(ggplot2)
library(edgeR)
library(limma)
library(stringr)
rawdata <- read.table("GeneCountIntersect210223_autosome.txt",header = T, row.names = 1)
colnames(rawdata) <- gsub("_[BCEFGHIJLMNO]_", "_", colnames(rawdata))

# 需要删除父母本存在严重偏离的样本、
## 去除部分样本后，需要重新分析的组别如下：
sample_del <- c("DB40_Pl1", "DB40_Pl2", "DB40_Pl3","DB70_Pl1", "DS168_Br1","DS40_Pl1", "LB70_Pl1","LS168_Br1", "LS168_Br2", "LS168_Mu1", "LS168_Mu2", "LS168_Li1", "LS168_Li2")

stage_tissue_del <- c("168_Br", "168_Li", "168_Mu", "40_Pl", "70_Pl")

stage_tissue <- c("40_Br","40_Mu","40_Li","70_Br","70_Mu","70_Li","115_Br","115_Mu","115_Li")

## 删除存在明显父母本偏离的样本
for(del_sample in sample_del){
  rawdata <- rawdata %>% dplyr::select(-contains(del_sample))
}

##调用函数Different_expression_gene_filtering()

###使用方式
# Different_expression_gene_filtering(rawdata, stage_tissue)
# Different_expression_gene_filtering(rawdata, stage_tissue_del)