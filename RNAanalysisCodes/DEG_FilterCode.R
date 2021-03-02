
### 筛选流程
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


### 数据清洗过程
imbalance_gene_filtering <- function(input_data, groups){
  for(i in groups){
    ## 根据分组从数据集中提取子集
    rawDat <- input_data %>% dplyr::select(grep(i, names(input_data)))
    y = DGEList(rawDat)
    
    ## 计算各组别中，D开头的与L开头的样本的个数
    dat_tpm <- y$samples
    dat_tpm$Rownames <- rownames(dat_tpm)
    a <- nrow(dat_tpm %>% dplyr::filter(grepl("^D", Rownames))) / 2 ##D开头的样本数
    b <- nrow(dat_tpm %>% dplyr::filter(grepl("^L", Rownames))) / 2 ##L开头的样本数
  
    # 根据样本个数，设置父母本重复数
    ref <- as.factor(c(rep(c("Maternal","Paternal"), times = a),rep(c("Paternal","Maternal"), times = b)))

    # 根据分组的情况，提取出样本的名称。
    if(grepl("168", i)){
        sample <- as.factor(substring(colnames(y), 1, 9))
    }else if(grepl("115", i)) {
        sample <- as.factor(substring(colnames(y), 1, 9))
    }else {
        sample <- as.factor(substring(colnames(y), 1, 8))
    }

    y$samples$ref <- ref
    y$samples$sample <- sample

    ### 首先对样本本身差异进行校准 (印记基因查找)
    design <- model.matrix(~sample+ref)
  
    ### 初步过滤(筛选 至少在3个样本中的count per million (CPM) 大于0的基因)
    keep <- rowSums(cpm(y) > 1) >= 3
    y <- y[keep, ,keep.lib.sizes=FALSE]

    ### 重新计算库大小
    y$samples$lib.size <- colSums(y$counts)

    ### 利用TMM方法进行文库大小的标准化校正
    y <- calcNormFactors(y, method = "TMM")
    ## 评估离散度
    y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    
    ## fits the negative binomial GLM
    fit <- glmQLFit(y, design)
    ## carry out the QL F-test
    qlf <- glmQLFTest(fit, coef = "refPaternal")
    
    result <- summary(decideTests(qlf))
    ## 提取表达量显著上调或下调的基因数量
    dif_num <- sum(result[1,1], result[3,1])
    
    ## 根据差异表达数量提取统计表格，并按照“logFC”排序
    MA_PA_DEG <- topTags(qlf, n = dif_num, sort.by = "logFC", p.value = 0.05)
    
    # 利用变化倍数进行初级过滤，logFC实际上为log2(FC),按照父母本比例大于6:4进行筛序
    filter1 <- as.data.frame(MA_PA_DEG) %>% filter(abs(logFC) > 6/4)
    

    #### 根据两组合中的偏离程度进行筛选
    MA_PA_DEG_COUNT <- y$counts[rownames(filter1),]
    MA_PA_DEG_COUNT <- as.data.frame(MA_PA_DEG_COUNT)
    
    ## 提取
    pair1_mat <- str_sub(colnames(MA_PA_DEG_COUNT)[1], -3,-1)
    pair1_pat <- str_sub(colnames(MA_PA_DEG_COUNT)[2], -3,-1)
    
    pair2_mat <-  str_sub(colnames(MA_PA_DEG_COUNT)[(a+b)*2], -3,-1)
    pair2_pat <-  str_sub(colnames(MA_PA_DEG_COUNT)[(a+b)*2-1], -3,-1)


    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    
    # 至少在一个组合里要大于60/40，且另一个组合大于55/45
    tpm_dat <- ((rowSums(tpm2)/ rowSums(tpm1) > 55/45 & rowSums(tpm4)/ rowSums(tpm3) > 60/40)| rowSums(tpm2)/rowSums(tpm1) > 60/40 & rowSums(tpm4)/rowSums(tpm3) > 55/45) | 
      ((rowSums(tpm1)/rowSums(tpm2) > 55/45 & rowSums(tpm3)/rowSums(tpm4) > 60/40)|
      (rowSums(tpm1)/rowSums(tpm2) > 60/40 & rowSums(tpm3)/rowSums(tpm4)> 55/45)) 
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter2 <- tpm_dat %>% filter(ImbalancePer == TRUE)


    #### 根据两组合中个体的父、母本是否相同进行筛选
    MA_PA_DEG_COUNT <- y$counts[rownames(filter2),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    tpm_dat <- merge(rowSums(tpm2 == tpm1), rowSums(tpm4 == tpm3), by= 0)
    filter3 <- tpm_dat %>% filter(!((x == a | y == b) | (x+y) > 0.25*(a+b))) 
    # 去除父母本相同过多的基因(5/25)
    
    
    #### 根据两组合中个体的父母本是否具有良好一致性进行初步筛选
    MA_PA_DEG_COUNT <- y$counts[filter3$Row.names,]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    tpm_dat <- merge(rowSums(tpm2 > tpm1), rowSums(tpm4 > tpm3), by= 0)
    filter4 <- tpm_dat %>% filter((x+y<=2) | (x+y >= a+b-2)) 
    # 去除两组合中优势父母本不一致的基因(3/20)
    
    
    #### 根据两组合中个体的父母本非一致性样本中，父母本差异程度进行筛选
    MA_PA_DEG_COUNT <- y$counts[filter4$Row.names,]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    ## 去除那些相反趋势，且在两个以及以上样本中存在相反且比值大于2的基因
    tpm_dat <- ((rowSums(tpm2)+rowSums(tpm4) > rowSums(tpm1) + rowSums(tpm3)) & (rowSums(tpm3/tpm4 > 2) + rowSums(tpm1/tpm2 > 2) <= 1) | 
    ((rowSums(tpm1)+rowSums(tpm3) > rowSums(tpm2) + rowSums(tpm4)) & (rowSums(tpm4/tpm3 > 2) + rowSums(tpm2/tpm1 > 2) <= 1))) 
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter5 <- tpm_dat %>% filter(is.na(ImbalancePer) | ImbalancePer) # 保留有NA的行
    
    
    #### 根据两组合中父母本一致性的样本数量进行筛选（要求一致性的个体至少大于2个或者1个大于且1个相等）
    MA_PA_DEG_COUNT <- y$counts[rownames(filter5),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    tpm_dat <- ((rowSums(tpm2)+rowSums(tpm4) > rowSums(tpm1) + rowSums(tpm3)) & (rowSums(tpm2 >= tpm1) >=2 & rowSums(tpm4>=tpm3) >= 2))|
      ((rowSums(tpm1)+rowSums(tpm3) > rowSums(tpm2) + rowSums(tpm4)) & (rowSums(tpm1 >= tpm2) >=2 & rowSums(tpm3>=tpm4) >= 2))
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter6 <- tpm_dat %>% filter(ImbalancePer)
    
    
    #### 根据两个组合中个体的父母本丰度平均值进行筛选
    MA_PA_DEG_COUNT <- y$counts[rownames(filter6),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    tpm_dat <- (rowMeans(tpm1) > 2 | rowMeans(tpm2) > 2) & (rowMeans(tpm3) > 2 | rowMeans(tpm4) > 2) # 去除在一个任何一个组合中丰度过低的基因
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter7 <- tpm_dat %>% filter(ImbalancePer)
    
    ### 导出
    MA_PA_DEG_COUNT <- y$counts[rownames(filter7),]
    write.csv(MA_PA_DEG_COUNT, file = paste("MA_PA_DEG_intersect_filtered",i,".csv", sep = ""), quote = F, row.names = T)
    }
}

###使用方式
# imbalance_gene_filtering(rawdata, stage_tissue)
# imbalance_gene_filtering(rawdata, stage_tissue_del)




### 文件的批量读取
path <- "D:\\Documents\\geneCount\\20210128intersect"

fileNames <- dir(pattern = "MA_PA_DEG_intersect_filtered")

filePath <- sapply(fileNames, function(x){
  paste(path, x, sep = '/')}) ## 生成读取文件路径

data <- lapply(filePath, function(x){
  read.csv(x, header = T)})  ## 读取数据，结果为列表


##将列表中每一个数据框的第一列输出并存在一个向量中。
DEG <- c()
for(i in 1:length(data)){
  DEG <- c(DEG,data[[i]]$X)
}
seed_genes <- as.data.frame(table(sort(DEG))) ## 统计每一个基因出现了多少次




##### 读入和调用画图代码进行批量画图
# 批量导入
path <- "D:\\Documents\\geneCount\\20210128intersect"
fileNames <- dir(pattern = "MA_PA_DEG_intersect_filtered")
filePath <- sapply(fileNames, function(x){
  paste(path, x, sep = '/')}) ## 生成读取文件路径
data <- lapply(filePath, function(x){
  read.csv(x, header = T, row.names = 1)})

# 开始画图
for(i in 1:length(data)){
  DEG <- data[[i]]
  File_name <- paste(str_sub(colnames(data[[i]])[1], -10,-9), str_sub(colnames(data[[i]])[1], 3,4), sep = "")
  DEG_FillBarPlot(DEG, File_name)
  DEG_BarPlot(DEG, File_name)
}
