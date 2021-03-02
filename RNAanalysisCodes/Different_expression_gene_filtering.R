Differet_expression_gene_filtering <- function(input_data, groups){
  for(i in sort(groups)){
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
  
  
  
    #### 1）利用变化倍数进行初级过滤，logFC实际上为log2(FC),按照父母本比例大于6:4进行筛序
    filter1 <- as.data.frame(MA_PA_DEG) %>% filter(abs(logFC) > log(6/4,2))
  
    ### 根据两组合中的偏离程度进行筛选
    MA_PA_DEG_COUNT <- y$counts[rownames(filter1),]
    MA_PA_DEG_COUNT <- as.data.frame(MA_PA_DEG_COUNT)
  
    ## 提取组合1的父母本
    pair1_mat <- str_sub(colnames(MA_PA_DEG_COUNT)[1], -3,-1)
    pair1_pat <- str_sub(colnames(MA_PA_DEG_COUNT)[2], -3,-1)
    ## 提取组合2的父母本
    pair2_mat <-  str_sub(colnames(MA_PA_DEG_COUNT)[(a+b)*2], -3,-1)
    pair2_pat <-  str_sub(colnames(MA_PA_DEG_COUNT)[(a+b)*2-1], -3,-1)
  
  
  
    #### 2) 根据两组合表达量比值进行筛选,选择包含各父母本的样本的基因表达量数据
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
  
    ## 两个组合中至少在一个组合里要大于60/40，且另一个组合大于55/45
    tpm_dat <- ((rowSums(tpm2)/ rowSums(tpm1) > 55/45 & rowSums(tpm4)/ rowSums(tpm3) > 60/40)| rowSums(tpm2)/rowSums(tpm1) > 60/40 & rowSums(tpm4)/rowSums(tpm3) > 55/45) | 
    ((rowSums(tpm1)/rowSums(tpm2) > 55/45 & rowSums(tpm3)/rowSums(tpm4) > 60/40)|
       (rowSums(tpm1)/rowSums(tpm2) > 60/40 & rowSums(tpm3)/rowSums(tpm4)> 55/45)) 
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter2 <- tpm_dat %>% filter(ImbalancePer == TRUE)
  
  
  
    #### 3) 根据两组合中个体的父、母本表达量相同的次数进行筛选,去除父母本表达量相同过多的基因
    MA_PA_DEG_COUNT <- y$counts[rownames(filter2),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    ## 判断每个组合中每个基因的表达量是否是父母本相等
    tpm_dat <- merge(rowSums(tpm2 == tpm1), rowSums(tpm4 == tpm3), by= 0)
    ## 在两个组合中父母本基因表达量相同的个体数不得超过两个组合总个数的4分之1,或者在一个组合中不得存在父母本都相等。
    filter3 <- tpm_dat %>% filter(!((x == a | y == b) | (x+y) > 0.25*(a+b))) 
    
    
    
    #### 4) 根据两组合中个体的父母本优势表达方是否一致进行初步筛选
    MA_PA_DEG_COUNT <- y$counts[filter3$Row.names,]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    tpm_dat <- merge(rowSums(tpm2 > tpm1), rowSums(tpm4 > tpm3), by= 0)
    # 保留两组合中优势父母本不一致样本数不超过2的基因
    filter4 <- tpm_dat %>% filter((x+y<=2) | (x+y >= a+b-2)) 
  

    
    #### 5）根据父母本差异程度，进一步对两组合中个体的父母本非一致性基因进行筛选
    MA_PA_DEG_COUNT <- y$counts[filter4$Row.names,]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    ## 去除那些相反趋势，且在两个以及以上样本中存在相反且比值大于2的基因
    tpm_dat <- ((rowSums(tpm2)+rowSums(tpm4) > rowSums(tpm1) + rowSums(tpm3)) & (rowSums(tpm3/tpm4 > 2) + rowSums(tpm1/tpm2 > 2) <= 1) | ((rowSums(tpm1)+rowSums(tpm3) > rowSums(tpm2) + rowSums(tpm4)) & (rowSums(tpm4/tpm3 > 2) + rowSums(tpm2/tpm1 > 2) <= 1))) 
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter5 <- tpm_dat %>% filter(is.na(ImbalancePer) | ImbalancePer)
    
    
    #### 6）根据两组合中父母本一致性的样本数进行筛选（要求一致性的个体至少大于2个或者1个大于且1个相等）
    MA_PA_DEG_COUNT <- y$counts[rownames(filter5),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    tpm_dat <- ((rowSums(tpm2)+rowSums(tpm4) > rowSums(tpm1) + rowSums(tpm3)) & (rowSums(tpm2 >= tpm1) >=2 & rowSums(tpm4>=tpm3) >= 2))|((rowSums(tpm1)+rowSums(tpm3) > rowSums(tpm2) + rowSums(tpm4)) & (rowSums(tpm1 >= tpm2) >=2 & rowSums(tpm3>=tpm4) >= 2))
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter6 <- tpm_dat %>% filter(ImbalancePer)
    
    
    
    #### 7）根据两个组合中个体的父母本表达量平均值进行筛选
    MA_PA_DEG_COUNT <- y$counts[rownames(filter6),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    # 去除在任何组合中父母本表达量均过低的基因，以平均表达量为2作为条件阈值。
    tpm_dat <- (rowMeans(tpm1) > 2 | rowMeans(tpm2) > 2) & (rowMeans(tpm3) > 2 | rowMeans(tpm4) > 2) 
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter7 <- tpm_dat %>% filter(ImbalancePer)
    
    ### 导出
    MA_PA_DEG_COUNT <- y$counts[rownames(filter7),]
    write.csv(MA_PA_DEG_COUNT, file = paste("MA_PA_DEG_intersect_filtered",i,".csv", sep = ""), quote = F, row.names = T)
  }
}