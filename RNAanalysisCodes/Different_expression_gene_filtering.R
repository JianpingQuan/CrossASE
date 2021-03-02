Differet_expression_gene_filtering <- function(input_data, groups){
  for(i in sort(groups)){
    ## ���ݷ�������ݼ�����ȡ�Ӽ�
    rawDat <- input_data %>% dplyr::select(grep(i, names(input_data)))
    y = DGEList(rawDat)
  
    ## ���������У�D��ͷ����L��ͷ�������ĸ���
    dat_tpm <- y$samples
    dat_tpm$Rownames <- rownames(dat_tpm)
    a <- nrow(dat_tpm %>% dplyr::filter(grepl("^D", Rownames))) / 2 ##D��ͷ��������
    b <- nrow(dat_tpm %>% dplyr::filter(grepl("^L", Rownames))) / 2 ##L��ͷ��������
  
    # �����������������ø�ĸ���ظ���
    ref <- as.factor(c(rep(c("Maternal","Paternal"), times = a),rep(c("Paternal","Maternal"), times = b)))
  
    # ���ݷ�����������ȡ�����������ơ�
    if(grepl("168", i)){
      sample <- as.factor(substring(colnames(y), 1, 9))
    }else if(grepl("115", i)) {
      sample <- as.factor(substring(colnames(y), 1, 9))
    }else {
      sample <- as.factor(substring(colnames(y), 1, 8))
    }
  
    y$samples$ref <- ref
    y$samples$sample <- sample
  
    ### ���ȶ����������������У׼ (ӡ�ǻ������)
    design <- model.matrix(~sample+ref)
  
    ### ��������(ɸѡ ������3�������е�count per million (CPM) ����0�Ļ���)
    keep <- rowSums(cpm(y) > 1) >= 3
    y <- y[keep, ,keep.lib.sizes=FALSE]
  
    ### ���¼�����С
    y$samples$lib.size <- colSums(y$counts)
  
    ### ����TMM���������Ŀ��С�ı�׼��У��
    y <- calcNormFactors(y, method = "TMM")
    ## ������ɢ��
    y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
  
    ## fits the negative binomial GLM
    fit <- glmQLFit(y, design)
    ## carry out the QL F-test
    qlf <- glmQLFTest(fit, coef = "refPaternal")
  
    result <- summary(decideTests(qlf))
    ## ��ȡ�����������ϵ����µ��Ļ�������
    dif_num <- sum(result[1,1], result[3,1])
  
    ## ���ݲ������������ȡͳ�Ʊ��񣬲����ա�logFC������
    MA_PA_DEG <- topTags(qlf, n = dif_num, sort.by = "logFC", p.value = 0.05)
  
  
  
    #### 1�����ñ仯�������г������ˣ�logFCʵ����Ϊlog2(FC),���ո�ĸ����������6:4����ɸ��
    filter1 <- as.data.frame(MA_PA_DEG) %>% filter(abs(logFC) > log(6/4,2))
  
    ### ����������е�ƫ��̶Ƚ���ɸѡ
    MA_PA_DEG_COUNT <- y$counts[rownames(filter1),]
    MA_PA_DEG_COUNT <- as.data.frame(MA_PA_DEG_COUNT)
  
    ## ��ȡ���1�ĸ�ĸ��
    pair1_mat <- str_sub(colnames(MA_PA_DEG_COUNT)[1], -3,-1)
    pair1_pat <- str_sub(colnames(MA_PA_DEG_COUNT)[2], -3,-1)
    ## ��ȡ���2�ĸ�ĸ��
    pair2_mat <-  str_sub(colnames(MA_PA_DEG_COUNT)[(a+b)*2], -3,-1)
    pair2_pat <-  str_sub(colnames(MA_PA_DEG_COUNT)[(a+b)*2-1], -3,-1)
  
  
  
    #### 2) ��������ϱ�������ֵ����ɸѡ,ѡ���������ĸ���������Ļ������������
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
  
    ## ���������������һ�������Ҫ����60/40������һ����ϴ���55/45
    tpm_dat <- ((rowSums(tpm2)/ rowSums(tpm1) > 55/45 & rowSums(tpm4)/ rowSums(tpm3) > 60/40)| rowSums(tpm2)/rowSums(tpm1) > 60/40 & rowSums(tpm4)/rowSums(tpm3) > 55/45) | 
    ((rowSums(tpm1)/rowSums(tpm2) > 55/45 & rowSums(tpm3)/rowSums(tpm4) > 60/40)|
       (rowSums(tpm1)/rowSums(tpm2) > 60/40 & rowSums(tpm3)/rowSums(tpm4)> 55/45)) 
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter2 <- tpm_dat %>% filter(ImbalancePer == TRUE)
  
  
  
    #### 3) ����������и���ĸ���ĸ����������ͬ�Ĵ�������ɸѡ,ȥ����ĸ����������ͬ����Ļ���
    MA_PA_DEG_COUNT <- y$counts[rownames(filter2),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    ## �ж�ÿ�������ÿ������ı������Ƿ��Ǹ�ĸ�����
    tpm_dat <- merge(rowSums(tpm2 == tpm1), rowSums(tpm4 == tpm3), by= 0)
    ## ����������и�ĸ�������������ͬ�ĸ��������ó�����������ܸ�����4��֮1,������һ������в��ô��ڸ�ĸ������ȡ�
    filter3 <- tpm_dat %>% filter(!((x == a | y == b) | (x+y) > 0.25*(a+b))) 
    
    
    
    #### 4) ����������и���ĸ�ĸ�����Ʊ��﷽�Ƿ�һ�½��г���ɸѡ
    MA_PA_DEG_COUNT <- y$counts[filter3$Row.names,]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    tpm_dat <- merge(rowSums(tpm2 > tpm1), rowSums(tpm4 > tpm3), by= 0)
    # ��������������Ƹ�ĸ����һ��������������2�Ļ���
    filter4 <- tpm_dat %>% filter((x+y<=2) | (x+y >= a+b-2)) 
  

    
    #### 5�����ݸ�ĸ������̶ȣ���һ����������и���ĸ�ĸ����һ���Ի������ɸѡ
    MA_PA_DEG_COUNT <- y$counts[filter4$Row.names,]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    ## ȥ����Щ�෴���ƣ����������Լ����������д����෴�ұ�ֵ����2�Ļ���
    tpm_dat <- ((rowSums(tpm2)+rowSums(tpm4) > rowSums(tpm1) + rowSums(tpm3)) & (rowSums(tpm3/tpm4 > 2) + rowSums(tpm1/tpm2 > 2) <= 1) | ((rowSums(tpm1)+rowSums(tpm3) > rowSums(tpm2) + rowSums(tpm4)) & (rowSums(tpm4/tpm3 > 2) + rowSums(tpm2/tpm1 > 2) <= 1))) 
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter5 <- tpm_dat %>% filter(is.na(ImbalancePer) | ImbalancePer)
    
    
    #### 6������������и�ĸ��һ���Ե�����������ɸѡ��Ҫ��һ���Եĸ������ٴ���2������1��������1����ȣ�
    MA_PA_DEG_COUNT <- y$counts[rownames(filter5),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    tpm_dat <- ((rowSums(tpm2)+rowSums(tpm4) > rowSums(tpm1) + rowSums(tpm3)) & (rowSums(tpm2 >= tpm1) >=2 & rowSums(tpm4>=tpm3) >= 2))|((rowSums(tpm1)+rowSums(tpm3) > rowSums(tpm2) + rowSums(tpm4)) & (rowSums(tpm1 >= tpm2) >=2 & rowSums(tpm3>=tpm4) >= 2))
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter6 <- tpm_dat %>% filter(ImbalancePer)
    
    
    
    #### 7��������������и���ĸ�ĸ��������ƽ��ֵ����ɸѡ
    MA_PA_DEG_COUNT <- y$counts[rownames(filter6),]
    tpm1 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_mat))
    tpm2 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair1_pat))
    tpm3 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_mat))
    tpm4 <- as.data.frame(MA_PA_DEG_COUNT) %>% select(contains(pair2_pat))
    # ȥ�����κ�����и�ĸ�������������͵Ļ�����ƽ��������Ϊ2��Ϊ������ֵ��
    tpm_dat <- (rowMeans(tpm1) > 2 | rowMeans(tpm2) > 2) & (rowMeans(tpm3) > 2 | rowMeans(tpm4) > 2) 
    tpm_dat <- as.data.frame(tpm_dat)
    colnames(tpm_dat) <- c("ImbalancePer")
    filter7 <- tpm_dat %>% filter(ImbalancePer)
    
    ### ����
    MA_PA_DEG_COUNT <- y$counts[rownames(filter7),]
    write.csv(MA_PA_DEG_COUNT, file = paste("MA_PA_DEG_intersect_filtered",i,".csv", sep = ""), quote = F, row.names = T)
  }
}