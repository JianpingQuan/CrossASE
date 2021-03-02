#################################
### 承接筛选得到的差异表达基因结果
#################################

##### 读入和调用画图代码进行批量画图
# 批量导入
path <- "D:\\Documents\\geneCount\\20210223\\Intersect"
fileNames <- dir(pattern = "MA_PA_DEG_intersect_filtered")
filePath <- sapply(fileNames, function(x){
  paste(path, x, sep = '/')}) ## 生成读取文件路径
data <- lapply(filePath, function(x){
  read.csv(x, header = T, row.names = 1)})

# 开始画图
for(i in 1:length(data)){
  DEG <- data[[i]]
  File_name <- paste(str_sub(colnames(data[[i]])[1], -10,-9), str_sub(colnames(data[[i]])[1], 3,4), sep = "")
  DEG_FillBarPlot(DEG, File_name) ##调用DEG_FillBarPlot()函数
  DEG_BarPlot(DEG, File_name) ##调用DEG_BarPlot()函数
}