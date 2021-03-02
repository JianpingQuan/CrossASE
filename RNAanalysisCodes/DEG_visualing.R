#################################
### �н�ɸѡ�õ��Ĳ�����������
#################################

##### ����͵��û�ͼ�������������ͼ
# ��������
path <- "D:\\Documents\\geneCount\\20210223\\Intersect"
fileNames <- dir(pattern = "MA_PA_DEG_intersect_filtered")
filePath <- sapply(fileNames, function(x){
  paste(path, x, sep = '/')}) ## ���ɶ�ȡ�ļ�·��
data <- lapply(filePath, function(x){
  read.csv(x, header = T, row.names = 1)})

# ��ʼ��ͼ
for(i in 1:length(data)){
  DEG <- data[[i]]
  File_name <- paste(str_sub(colnames(data[[i]])[1], -10,-9), str_sub(colnames(data[[i]])[1], 3,4), sep = "")
  DEG_FillBarPlot(DEG, File_name) ##����DEG_FillBarPlot()����
  DEG_BarPlot(DEG, File_name) ##����DEG_BarPlot()����
}