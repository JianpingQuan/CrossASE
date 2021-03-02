
## 画条形堆叠图（fill）的代码
## 批量画条形堆叠图代码(position = "fill")
DEG_FillBarPlot <- function(DEG_dataframe, output_name){
  library(stringr) 
  library(ggplot2)
  DEG_t <- t(DEG_dataframe)
  pair1 <- sum(grepl(pattern = "^D", x = rownames(DEG_t))) / 2
  pair2 <- sum(grepl(pattern = "^L", x = rownames(DEG_t))) / 2
  Parent_of_Origin <- c(rep(c("Maternal","Paternal"), times = pair1),
                        rep(c("Paternal","Maternal"), times = pair2))  
  
  pdf(paste(output_name, "_fill.pdf", sep = ""))
  
  for(j in seq(1, ncol(DEG_t), by = 1)){
    Samples <- str_sub(rownames(DEG_t), 1, -8)
    Gene_name <- colnames(DEG_t)[j]
    Mapped_reads_count <- DEG_t[,j]
    data_draw <- data.frame(Samples, Mapped_reads_count, Parent_of_Origin)
    p <- ggplot(data = data_draw, mapping = aes(x = Samples, y = Mapped_reads_count, fill = Parent_of_Origin))+
      geom_bar(stat = "identity", position = "fill", width = 0.8, alpha = 0.8, color = "white") + 
      scale_fill_manual(values = c("#5CB85CFF","#9632B8FF"))+
      labs(title = Gene_name) +ylab("Mapped Reads Count") +theme_set(theme_bw())+ 
      theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = NULL, size = 9)) + 
      theme(panel.grid.major = element_line(linetype = "dotted"), panel.grid.minor = element_blank())+
      guides(fill=guide_legend(title="Parent of Origin"))
    print(p)
  }
  dev.off()
}
