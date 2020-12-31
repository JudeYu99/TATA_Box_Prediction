# 读取预测结果文件
data <- read.table("TATA_distance.gff")
colnames(data) <- c("chr_ID", "Start", "End", "Score", "p_value", "Chain", "Distance")

write.table(data, "Plot_data.csv", quote = F, row.names = F, sep = "\t")

# 直方图绘制
pdf("Distance_Frequency.pdf") 
hist(data$Distance, xlim = c(0, 20000), freq = T, main = "Distance Frequency Distribution", xlab = "Distance") 
dev.off()

# 散点图绘制
library(ggplot2)

title <- "Saccharomyces cerevisiae YJM993 all chromosomes"
p <- ggplot(data, aes(x = Distance, y = Score, color = Chain)) + 
            geom_point(alpha = 1/3) + 
            labs(title = title) +
            theme(plot.title = element_text(hjust = 0.5)) +
            xlab("Distance") +
            ylab("Score")
ggsave(p, filename = "YJM993_Distance_Score.pdf", width = 7, height = 6)
