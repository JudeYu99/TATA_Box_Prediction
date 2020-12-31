
################################### Perl 处理完 ######################################
library(pROC)

Dis50 <-  read.csv("Dis50.gff", sep = "\t", header = F)
temp <- subset(Dis50, select = c("V8", "V4"))
ROC_data <- data.frame(cbind(1:nrow(temp), rep("TATA"), temp))
colnames(ROC_data) <- c("ID", "Type",  "Attr", "ROX_Ct")

r1 <- roc(Attr~ROX_Ct, ROC_data)
r1
pdf("ROC_Dis50.pdf")
roc(ROC_data$Attr, ROC_data$ROX_Ct, plot = T, print.thres = T, print.auc = T)
dev.off()

################################################################################################################################################

Dis100 <-  read.csv("Dis100.gff", sep = "\t", header = F)
temp <- subset(Dis100, select = c("V8", "V4"))
ROC_data <- data.frame(cbind(1:nrow(temp), rep("TATA"), temp))
colnames(ROC_data) <- c("ID", "Type",  "Attr", "ROX_Ct")

r1 <- roc(Attr~ROX_Ct, ROC_data)
r1
pdf("ROC_Dis100.pdf")
roc(ROC_data$Attr, ROC_data$ROX_Ct, plot = T, print.thres = T, print.auc = T)
dev.off()

################################################################################################################################################

Dis150 <-  read.csv("Dis150.gff", sep = "\t", header = F)
temp <- subset(Dis150, select = c("V8", "V4"))
ROC_data <- data.frame(cbind(1:nrow(temp), rep("TATA"), temp))
colnames(ROC_data) <- c("ID", "Type",  "Attr", "ROX_Ct")

r1 <- roc(Attr~ROX_Ct, ROC_data)
r1
pdf("ROC_Dis150.pdf")
roc(ROC_data$Attr, ROC_data$ROX_Ct, plot = T, print.thres = T, print.auc = T)
dev.off()

################################################################################################################################################

Dis200 <-  read.csv("Dis200.gff", sep = "\t", header = F)
temp <- subset(Dis200, select = c("V8", "V4"))
ROC_data <- data.frame(cbind(1:nrow(temp), rep("TATA"), temp))
colnames(ROC_data) <- c("ID", "Type",  "Attr", "ROX_Ct")

r1 <- roc(Attr~ROX_Ct, ROC_data)
r1
pdf("ROC_Dis200.pdf")
roc(ROC_data$Attr, ROC_data$ROX_Ct, plot = T, print.thres = T, print.auc = T)
dev.off()

################################################################################################################################################

Dis250 <-  read.csv("Dis250.gff", sep = "\t", header = F)
temp <- subset(Dis250, select = c("V8", "V4"))
ROC_data <- data.frame(cbind(1:nrow(temp), rep("TATA"), temp))
colnames(ROC_data) <- c("ID", "Type",  "Attr", "ROX_Ct")

r1 <- roc(Attr~ROX_Ct, ROC_data)
r1
pdf("ROC_Dis250.pdf")
roc(ROC_data$Attr, ROC_data$ROX_Ct, plot = T, print.thres = T, print.auc = T)
dev.off()


