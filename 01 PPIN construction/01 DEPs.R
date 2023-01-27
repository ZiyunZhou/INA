## The differential expression analysis of primary and recurrent proteome.

library(pheatmap)
library(readxl)
library(ggplot2)

data1_heat <- read_excel("./proteome.xlsx",col_names = TRUE)

m <- nrow(data1_heat)
data2_heat <- data.frame()
for(i in 1:m) {
  if(mean(data1_heat[i,3:32] != 0)) {
    data2_heat <- rbind(data2_heat,data1_heat[i,])
  }
}
colnames(data2_heat) <- c("UniprotID", "GeneID", paste("TET_C", 1:15, sep = ""),paste("TET_P", 1:15, sep = ""))


p_value <- c()
logFC <- c()
absFC <- c()
fc <- c()


for(i in 1:nrow(data2_heat)) {
  # T-test
  test <- t.test(data2_heat[i,3:17], data2_heat[i,18:32])
  p_value[i] <- test$p.value
  # judge the P-value
  FC <- mean(as.matrix(data2_heat[i,3:17])) / mean(as.matrix(data2_heat[i,18:32]))
  absFC[i] <- abs(mean(as.matrix(data2_heat[i,3:17])) / mean(as.matrix(data2_heat[i,18:32])))
  logFC[i] <- log2(FC)
}
p_adjust <- p.adjust(p_value, method = "BH", n = length(p_value))
data2_heat$p_value <- p_value
data2_heat$p_adjust <- p_adjust
data2_heat$logp <- -log10(data2_heat$p_value)
data2_heat$logFC <- logFC
data2_heat$absFC <- absFC

data_cancer_noncancer_menu2_up <- data2_heat[p_value < 0.05 & absFC > 1.5,]  #1222
data_cancer_noncancer_menu2_down <- data2_heat[p_value < 0.05 & absFC < 2/3,]  #122

diff <- rbind(data_cancer_noncancer_menu2_up,data_cancer_noncancer_menu2_down)

write.table(diff[,c(2,33:37)],"./primary_diff_protein.txt", quote = F, col.names = T, row.names = F, sep = "\t")


# heatmap
pheatmap(as.matrix(rbind(data_cancer_noncancer_menu2_up[,3:32],data_cancer_noncancer_menu2_down[,3:32])),
         scale="row",
         treeheight_row=100,
         cluster_cols=FALSE,
         color=colorRampPalette(c("green","black","red"))(1000),
         border_color=NA,
         fontsize_row=2,
         fontsize_col=8,
         file = './primary_diff_heatmap.pdf')

#volvano-plot
# ddfine a value to judge the colour
data2 = data2_heat[,c(35,36)]
data2$threshold[data2[,2] > log2(1.5) & data2[,1] > -log10(0.05) ] = "up"
data2$threshold[data2[,2] < log2(2/3) & data2[,1] > -log10(0.05) ] = "down"
data2$threshold[is.na(data2$threshold)] <- "non"
data2 <- as.data.frame(data2)

pdf("primary_diff_volvano.pdf")
ggplot(data = data2, aes(x = data2[,2], y = data2[,1], colour = threshold))+
  xlab("log2(Fold Change)") +
  ylab("-log10(P_Value)") +
  theme_bw()+
  geom_point(alpha = 2, size = 1.6) +
  labs(title = "Cancer-Noncancer") +
  scale_color_manual(values = c("non" = "grey", "up" = "red", "down" = "blue"))
dev.off()

m <- nrow(data1_heat)
data2_heat <- data.frame()
for(i in 1:m) {
  if(sd(data1_heat[i,3:7]) != 0 || sd(data1_heat[i,8:17]) != 0) {
    data2_heat <- rbind(data2_heat,data1_heat[i,])
  }
}
p_value <- c()
logFC <- c()
absFC <- c()

for(i in 1:nrow(data2_heat)) {
  # T-test
  test <- t.test(data2_heat[i,3:7], data2_heat[i,8:17])
  p_value[i] <- test$p.value
  # judge the P-value
  FC <- mean(as.matrix(data2_heat[i,3:7])) / mean(as.matrix(data2_heat[i,8:17]))
  absFC[i] <- abs(mean(as.matrix(data2_heat[i,3:7])) / mean(as.matrix(data2_heat[i,8:17])))
  logFC[i] <- log2(FC)
}
p_adjust <- p.adjust(p_value, method = "BH", n = length(p_value))
data2_heat$p_value <- p_value
data2_heat$p_adjust <- p_adjust
data2_heat$logp <- -log10(data2_heat$p_value)
data2_heat$logFC <- logFC
data2_heat$absFC <- absFC

data_cancer_recurrence_menu2_up <- data2_heat[p_value < 0.05 & absFC > 1.5,]  #39
data_cancer_recurrence_menu2_down <- data2_heat[p_value < 0.05 & absFC < 2/3,]   #32

diff <- rbind(data_cancer_recurrence_menu2_up,data_cancer_recurrence_menu2_down)

write.table(diff[,c(2,33:37)],"./recurrent_diff_protein.txt", quote = F, col.names = T, row.names = F, sep = "\t")

pheatmap(as.matrix(rbind(data_cancer_recurrence_menu2_up[,3:17],data_cancer_recurrence_menu2_down[,3:17])),
         scale="row",
         treeheight_row=100,
         cluster_cols=FALSE,
         color=colorRampPalette(c("green","black","red"))(1000),
         border_color=NA,
         fontsize_row=2,
         fontsize_col=8,
         file = './recurrent_diff_heatmap.pdf')

#volvano-plot
# ddfine a value to judge the colour
data2 = data2_heat[,c(35,36)]
data2$threshold[data2[,2] > log2(1.5) & data2[,1] > -log10(0.05) ] = "up"
data2$threshold[data2[,2] < log2(2/3) & data2[,1] > -log10(0.05) ] = "down"
data2$threshold[is.na(data2$threshold)] <- "non"
data2 <- as.data.frame(data2)
#data2$threshold[data2[,1] < 0.05 & (data2[,2] >= log2(1.5) | data2[,2] <= log2(2/3)) ] = "non"

pdf("recurrent_diff_volvano.pdf")
ggplot(data = data2, aes(x = data2[,2], y = data2[,1], colour = threshold))+
  xlab("log2(Fold Change)") +
  ylab("-log10(P_Value)") +
  theme_bw()+
  geom_point(alpha = 2, size = 1.6) +
  labs(title = "cancer_recurrence") +
  scale_color_manual(values = c("non" = "grey", "up" = "red", "down" = "blue"))
dev.off()



