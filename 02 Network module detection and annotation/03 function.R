setwd("G:/TET/网络重建")


####### GO富集

library(org.Hs.eg.db) #人类注释数据库
library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)
library(enrichplot)

f <- read.table("./5 模块划分/primary_node_Module.txt", header=T, sep='\t')
m <- f[which(f[,2]=="16"),1]

go <- enrichGO(m, OrgDb = org.Hs.eg.db, ont='ALL', 
               pAdjustMethod = 'BH', pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2, keyType = 'SYMBOL')
write.csv(go, "./6 富集分析/primary-go16.csv")

pdf("./6 富集分析/primary-go16.pdf", width = 10, height = 10)
dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~.,scale="free")
dev.off()



####### KEGG富集
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(stringr)
# library(R.utils)
# R.utils::setOption("clusterProfiler.download.method","auto")

m <- f[f$module == 16,]
# m$name <- as.character(m)
test1 = bitr(m$name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ekk <- enrichKEGG(gene = test1$ENTREZID, organism  = 'hsa', 
                  pAdjustMethod = 'BH', pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2, keyType = 'kegg') 
ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(ekk, file="./6 富集分析/primary-kegg16.csv")
pdf("./6 富集分析/primary-kegg16.pdf", height = 6, width = 6)
dotplot(ekk, font.size = 10)
dev.off()































setwd("G:/TET/网络重建")


# GO富集

library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(openxlsx)

f <- read.xlsx("./8 距离/seed_module_vertex.xlsx", colNames=T, sep='\t', sheet = "3module_all")

go <- enrichGO(f$geneSymbol, OrgDb = org.Hs.eg.db, ont='ALL', 
               pAdjustMethod = 'BH', pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2, keyType = 'SYMBOL')
a <- as.data.frame(go)
write.csv(go, "./8 距离/3module_all-go.csv")

pdf("./8 距离/3module_all-go.pdf", width = 10, height = 10)
dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~.,scale="free")
dev.off()



# KEGG富集
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(stringr)

# m$name <- as.character(m)
test1 = bitr(f$geneSymbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ekk <- enrichKEGG(gene = test1$ENTREZID, organism  = 'hsa', 
                  pAdjustMethod = 'BH', pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2, keyType = 'kegg') 
ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
b <- as.data.frame(ekk)
write.csv(ekk, file="./8 距离/3module_all-kegg.csv")
pdf("./8 距离/3module_all-kegg.pdf", height = 6, width = 6)
dotplot(ekk, font.size = 10)
dev.off()



















setwd("G:/TET/网络重建")


####### GO富集

library(org.Hs.eg.db) #人类注释数据库
library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)
library(enrichplot)

f <- read.table("./5 模块划分/primary_node_Module_5.txt", header=T, sep='\t')
m <- f[,1]

go <- enrichGO(m, OrgDb = org.Hs.eg.db, ont='ALL', 
               pAdjustMethod = 'BH', pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2, keyType = 'SYMBOL')
write.csv(go, "./6 富集分析/primary5module-go.csv")

pdf("./6 富集分析/primary5module-go.pdf", width = 10, height = 10)
dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~.,scale="free")
dev.off()



####### KEGG富集
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(stringr)
# library(R.utils)
# R.utils::setOption("clusterProfiler.download.method","auto")

m <- f
# m$name <- as.character(m)
test1 = bitr(m$name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ekk <- enrichKEGG(gene = test1$ENTREZID, organism  = 'hsa', 
                  pAdjustMethod = 'BH', pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2, keyType = 'kegg') 
ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(ekk, file="./6 富集分析/primary5module-kegg.csv")
pdf("./6 富集分析/primary5module-kegg.pdf", height = 6, width = 6)
dotplot(ekk, font.size = 10)
dev.off()




