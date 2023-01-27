
library(igraph)
library(GOSemSim)
library(org.Hs.eg.db)
library(openxlsx)

net <- read.xlsx("./ppi_pearson_600_0.3.xlsx", colNames = T,rowNames = F)

net_edges <- net[,1:2]

# Create the network
g = graph_from_data_frame(net_edges, directed=FALSE)
net_nodes = as.data.frame(rbind(as.matrix(net_edges[,1]),as.matrix(net_edges[,2])))
net_nodes = unique(net_nodes)

#Calculate the edge betweenness
net_bet <- edge_betweenness(g, e = E(g), directed = F)
net_bet <- cbind(net_edges, net_bet)
colnames(net_bet)=c("name1","name2","bet")
data1 <- merge(net, net_bet)


# Calculate functional similarity
hsGO <- godata('org.Hs.eg.db', keytype = "SYMBOL",
               ont="BP", computeIC=FALSE)
funcsim <- mgeneSim(net_nodes[,1], semData=hsGO, measure="Wang",
                    combine="BMA", verbose=FALSE)
upsim <- upper.tri(funcsim)
mergesim <- data.frame(row = rownames(funcsim)[row(funcsim)[upsim]],
                       column = rownames(funcsim)[col(funcsim)[upsim]],
                       cor =(funcsim)[upsim])
temp = mergesim
temp[,1] = mergesim[,2]
temp[,2] = mergesim[,1]
finalsim <- rbind(mergesim, temp)
colnames(finalsim) = c("name1", "name2", "cor")
colnames(net_edges) = c("name1", "name2")
GO <- merge(finalsim, net_edges,by=c("name1", "name2"), sort=F)
 #TFC of edges without functional similarity to be zero
nasim=dplyr::anti_join(net_edges[,1:2], GO[,1:2],by = c("name1", "name2"))
nasim$cor=0
GO=rbind(GO,nasim)

data2 <- merge(data1, GO)

setwd("./wf-frank2019-ne.PCA-2a307ce")
library(devtools)
devtools::load_all()
library(ne.PCA)

## core network
ecount<-nrow(net_edges)

NET <- data2[,c(1,2,7,8)]
colnames(NET) <- c("name1", "name2", "bet", "cor")

scaleA<-NET$bet-min(NET$bet)
scaleB<-max(NET$bet)-min(NET$bet)
#Edge betweenness normalization
NET$scaleTN <- scaleA/scaleB
TF<-NET$scaleTN+NET$cor
NET$TFC  <- 100*TF/(2-TF)
TFC=NET[,1:2]
TFC$TFCs=NET$TFC

data3 <- merge(data2, TFC)

#Export TFC score
write.xlsx(data3, 'all_interaction_TFC.xlsx', sep="\t", quote=F, rowNames = F)


