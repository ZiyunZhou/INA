# @ author: Ziyun Zhou
# @ Email: zhouziyun1900@hotmail.com

##############################################################################
#
#   @> Step 0. Import R packages and data in need
#
##############################################################################

setwd("./INA/example")
library(INA)
library(openxlsx)
# load ppi database
data("ppi_database")

##############################################################################
#
#   @> Step 1. Generate PPI netwrok file.
#
##############################################################################

expression <- read.xlsx("./Data/proteome.xlsx", colNames = T)
group <- read.table("./Data/group.txt", header = T)
expression <- expr_filter(expression, filter = 0.5)
proteome_net <- proteome_net_construction(expression, ppi_database, pcc_cut = 0.3)
DEPs_protein <- DEPs(expression, group, fc_cutoff = 1.5, p_cutoff = 0.05)
DEPs_net <- DEPs_net_construction(DEPs_protein, proteome_net)
ranki <- ranking(expression, group)
# Import the whole network into cytoscape to calculate the shortest path and degree, and read the result into R
topology_data <- read.csv("./Data/Sheet 1 default node.csv", header = T, stringsAsFactors = F)
RN_score <- RNs(topology_data, ranki)
RNs_net <- RNs_net_construction(RN_score, proteome_net, top_cut = 50)


##############################################################################
#
#   @> Step 2. Network module identification and annotations.
#
##############################################################################

DEPs_module <- module(DEPs_net, module_cut = 9)


##############################################################################
#
#   @> Step 3. Shortest path-based prioriting.
#
##############################################################################

# example: select module 3 to calculate psv2n scores
select_vertex <- read.xlsx('./Data/vertex.xlsx', colNames = T)
seed2DEPs <- Shortest_path(select_vertex, proteome_net)
