# The construction of PPIN based on RNs.

library("openxlsx")
library("e1071")

svmrfeFeatureRanking = function(x,y){
  
  #Checking for the variables
  stopifnot(!is.null(x) == TRUE, !is.null(y) == TRUE)
  
  n = ncol(x)
  survivingFeaturesIndexes = seq_len(n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  
  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,
                   scale=FALSE, type="C-classification", kernel="linear" )
    
    #compute the weight vector
    w = t(svmModel$coefs)%*%svmModel$SV
    
    #compute ranking criteria
    rankingCriteria = w * w
    
    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)$ix
    
    #update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1
    cat(paste0(rankedFeatureIndex,"\r"))
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
  }
  
  return (featureRankedList)
}

data <- read.xlsx("proteome.xlsx", sheet = "rec")

y <- data[1,seq(3,ncol(data))]
data<-data[-1,]
proteinid <- data[,1]
geneid <- data[,2]
data <- data[,c(-1,-2)]

rownames(data) <- 1:nrow(data)
ranki <- svmrfeFeatureRanking(x=t(data),y=as.character(y))
ranki
n = length(ranki)
li = data.frame()
for (i in 1:n){
	li[i,1] <- (n + 1 - ranki[i]) / n

}
data_Rs <- data.frame(geneid,proteinid,ranki,li)
colnames(data_Rs)[4] <- 'Rs'

data_network <- read.csv("ppi_pearson_600_0.3.txt default node.csv",header = T,stringsAsFactors = F)
colnames(data_network)[8] <- 'proteinid'
data_network_new <- data_network[,c(8,1,5)]
data_all <- merge(data_Rs,data_network_new,by = "proteinid",all = T)
for (i in 1:nrow(data_all)){
  data_all[i,7] <- data_all[i,4] * data_all[i,6] / data_all[i,5]
}
colnames(data_all)[7] <- 'RN'
data_all = data_all[order(data_all[,7],decreasing=T),]

for (i in 1:nrow(data_all)){
  data_all[i,8] <- i
}
colnames(data_all)[8] <- 'rank_by_RN'


