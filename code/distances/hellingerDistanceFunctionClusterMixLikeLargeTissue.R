library(readr)
library(dplyr)
library(magrittr)
#Function to calculate Gene Expression
getGeneExpression <- function(cts, gene_id) {
  gene.cts <-rowsum(cts,  gene_id)
  total.cts <- gene.cts[match(gene_id, rownames(gene.cts)),]
  return(total.cts)
}
hellingerFunc=function(vec){
  return(sqrt(2*sum(vec))   )
}



makeDistancePlot=function(sampData){
  numSamp=ncol(sampData)-2
  gtexDataRed=gtexData[gtexData$gene_id%in%sampData$gene_id,]
  numGtexSamp=ncol(gtexDataRed)-2
  fullMergeSet=merge(sampData,gtexDataRed,by=c("gene_id","feature_id"),all.x=TRUE,all.y=TRUE)
  fullMergeSet[is.na(fullMergeSet)]=0
  fullMergeProp=data.frame(fullMergeSet[,1:2],sqrt(fullMergeSet[,-(1:2)]/getGeneExpression(fullMergeSet[,-(1:2)],fullMergeSet$gene_id)  ))
  fullMergeProp[is.na(fullMergeProp)]=0
  distanceMatrix=matrix(NA,ncol=numSamp,nrow=numGtexSamp)
  colnames(distanceMatrix)=colnames(fullMergeProp)[3:(numSamp+2)]
  rownames(distanceMatrix)=colnames(fullMergeProp)[(numSamp+3):ncol(fullMergeProp)]
  for(i in 1:numSamp){
    message(paste("Calculating Distance for Sample ",colnames(distanceMatrix)[i]))
    tempSamp=fullMergeProp[,i+2]
    tempDf=data.frame(fullMergeProp[,1:2],(tempSamp-fullMergeProp[,-(1:(numSamp+2))])^2)
    tempDf%>%group_by(gene_id)%>%summarize_at(vars(-feature_id),hellingerFunc)->outDist
    distanceMatrix[,i]=colSums(outDist[,-1])

  }

  return(distanceMatrix)
}


gtexData=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/reducedGTExFullGeneListFromModel.txt",delim="\t")
colnames(gtexData)[1]="feature_id"


#####Mix Analysis####
#Load in experimental data
sampData=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/sampData/liverSpleenGeneSet.txt",delim="\t")


mixGeneList=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/liverSpleenGeneListModel.txt",delim="\t")

sampData=sampData[sampData$gene_id%in%mixGeneList$gene_id,]

mixDistances=makeDistancePlot(sampData)
mixOut=data.frame("GTExID"=rownames(mixDistances),mixDistances)
write_delim(mixOut,"/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/outDist/mixDistancesLargeTissuePaper.txt",delim="\t")



