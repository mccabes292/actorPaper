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

gtexDemo=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/fullDemo_GeneListDistancesLargeTissue.txt",delim="\t")
tissueList=c("Brain - Cerebellar Hemisphere", "Liver","Pancreas","Nerve - Tibial","Brain - Spinal cord (cervical c-1)")
gtexDemoRed=gtexDemo[gtexDemo$SMTSD%in%tissueList,]

gtexData=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/reducedGTExFullGeneListFromModel.txt",delim="\t")
gtexData=gtexData[,c(TRUE,TRUE,colnames(gtexData)[-(1:2)]%in%gtexDemoRed$sample_id)]
colnames(gtexData)[1]="feature_id"


#####Sim Analysis####
#Load in experimental data
sampData1=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/sampData/sim_NoIntersection_Counts.txt",delim="\t")
sampData2=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/sampData/sim_BrainMix_Counts.txt",delim="\t")
sampData=rbind(sampData1,sampData2)

simGeneList=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/simGeneListModel.txt",delim="\t")

sampData=sampData[sampData$gene_id%in%simGeneList$gene,]

simDistances=makeDistancePlot(sampData)
simOut=data.frame("GTExID"=rownames(simDistances),simDistances)
write_delim(simOut,"/proj/milovelab/mccabe/proj/GTEx/data/paperReproduce/outDist/simDistancesLargeTissuePaper.txt",delim="\t")



