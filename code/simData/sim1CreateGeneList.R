###Obtain genes for simulation


library(readr)
library(dplyr)
library(magrittr)
library(dplyr)
library(pheatmap)
library(matrixStats)
minAlpha=0.00005
set.seed(1992)

#Function to subset gene names
strp=function(x){
  substr(x,1,15)
}
#Function to calculate Gene Expression
getGeneExpression <- function(cts, gene_id) {
  gene.cts <-rowsum(cts,  gene_id)
  total.cts <- gene.cts[match(gene_id, rownames(gene.cts)),]
  return(total.cts)
}
##Function to calculate the parameters of the BB distribution. estvec is the proportion estimates and precVec is the precision
calculateDirichletAlpha=function(estVec,precVec  ){
  cn=colnames(estVec)
  estTemp=precVec[match(estVec$gene_id,precVec$gene_id),]
  dTempAlpha1=data.frame("gene_id"=estVec$gene_id,"feature_id"=estVec$feature_id,  estTemp[,-1]*(estVec[,-(1:2)]))
  colnames(dTempAlpha1)=cn
  return(dTempAlpha1)
}




#Function to calculate Gene Expression
getGeneExpression2 <- function(cts, gene_id) {
  gene.cts <-rowsum(cts,  gene_id)
  return(gene.cts)

}


load("C:/Sean/UNC/Mike Love/GTEx/rPackages/director/R/sysdata.rda")
tissueAlpha=alphaDirSet
numIso=table(tissueAlpha$gene_id)
keepGeneTissue=names(numIso)[numIso<=5]
tissueList=c("gene_id","feature_id","Brain_Cerebellum", "Liver","Pancreas","Nerve_Tibial","Brain_Spinal_cord_(cervical_c_1)")
tissueAlphaRed=tissueAlpha[tissueAlpha$gene_id%in%keepGeneTissue,tissueList]
tissueList=c("gene_id","feature_id","Brain_Cerebellum", "Liver","Pancreas","Nerve_Tibial","Brain_Spinal_Cord")
colnames(tissueAlphaRed)=tissueList
tissueAlphaRed=tissueAlphaRed[apply(!is.na(tissueAlphaRed[,-(1:2)]),1,sum)>0,]
#tissueAlphaRed[tissueAlphaRed<minAlpha]=NA
minAlphaFunc=function(vec){
  if(all(is.na(vec))){
    return(vec)
  }
  vec[is.na(vec)]=minAlpha
  if(all(vec<minAlpha)){
    return(rep(NA,length(vec)))
  }
  return(vec)
}

tissueAlphaRed%>%group_by(gene_id)%>%mutate_at(vars(-gene_id,-feature_id),minAlphaFunc)->tissueAlphaRed



numIsNA=data.frame("gene_id"=tissueAlphaRed$gene_id,"numNA"=rowSums(!is.na(tissueAlphaRed[,-(1:2)])))
numIsNA%>%group_by(gene_id)%>%summarize(maxNA=max(numNA))->maxNAVal

tissueAlphaRedFull=tissueAlphaRed[tissueAlphaRed$gene_id%in% (maxNAVal$gene_id[maxNAVal$maxNA>2]),]
sumAlpha=getGeneExpression(tissueAlphaRedFull[,-(1:2)],tissueAlphaRedFull$gene_id)
tissueEst=data.frame("gene_id"=tissueAlphaRedFull$gene_id,"feature_id"=tissueAlphaRedFull$feature_id,tissueAlphaRedFull[,-(1:2)]/sumAlpha )

sumAlpha[sumAlpha>50]=50
tissueAlphaRedFull=data.frame("gene_id"=tissueAlphaRedFull$gene_id,"feature_id"=tissueAlphaRedFull$feature_id,tissueEst[,-(1:2)]*sumAlpha)
tissuePrec=sumAlpha[!duplicated(tissueAlphaRedFull$gene_id),]


hellingerFunc=function(vec){
  return(sqrt(2*sum(vec)))
}


numTissue=ncol(tissueEst)-2
minHellingerMatrix=matrix(NA,ncol=numTissue,nrow=nrow(tissuePrec))
rownames(minHellingerMatrix)=rownames(tissuePrec)
colnames(minHellingerMatrix)=colnames(tissuePrec)


for(i in 1:numTissue){
  tissueIndex=i+2
  diffTemp=data.frame("gene_id"=tissueEst$gene_id,(sqrt(tissueEst[,-c(1,2,tissueIndex)])-sqrt(tissueEst[,tissueIndex]))^2)
  diffTemp%>%group_by(gene_id)%>%summarize_all(hellingerFunc)->hellingerDiff
  minHellinger=apply(hellingerDiff[,-1],1,min,na.rm=TRUE)
  minHellinger[is.infinite(minHellinger)]=0
  minHellingerMatrix[,i]=minHellinger
  if(i==1){
    brainDiff=unlist(hellingerDiff[,ncol(hellingerDiff)])
    brainMin=apply(hellingerDiff[,2:4],1,min,na.rm=TRUE)
  }


}
brainMat=data.frame("BrainDiff"=brainDiff,"nextMin"=brainMin)
rownames(brainMat)=rownames(tissuePrec)
brainMat$precDiff=tissuePrec[,1]-tissuePrec[,5]

maxDiff=apply(minHellingerMatrix,1,min)
whichMaxDiff=apply(minHellingerMatrix,1,which.max)

brainGenes=(maxDiff)[maxDiff>0.1&whichMaxDiff==1]
liverGenes=(maxDiff)[maxDiff>0.1&whichMaxDiff==2]
pancreasGenes=(maxDiff)[maxDiff>0.1&whichMaxDiff==3]
nerveGenes=(maxDiff)[maxDiff>0.1&whichMaxDiff==4]
brain2Genes=(maxDiff)[maxDiff>0.1&whichMaxDiff==5]

brainGenesFinal=names(brainGenes[order(brainGenes,decreasing=TRUE)[1:100]])
liverGenesFinal=names(liverGenes[order(liverGenes,decreasing=TRUE)[1:100]])
pancreasGenesFinal=names(pancreasGenes[order(pancreasGenes,decreasing=TRUE)[1:25]])
nerveGenesFinal=names(nerveGenes[order(nerveGenes,decreasing=TRUE)[1:50]])
brain2GenesFinal=names(brain2Genes[order(brain2Genes,decreasing=TRUE)[1:50]])


finalGeneList=data.frame("gene"=c(brainGenesFinal,liverGenesFinal,pancreasGenesFinal,nerveGenesFinal,brain2GenesFinal),
                         "Tissue"=c(rep("Brain_Cerebellum",length(brainGenesFinal)),rep("Liver",length(liverGenesFinal)),rep("Pancreas",length(pancreasGenesFinal)),
                               rep("Nerve_Tibial",length(nerveGenesFinal)),rep("Brain_Spinal_Cord",length(brain2GenesFinal) )  ))

write_delim(finalGeneList,"C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/data/geneLists/sim_NoIntersection.txt",delim="\t")


#Reduce to tissue with a small precision difference
brainMatRed=brainMat[abs(brainMat$precDiff)<5&(!is.na(brainMat$precDiff)),]
brainMatRedFilt=brainMatRed[which(brainMatRed$BrainDiff<0.04&brainMatRed$nextMin>0.08),]
brainMatRedFilt=brainMatRedFilt[order(brainMatRedFilt$BrainDiff,decreasing=FALSE)[1:50],]
brainMixDF=data.frame("gene"=rownames(brainMatRedFilt),"TissueGroup"=rep("MixBrain",nrow(brainMatRedFilt)),"Tissue"=rbinom(nrow(brainMatRedFilt),1,0.5  ))

write_delim(brainMixDF,"C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/data/geneLists/sim_BrainMix.txt",delim="\t")

