##Simulate data from gene list
##

#Libraries to load
library(readr)
library(MCMCpack)
library(magrittr)
library(dplyr)
#Set seed for reproducible simulations
set.seed(1992)


#Functions to be used
calculateDMAlpha=function(tissueEstGeneList,tissuePrec){
  precision=tissuePrec[match(tissueEstGeneList$gene_id,tissuePrec$gene_id),]
  tissueAlpha=tissueEstGeneList[,-(1:2)]*precision[,-1]
  df=data.frame(tissueEstGeneList[,1:2],tissueAlpha )
  return(df)
}

#Function to subset gene names
strp=function(x){
  substr(x,1,15)
}



#Function taking alpha values and provides iso expression. Function is for one gene and will be
#replicated using dplyr later
isoExpr=function(alpha){

  geneExp=rnbinom(1,size=20,mu=200)
  pDir=rdirichlet(1,alpha)
  return(rmultinom(1,size=geneExp,prob=pDir))

}


#Vectorize multinomial function to determine tissue membership
#Returns list of dimension equal to the number of samples with multinomial draws for each gene
tissueMemFunc=function(nSamp,probMat){
  # tissue=NULL
  # for(i in 1:nSamp){
  #   tissue[[i]]=apply(probMat,1,function(x){rmultinom(1,1,x)})
  #   rownames(tissue[[i]])=colnames(probMat)
  # }
  # return(tissue)
  tissue=apply(probMat,1,function(x){rmultinom(1,1,x)})
  rownames(tissue)=colnames(probMat)
  return(tissue)
}
#Number of Samples to be simulated
nSamp<<-50
idList=paste("id_",1:nSamp,sep="")

#Maximum number of isoforms to be expressed
maxIso=4

#Value to be plugged in for NA Alpha Values in the DM
minAlpha=0.00005

#Read in GTEx estimated estimates and precision
geneListMat=read_delim("C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/data/geneLists/sim_NoIntersection.txt",delim="\t")


load("C:/Sean/UNC/Mike Love/GTEx/rPackages/director/R/sysdata.rda")
tissueAlpha=alphaDirSet
tissueList=c("gene_id","feature_id","Brain_Cerebellum", "Liver","Pancreas","Nerve_Tibial","Brain_Spinal_cord_(cervical_c_1)")

tissueAlphaRed=tissueAlpha[tissueAlpha$gene_id%in%geneListMat$gene,tissueList]
tissueAlphaRed[is.na(tissueAlphaRed)]=minAlpha
tissueList=c("gene_id","feature_id","Brain_Cerebellum", "Liver","Pancreas","Nerve_Tibial","Brain_Spinal_Cord")
colnames(tissueAlphaRed)=tissueList

##########End Simulation Specifications#######
set.seed(1992)
geneListOrder=geneListMat[match(tissueAlphaRed$gene_id,geneListMat$gene),]
modelMat=model.matrix(~factor(geneListOrder$Tissue,levels = tissueList[-(1:2)])-1)
colnames(modelMat)=tissueList[-(1:2)]
modelAlpha=data.frame(tissueAlphaRed[,1:2],"alpha"=matrix(rep(rowSums(tissueAlphaRed[,-(1:2)]*modelMat),nSamp ),ncol=nSamp)  )
colnames(modelAlpha)=c("gene_id","feature_id",paste("id_",1:nSamp,sep=""))
modelAlpha%>%group_by(gene_id)%>%mutate_at(vars(-gene_id,-feature_id),isoExpr)->outIso

write_delim(outIso,"C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/data/simData/sim_NoIntersection_Counts.txt",delim="\t")
