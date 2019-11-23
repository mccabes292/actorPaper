library(readr)
library(magrittr)
library(DirichletMultinomial)
args=commandArgs(TRUE)
index=as.numeric(args[1])

calcDMMLE=function(gene){
  print(gene)
  ctsTemp=gtex[gtex$gene_id==gene,]
  cs=colSums(ctsTemp[,-(1:2)])
  ctsTemp=ctsTemp[,c(TRUE,TRUE,cs!=0)]
  fit=DirichletMultinomial::dmn(t(ctsTemp[,-(1:2)]),k=1,verbose=TRUE )
  est=fit@fit$Estimate
  rownames(est)=(ctsTemp$feature_id)
  return(est)

}

sampDat=read_tsv("/proj/milovelab/mccabe/proj/GTEx/data/GTEx_v7_Annotations_SampleAttributesDS.txt")
colN=read_table("/proj/milovelab/mccabe/proj/GTEx/data/colNames.txt")
colN2=unlist(colN[-(1:2),1])
matchCol=sampDat$SAMPID%in%colN2
strp=function(x) substr(x,1,15)

sampDatFull=sampDat[matchCol,]
list1=names(table(sampDatFull$SMTSD))
tissueType=list1[index]
sampDatRed=sampDatFull[sampDatFull$SMTSD==tissueType,]
dim(sampDatRed)
print(tissueType)

tissueName=paste(strsplit(gsub("-"," ",tissueType),"\\s+")[[1]],collapse="_")
print(tissueName)


#isoList=read.table("/proj/milovelab/mccabe/proj/GTEx/data/GTExSuffExprIsos.txt")
#isoList2=unlist(isoList)

gtex=read_tsv("/proj/milovelab/mccabe/proj/GTEx/data/GTExScaleTPM.txt")
dim(gtex)
gtex=gtex[,c("transcript_id","gene_id",sampDatRed$SAMPID) ]

dim(gtex)
intercept=data.frame("sample_id"=sampDatRed$SAMPID, "int"=rep(1,nrow(sampDatRed)))
#Find precision estimates using DRIMSeq
cn=colnames(gtex)
cn[1]="feature_id"
gtex=data.frame(gtex)
colnames(gtex)=cn
gtex=gtex[,c(1,2,3:ncol(gtex))]
dim(gtex)
rs=apply(gtex[,-(1:2)],1,sum)
gtex=gtex[rs>0,]
dim(gtex)
gene.cts <-rowsum(gtex[,-(1:2)],  gtex$gene_id)
gene.cts.sum=apply(gene.cts>10,1,sum)
gtex=gtex[gtex$gene_id%in%names(gene.cts.sum[gene.cts.sum>0.7*nrow(intercept)]),]
numGenes=table(gtex$gene_id)
keepGenes=names(numGenes)[numGenes>1]
gtex=gtex[gtex$gene_id%in%keepGenes,]
geneList=unique(gtex$gene_id)
precVal=lapply(geneList,calcDMMLE)
precDF=data.frame("gene_id"=gtex$gene_id,"feature_id"=gtex$feature_id,"prec"=unlist(precVal))

write_delim(precDF,paste("/proj/milovelab/mccabe/proj/GTEx/data/precisionNoIsoFilt_MLE/",tissueName,"AlphaNoIsoFilt.txt",sep=""),delim="\t")
