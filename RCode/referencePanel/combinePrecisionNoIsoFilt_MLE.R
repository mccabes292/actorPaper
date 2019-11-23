library(readr)
library(magrittr)
library(DRIMSeq)
library(dplyr)
library(GenomicFeatures)

sampDat=read_tsv("/proj/milovelab/mccabe/proj/GTEx/data/GTEx_v7_Annotations_SampleAttributesDS.txt")
colN=read_table("/proj/milovelab/mccabe/proj/GTEx/data/colNames.txt")
colN2=unlist(colN[-(1:2),1])
matchCol=sampDat$SAMPID%in%colN2
strp=function(x) substr(x,1,15)

sampDatFull=sampDat[matchCol,]
tab1=table(sampDatFull$SMTSD)
list1=names(tab1)


makeTissueName=function(index){
  tissueName=paste(strsplit(gsub("-"," ",list1[index]),"\\s+")[[1]],collapse="_")
  return(tissueName)
}

n1=length(list1)

tissueNames=sapply(1:n1,makeTissueName)
smallSampleTissue=c("Cervix_Endocervix","Cervix_Ectocervix","Fallopian_Tube","Bladder","Kidney_Cortex","Whole_Blood","Testis")
tissueNames=tissueNames[!(tissueNames%in%smallSampleTissue)]

numDF=data.frame(tab1)
colnames(numDF)=c("tissue","samples")

write.table("/proj/milovelab/mccabe/proj/GTEx/data/tissueTable.txt",quote=F,row.names=F,sep="\t")


mergeAlpha=read_delim(paste("/proj/milovelab/mccabe/proj/GTEx/data/precisionNoIsoFilt_MLE/", tissueNames[1],"AlphaNoIsoFilt.txt",sep=""),delim="\t")
colnames(mergeAlpha)=c("gene_id","feature_id",tissueNames[1])

for(i in 2:length(tissueNames)){

  tempAlpha=read_delim(paste("/proj/milovelab/mccabe/proj/GTEx/data/precisionNoIsoFilt_MLE/", tissueNames[i],"AlphaNoIsoFilt.txt",sep=""),delim="\t")
  tempCol1=colnames(mergeAlpha)


  mergeAlpha=merge(mergeAlpha,tempAlpha,by=c("gene_id","feature_id"),all=TRUE)


  colnames(mergeAlpha)=c(tempCol1,tissueNames[i])
}

ensdb=loadDb("/proj/milovelab/mccabe/proj/GTEx/data/gencodeV19/gencode.v19.annotation.sqlite")
keyEns=AnnotationDbi::select(ensdb,keys=unique(mergeAlpha$gene_id),keytype = "GENEID",columns=c("GENEID","TXNAME"))
numIso=table(keyEns$GENEID)
singleIsoGenes=numIso[numIso==1]
mergeAlpha=mergeAlpha[!(mergeAlpha$gene_id%in%names(singleIsoGenes)),]


write_delim(mergeAlpha,"/proj/milovelab/mccabe/proj/GTEx/data/precisionNoIsoFilt_MLE/fullTissueAlphaNoIsoFilt_MLE.txt",delim="\t")

