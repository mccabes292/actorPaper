library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)


#Function to subset gene names
strp=function(x){
  substr(x,1,15)
}

set.seed(1992)
sampDat=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/GTEx_v7_Annotations_SampleAttributesDS.txt",delim="\t")
colN=read_table("/proj/milovelab/mccabe/proj/GTEx/data/colNames.txt")
colN2=unlist(colN[-(1:2),1])
matchCol=sampDat$SAMPID%in%colN2
sampDatFull=sampDat[matchCol,]


sampDatFull%>%group_by(SMTSD)%>%nest()%>%mutate(n=map_dbl(data,nrow))->sampDatRed
sampDatRed$nSamp=ifelse(sampDatRed$n<50,0,pmin(sampDatRed$n,50))
sampDatRed%>%mutate(samp=map2(data,nSamp,sample_n))%>%select(SMTSD,samp)%>%unnest()->sampDatRed

outMat=data.frame("sample_id"=sampDatRed$SAMPID,"SMTS"=sampDatRed$SMTS,"SMTSD"=sampDatRed$SMTSD)

write_delim(outMat,"/proj/milovelab/mccabe/proj/GTEx/data/actorPaper/fullDemo_GeneListDistancesLargeTissue.txt",delim="\t")

geneList=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/actorPaper/fullGeneListFromModel.txt",delim="\t")



gtex=read_tsv("/proj/milovelab/mccabe/proj/GTEx/data/gtex.txt.gz")
dim(gtex)
gtex=gtex[(gtex$gene_id)%in%as.character(unlist(geneList$gene)),c("transcript_id","gene_id",sampDatRed$SAMPID) ]
dim(gtex)
write_delim(gtex,"/proj/milovelab/mccabe/proj/GTEx/data/actorPaper/reducedGTExFullGeneListFromModel.txt",delim="\t")
