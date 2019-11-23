library(readr)
library(magrittr)



sampDat=read_tsv("/proj/milovelab/mccabe/proj/GTEx/data/GTEx_v7_Annotations_SampleAttributesDS.txt")
colN=read_table("/proj/milovelab/mccabe/proj/GTEx/data/colNames.txt")
colN2=unlist(colN[-(1:2),1])
matchCol=sampDat$SAMPID%in%colN2
strp=function(x) substr(x,1,15)



sampDatFull=sampDat[matchCol,]


smallSampleTissue=c("Cervix_Endocervix","Cervix_Ectocervix","Fallopian_Tube","Bladder","Kidney_Cortex")

sampDatRed=sampDatFull[!(sampDatFull$SMTSD%in%smallSampleTissue),]
dim(sampDatRed)

gtex=read_tsv("/proj/milovelab/mccabe/proj/GTEx/data/gtex.txt")
dim(gtex)
gtex=gtex[,c("transcript_id","gene_id",sampDatRed$SAMPID) ]

dim(gtex)


cs=apply(gtex[,-(1:2)],2,sum)
cs=data.frame("sample_id"=colnames(gtex)[-(1:2)],"libSize"=cs)

write.table(cs,"/proj/milovelab/mccabe/proj/GTEx/data/GTExLibrarySize.txt",sep="\t",quote=FALSE,row.names = FALSE)
