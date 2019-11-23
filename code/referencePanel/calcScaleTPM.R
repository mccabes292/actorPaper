library(readr)
library(magrittr)



sampDat=read_tsv("/proj/milovelab/mccabe/proj/GTEx/data/GTEx_v7_Annotations_SampleAttributesDS.txt")
colN=read_table("/proj/milovelab/mccabe/proj/GTEx/data/colNames.txt")
colN2=unlist(colN[-(1:2),1])
matchCol=sampDat$SAMPID%in%colN2
strp=function(x) substr(x,1,15)



sampDatFull=sampDat[matchCol,]


smallSampleTissue=c("Cervix - Endocervix","Cervix - Ectocervix","Fallopian Tube","Bladder","Kidney - Cortex")

sampDatRed=sampDatFull[!(sampDatFull$SMTSD%in%smallSampleTissue),]
dim(sampDatRed)

gtex=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz",delim="\t")
dim(gtex)
gtex=gtex[,c("transcript_id","gene_id",sampDatRed$SAMPID) ]

dim(gtex)

libSize=read_delim("/proj/milovelab/mccabe/proj/GTEx/data/GTExLibrarySize.txt",delim="\t")
all(libSize$sample_id==colnames(gtex)[-(1:2)])

scale=1e6

scaleTPM=data.frame(gtex$transcript_id,gtex$gene_id,t(t(gtex[-(1:2)])*libSize/1e6 ))
colnames(scaleTPM)=colnames(gtex)


write.table(scaleTPM,"/proj/milovelab/mccabe/proj/GTEx/data/GTExScaleTPM.txt",sep="\t",quote=FALSE,row.names = FALSE)
