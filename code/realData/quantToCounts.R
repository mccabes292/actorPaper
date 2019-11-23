library(GenomicFeatures)
library(tximport)
txdb.filename <- "/proj/milovelab/mccabe/proj/GTEx/data/gencodeV19/gencode.v19.annotation.sqlite"
txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]


sampList=read.table("SRRList.txt",header=FALSE)
files=file.path(paste("/proj/milovelab/mccabe/fastq/motorNeurons/fastq/",sampList[,1],"_quant/quant.sf",sep=""))
names(files)=sampList[,1]
txi <- tximport(files, type="salmon", txOut=TRUE,countsFromAbundance="scaledTPM")
cts <- txi$counts
dim(cts)
cts <- cts[rowSums(cts) > 0,]
dim(cts)

strp1=function(x) substr(x,1,15)
strp2=function(x) substr(x,19,19+15-1)


txdf=txdf[match(rownames(cts),txdf$TXNAME),]


countMat=data.frame("gene_id"=txdf$GENEID,"feature_id"=txdf$TXNAME,cts   )

write.table(countMat,"/proj/milovelab/mccabe/proj/GTEx/data/otherDataSets/motorNeurons/motorNeuronsCts.txt",sep="\t",row.names = F,quote=F)

