library(readr)
library(GenomicFeatures)


alphaDirSet=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/referencePanel/fullTissueAlphaNoIsoFilt_MLE.txt",delim="\t")
ensdb=loadDb("C:/Sean/UNC/Mike Love/GTEx/RCode/GTEx/subset/gencode_v19.sqlite")
keyEns=AnnotationDbi::select(ensdb,keys=unique(alphaDirSet$gene_id),keytype = "GENEID",columns=c("GENEID","TXNAME"))
numIso=table(keyEns$GENEID)

keepGenes=names(numIso)[numIso<6]

alphaDirSetRed=alphaDirSet[alphaDirSet$gene_id%in%keepGenes,]

write_delim(alphaDirSetRed,"C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/referencePanel/finalFullTissueAlpha.txt",delim="\t")
