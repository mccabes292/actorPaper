library(readr)
gtex=read_tsv("/proj/milovelab/mccabe/proj/GTEx/data/gtex.txt")
cn=colnames(gtex)
write.table(cn,"/proj/milovelab/mccabe/proj/GTEx/data/colNames.txt",quote=F,row.names=F)
