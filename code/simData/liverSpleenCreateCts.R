library(readr)
liver=read_delim("C:/Sean/UNC/Mike Love/GTEx/data/subset/fullLiverCts.txt",delim="\t")
colnames(liver)[1]="feature_id"
spleen=read_delim("C:/Sean/UNC/Mike Love/GTEx/data/subset/fullSpleenCts.txt",delim="\t")
colnames(spleen)[1]="feature_id"



load("C:/Sean/UNC/Mike Love/GTEx/rPackages/director/R/sysdata.rda")
#Check if isoforms are aligned
all(liver$feature_id==spleen$feature_id)
numIso=table(liver$gene_id)
fullGeneList=unique((liver$gene_id)[liver$gene_id%in%alphaDirSet$gene_id])
set.seed(1992)
geneSample=sample(fullGeneList,length(fullGeneList)/2)
isoSample=liver$feature_id[liver$gene_id%in%geneSample]
liver[liver$feature_id%in%isoSample,-(1:2)]=0
spleen[!(spleen$feature_id%in%isoSample),-(1:2)]=0
liverSamp=sample(1:(ncol(liver)-2),100)
spleenSamp=sample(1:(ncol(spleen)-2),100)
mixSet=data.frame("gene_id"=liver$gene_id,"feature_id"=liver$feature_id,data.matrix(liver[,liverSamp+2])+data.matrix(spleen[,spleenSamp+2]))
write_delim(mixSet,"C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/data/simData/liverSpleenGeneSet.txt",delim="\t")
geneDF=data.frame("gene_id"=fullGeneList,"Tissue"=rep(NA,length(fullGeneList)))
geneDF$Tissue=ifelse(geneDF$gene_id%in%geneSample,"Spleen","Liver")
write_delim(geneDF,"C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/data/geneLists/liverSpleenGenesSet_GeneSim.txt",delim="\t")
