library(readr)
library(ldaAlign)
library(ggplot2)
library(usedist)
library(pheatmap)
distFunc=function(x,y){
  sqrt(2*sum((sqrt(x)-sqrt(y))^2))
}

#Load Isoform Data
isoDataOrig=read_delim("C:/Sean/UNC/Mike Love/GTEx/data/newData/motorNeurons/motorNeuronsCts.txt",delim="\t")
condition=read_delim("C:/Sean/UNC/Mike Love/GTEx/data/newData/motorNeurons/condition.txt",delim="\t")
isoDataOrig=isoDataOrig[,c(TRUE,TRUE,condition$condition=="Control")]

m1=makeDataObject(isoDataOrig,reduceTissue = TRUE)
fit1=fitModel(m1)

distMat=dist_make((fit1$tissuePhi),distFunc)
p1=pheatmap(t(fit1$tissuePhi),cluster_cols = hclust(distMat),show_colnames = FALSE)
ggsave(filename="brainHeat.png",plot=p1,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/brainPlots/",dpi=400,width=8,height=6)




p2=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000064419.9")
ggsave(filename="goodLiver.png",plot=p2,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/",dpi=400,width=7,height=8  )
p2=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000047932.9")
ggsave(filename="goodLiver2.png",plot=p2,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/",dpi=400,width=7,height=8   )
p2=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000077721.11")
ggsave(filename="goodLiver3.png",plot=p2,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/",dpi=400,width=7,height=8 )


p3=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000101474.7")
ggsave(filename="goodHeart.png",plot=p3,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/" ,dpi=400,width=7,height=8 )
p3=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000106244.8")
ggsave(filename="goodHeart2.png",plot=p3,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/" ,dpi=400,width=7,height=8 )
p3=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000106261.12")
ggsave(filename="goodHeart3.png",plot=p3,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/" ,dpi=400,width=7,height=8 )



p4=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000112584.9")
ggsave(filename="goodMusc.png",plot=p4,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/",dpi=400,width=7,height=8   )
p4=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000114999.7")
ggsave(filename="goodMusc2.png",plot=p4,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/",dpi=400,width=7,height=8   )
p4=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000117475.9")
ggsave(filename="goodMusc3.png",plot=p4,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/",dpi=400,width=7,height=8   )


p5=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000108797.7")
ggsave(filename="goodAmyg.png",plot=p5,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/",dpi=400,width=7,height=8   )


p6=makeGenePlots(fit1,isoDataOrig,geneList="ENSG00000122484.8")
ggsave(filename="goodCereb.png",plot=p6,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/",dpi=400,width=7,height=8   )


medData=read_delim("C:/Sean/UNC/Mike Love/GTEx/data/medioid/fullTissueMedioid.txt",delim="\t")
medData=medData[medData$gene_id%in%rownames(fit1$likelihoodSum),c("gene_id","feature_id",colnames(fit1$likelihoodSum))]
medGeneExp=data.frame(getGeneExpression2(medData[,-(1:2)],medData$gene_id))
medGeneExpOrd=medGeneExp[rownames(fit1$likelihoodSum),]


pMed=pheatmap(t(log(medGeneExpOrd[p1$tree_col$order,p1$tree_row$order]+1)),cluster_cols = FALSE,cluster_rows = FALSE,show_colnames = FALSE,breaks=seq(0,8,0.08))
ggsave(filename="geneExpHeat.png",plot=pMed,path="C:/Sean/UNC/Mike Love/GTEx/PaperReproduce/results/brainPlots/",dpi=400,width=8,height=6)

