library(actor)
library(readr)
library(ggplot2)
library(pheatmap)
set.seed(1992)



#Sim1
#Load Isoform Data
isoDataOrig1=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/simData/sim_NoIntersection_Counts.txt",delim="\t")
isoDataOrig2=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/simData/sim_BrainMix_Counts.txt",delim="\t")
isoDataOrig=rbind(isoDataOrig1,isoDataOrig2)
simGeneList1=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/sim_NoIntersection.txt",delim="\t")
simGeneList2=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/sim_BrainMix.txt",delim="\t")
colnames(simGeneList2)=c("gene","Tissue","SimTissue")
simGeneList=rbind(simGeneList1,simGeneList2[,-3])

tissueList=c("Brain_Cerebellum", "Liver","Pancreas","Nerve_Tibial","Brain_Spinal_cord_(cervical_c_1)")

#Run Model for Sim 1
m1=actor(isoDataOrig,reduceTissue = TRUE,tissueList=tissueList)
fit1=actorFit(m1)

simGeneListFinal=simGeneList[match(rownames(fit1$likelihoodSum),simGeneList$gene),]
colnames(simGeneListFinal)[2]="Sim_Tissue"

  p1=plot(fit1,plotType = "TissueMem",addAnnot=TRUE,annotMat_Col = simGeneListFinal[,-1])
  ggsave(p1,file="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/sim1HeatGG.png",width=15,height=10,dpi=400)

fit1$dirEta

plot(fit1,plotType = "ClassMem",addAnnot = TRUE)

plot(fit1,plotType = "GeneClass")


write_delim(simGeneListFinal,'C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/simGeneListModel.txt',delim="\t")


simTissue=apply(fit1$classLambda,1,which.max)
table(simTissue,simGeneListFinal$Sim_Tissue)

#Sim 2 - Spleen/Liver


livSpleenCts=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/simData/liverSpleenGeneSet.txt",delim="\t")
livSpleenTissue=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/liverSpleenGenesSet_GeneSim.txt",delim="\t")

simGeneListLiv=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/liverSpleenGenesSet_GeneSim.txt",delim="\t")
set.seed(1993)
m2=actor(livSpleenCts,reduceTissue = TRUE)
fit2=actorFit(m2)
simGeneListLivFinal=simGeneListLiv[simGeneListLiv$gene_id%in%rownames(fit2$likelihoodSum),]

p2=plot(fit2,plotType="TissueMem",addAnnot = TRUE,annotMat_Col = simGeneListLivFinal[,-1])
ggsave(p2,file="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/sim2HeatGG.png",width=15,height=10,dpi=400)

fit2$dirEta

write_delim(simGeneListLivFinal,"C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/liverSpleenGeneListModel.txt",delim="\t")

simTissueLiv=apply(fit2$classLambda,1,which.max)
table(simTissueLiv,simGeneListLivFinal$Tissue)


#Motor Neuron Set
#FASTQ files available from https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP067645

isoDataOrig=read_delim("C:/Sean/UNC/Mike Love/GTEx/data/newData/motorNeurons/motorNeuronsCts.txt",delim="\t")
condition=read_delim("C:/Sean/UNC/Mike Love/GTEx/data/newData/motorNeurons/condition.txt",delim="\t")
isoDataOrig=isoDataOrig[,c(TRUE,TRUE,condition$condition=="Control")]

set.seed(1994)
#Fit Neuron data model
m3=actor(isoDataOrig,reduceTissue = TRUE)
fit3=actorFit(m3)

p3=plot(fit3,plotType = "TissueMem")
ggsave(p3,file="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/brainHeatGG.png",width=15,height=10,dpi=400)

fit3$dirEta

brainList=data.frame("gene"=rownames(fit3$likelihoodSum))
write_delim(brainList,"C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/motorNeuronGeneListModel.txt",delim="\t")


#Function to calculate Gene Expression
getGeneExpression2 <- function(cts, gene_id) {
  gene.cts <-rowsum(cts,  gene_id)
  return(gene.cts)
  
}

###Supplementary Figures

gp1=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000064419.9")
ggsave(filename="goodLiver.png",plot=gp1,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=7,height=8  )
gp2=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000047932.9")
ggsave(filename="goodLiver2.png",plot=gp2,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=7,height=8   )
gp3=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000077721.11")
ggsave(filename="goodLiver3.png",plot=gp3,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=7,height=8 )


gp4=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000101474.7")
ggsave(filename="goodHeart.png",plot=gp4,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/" ,dpi=400,width=7,height=8 )
gp5=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000106244.8")
ggsave(filename="goodHeart2.png",plot=gp5,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/" ,dpi=400,width=7,height=8 )
gp6=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000106261.12")
ggsave(filename="goodHeart3.png",plot=gp6,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/" ,dpi=400,width=7,height=8 )



gp7=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000112584.9")
ggsave(filename="goodMusc.png",plot=gp7,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=7,height=8   )
gp8=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000114999.7")
ggsave(filename="goodMusc2.png",plot=gp8,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=7,height=8   )
gp9=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000117475.9")
ggsave(filename="goodMusc3.png",plot=gp9,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=7,height=8   )


gp10=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000108797.7")
ggsave(filename="goodAmyg.png",plot=gp10,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=7,height=8   )


gp11=makeGenePlots(fit3,isoDataOrig,geneList="ENSG00000122484.8")
ggsave(filename="goodCereb.png",plot=gp11,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=7,height=8   )




medData=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/referencePanel/medioid/fullTissueMedioid.txt",delim="\t")
medData=medData[medData$gene_id%in%rownames(fit3$likelihoodSum),c("gene_id","feature_id",colnames(fit3$likelihoodSum))]
medGeneExp=data.frame(getGeneExpression2(medData[,-(1:2)],medData$gene_id))
medGeneExpOrd=medGeneExp[rownames(fit3$likelihoodSum),]


pMed=pheatmap(t(log(medGeneExpOrd[p3$tree_col$order,p3$tree_row$order]+1)),cluster_cols = FALSE,cluster_rows = FALSE,show_colnames = FALSE,breaks=seq(0,8,0.08))
ggsave(filename="geneExpHeat.png",plot=pMed,path="C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/supp/",dpi=400,width=8,height=6)





