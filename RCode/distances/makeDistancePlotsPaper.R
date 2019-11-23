library(readr)
library(stringr)
library(ggplot2)
library(scales)
library(cowplot)
convertId=function(x){
  len=str_length(x)
  if(substr(x,len-1,len  )==".1"){
    return(str_replace_all(substr(x,1,len-2),"\\.","-"  ) )
  }else{
    return(str_replace_all(x,"\\.","-"  ) )
  }
}

makeFullDistancePlots=function(distMat,numGenes,gtexDemoOrd,fileName,filePath,pointSize=1,plotWidth=18,addLegend=TRUE,makeIndivPlots=TRUE,fontSize=21){
  if(all(distMat$GTExID==gtexDemoOrd$sample_id)){
    message("Samples Match")
  }else{
    warning("Error: Samples do not match")
  }
  
  n <- length(levels(factor(gtexDemoOrd$SMTSD))) # number of colors
  
  cols <- hue_pal(h = c(0, 360) + 15,
                  c = 100, l = 65,
                  h.start = 0, direction = 1)(n)[order(sample(1:n, n))]
  if(makeIndivPlots==TRUE){
    for(i in 2:ncol(distMat)){
      df=data.frame("x"=1:nrow(distMat), "dist"=unlist(distMat[,i])/numGenes,"Tissue"=factor(gtexDemoOrd$SMTSD))
      p1=ggplot(df,aes(x=x,y=dist,col=Tissue))+geom_point(shape=pointSize)+xlab(colnames(distMat)[i])+ylab("Hellinger Distance")+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x=element_blank(),text=element_text(size=21),text=element_text(size=fontSize))+scale_color_manual(values = cols)
      if(addLegend==TRUE){
        p1=p1+theme(legend.position = "none")
      }
      ggsave(paste(filePath,fileName,i,".png",sep=""),plot=p1,width = plotWidth,dpi=400)
    }
  }
  df=data.frame("x"=1:nrow(distMat),"avgDist"=rowMeans( distMat[,-1]/numGenes  ),"Tissue"=factor(gtexDemoOrd$SMTSD))
  p1=ggplot(df,aes(x=x,y=avgDist,col=Tissue))+geom_point(shape=pointSize)+xlab("")+ylab("Hellinger Distance")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x=element_blank(),text=element_text(size=fontSize))+scale_color_manual(values = cols)
  
  if(addLegend==FALSE){
  
    p1=p1+theme(legend.position = "none")
  }
  ggsave(paste(filePath,fileName,"Average.png",sep=""),plot=p1,width = plotWidth,dpi=400)
  
  
}






#Update demo pull from longleaf
gtexDemo=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/GTEx/fullDemo_GeneListDistancesLargeTissue.txt",delim="\t")
gtexDemoOrd=gtexDemo[order(gtexDemo$SMTSD,decreasing=FALSE),]
tissueList=c("Brain - Cerebellar Hemisphere", "Liver","Pancreas","Nerve - Tibial","Brain - Spinal cord (cervical c-1)")
gtexDemoOrdSim=gtexDemoOrd[gtexDemoOrd$SMTSD%in%tissueList,]
#Sim 1 Distances
sim1Dist=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/distances/outDist/simDistancesLargeTissuePaper.txt",delim="\t")
sim1Dist$GTExID=str_replace_all(sim1Dist$GTExID,"\\.","-")
sim1GeneList=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/simGeneListModel.txt",delim="\t")
numGenes=nrow(sim1GeneList)
sim1DistOrd=sim1Dist[match(gtexDemoOrdSim$sample_id,sim1Dist$GTExID),]
set.seed(1990)

makeFullDistancePlots(distMat=sim1DistOrd,numGenes=numGenes,gtexDemoOrd = gtexDemoOrdSim,fileName= "sim1Dist",
                      filePath =  "C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/distances/sim1/",makeIndivPlots = FALSE,
                      pointSize=16,plotWidth=10,fontSize=21)


#Liver Spleen Distances
livSpleenDist=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/distances/outDist/mixDistancesLargeTissuePaper.txt",delim="\t")
livSpleenDist$GTExID=sapply(livSpleenDist$GTExID,convertId)
livSpleenGeneList=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/liverSpleenGeneListModel.txt",delim="\t")
numGenesLivSpleen=nrow(livSpleenGeneList)
livSpleenDistOrd=livSpleenDist[match(gtexDemoOrd$sample_id,livSpleenDist$GTExID),]
set.seed(1991)
makeFullDistancePlots(distMat=livSpleenDistOrd,numGenes=numGenesLivSpleen,gtexDemoOrd = gtexDemoOrd,fileName= "livSpleenDist",
                      filePath =  "C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/distances/livSpleen/",
                      pointSize=1,plotWidth=10,makeIndivPlots=FALSE,addLegend=FALSE,fontSize=21)
set.seed(1991)
makeFullDistancePlots(distMat=livSpleenDistOrd,numGenes=numGenesLivSpleen,gtexDemoOrd = gtexDemoOrd,fileName= "livSpleenDistLegend",
                      filePath =  "C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/distances/livSpleen/",
                      pointSize=1,plotWidth=20,makeIndivPlots=FALSE,addLegend=TRUE,fontSize=10)




#Brain Distances
brainDist=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/distances/outDist/brainDistancesLargeTissuePaper.txt",delim="\t")
brainGeneList=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/motorNeuronGeneListModel.txt",delim="\t")
numGenesBrain=nrow(brainGeneList)
brainDist$GTExID=sapply(brainDist$GTExID,convertId)
brainDistOrd=brainDist[match(gtexDemoOrd$sample_id,brainDist$GTExID),]
set.seed(1992)
makeFullDistancePlots(distMat=brainDistOrd,numGenes=numGenesBrain,gtexDemoOrd = gtexDemoOrd,fileName= "brainDist",
                      filePath =  "C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/distances/brain/",
                      pointSize=1,plotWidth=10,makeIndivPlots=FALSE,addLegend=FALSE,fontSize=21)
set.seed(1992)
makeFullDistancePlots(distMat=brainDistOrd,numGenes=numGenesBrain,gtexDemoOrd = gtexDemoOrd,fileName= "brainDistLegend",
                      filePath =  "C:/Sean/UNC/Mike Love/GTEx/actorPaper/results/distances/brain/",
                      pointSize=1,plotWidth=20,makeIndivPlots=FALSE,addLegend=TRUE,fontSize=10)
