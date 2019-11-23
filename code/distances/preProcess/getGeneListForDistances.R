library(readr)

mixList=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/liverSpleenGeneListModel.txt",delim="\t")
brainList=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/motorNeuronGeneListModel.txt",delim="\t")
simList=read_delim("C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/simGeneListModel.txt",delim="\t")


geneList=data.frame("gene"=unique(c(mixList$gene_id,brainList$gene,simList$gene)))

write_delim(geneList,"C:/Sean/UNC/Mike Love/GTEx/actorPaper/data/geneLists/fullGeneListFromModel.txt",delim="\t")

