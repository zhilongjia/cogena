#! /usr/bin/Rscript


#require(clena)

require(devtools)
load_all("~/clena/")

data("PD")

#annotaion
annoGMT <- "c2.chemp.v4.0.symbols.gmt"
annoGMT <- system.file("data", annoGMT, package="clena")

#annoGMT <- "..../c2.cp.biocarta.v4.0.symbols.gmt"

anno = gene2set(anno=annoGMT, difSAM_GSE7621_DEG, TermFreq=0)
annotationGenesPop = gene2set(anno=annoGMT, rownames(GSE7621.filtered.expr), TermFreq=0)
annotationGenesPop <- annotationGenesPop[,colnames(anno)]


#GSE7621 
nClust <- 7:20
ncore <- 7
clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")

GSE7621.clena.cluster <- clena(GSE7621.SAM.DEG.expr, 
                               nClust=nClust, 
                               clMethods=clMethods, 
                               metric="correlation", 
                               method="complete",  
                               annotation=anno,  
                               sampleLabel=sampleLabel.GSE7621, 
                               ncore=ncore, 
                               annotationGenesPop=annotationGenesPop, 
                               verbose=TRUE)
save.image("test_clena_chemp.RData")
#save(difSAM_GSE7621_DEG, GSE7621.filtered.expr, GSE7621.SAM.DEG.expr, sampleLabel.GSE7621, file="~/clena/data/PD.RData")

# summary(GSE7621.clena.cluster)
# clusterMethods(GSE7621.clena.cluster)
# nClusters(GSE7621.clena.cluster)
# clusters(GSE7621.clena.cluster, "hierarchical")
# mat(GSE7621.clena.cluster)
# enrichment(GSE7621.clena.cluster, "som", "3")



