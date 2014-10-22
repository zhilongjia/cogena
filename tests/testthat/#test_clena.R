#! /usr/bin/Rscript

#arg options
#c2.cgp.v4.0.symbols.gmt          c2.cp.v4.0.symbols.gmt   mouse.rda
#c2.cp.kegg.v4.0.symbols.gmt      c5.bp.v4.0.symbols.gmt   msigdb.v4.0.symbols.gmt
#c2.cp.reactome.v4.0.symbols.gmt  c6.all.v4.0.symbols.gmt

setwd("/home/zjia/data/Jobbing/PD/src")
load("../result/PD_DEG_analysis.RData")

#require(clena)
require(devtools)
load_all("~/clena/")


# args<-commandArgs(T)
# annoGMT <- args[1]
# if (!exists("annoGMT")) {annoGMT <- "c2.cp.v4.0.symbols.gmt"}

nClust <- 2:4
ncore <- 3
annoGMT <- "c2.cp.v4.0.symbols.gmt"
anno = gene2set(anno=annoGMT, GSE7621.DEG, TermFreq=0)
annotationGenesPop = gene2set(anno=annoGMT, rownames(GSE7621.Explist.filtered$GSE7621$x), TermFreq=0)
annotationGenesPop <- annotationGenesPop[,colnames(anno)]

clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")

#GSE7621 (limma)
sampleLabel.GSE7621 <- apply(GSE7621label[,"title", drop=FALSE],1, c )
GSE7621.clena.cluster <- clena(GSE7621.DEG.expr, nClust=nClust, clMethods=clMethods, metric="correlation", method="complete", validation="enrichment", annotation=anno, maxitems=nrow(GSE7621.DEG.expr), sampleLabel=sampleLabel.GSE7621, ncore=ncore, annotationGenesPop=annotationGenesPop, verbose=TRUE)
