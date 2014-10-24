#！ /usr/bin/Rscript

setwd("/home/zjia/data/Jobbing/PD/src")

#require(clena)
require(devtools)
load_all("~/clena")

load("~/clena/tests/testthat/test_clena.RData")
#load("~/clena/tests/testthat/test_clena_chemp.RData")


#load("../result/PD_clena_analysis.c2.cp.v4.0.symbols.gmt.RData")
load("/home/zjia/data/Jobbing/PD/result/PD_clena_analysis.c2.cp.kegg.v4.0.symbols.gmt.RData")

summary(GSE7621.clena.cluster)
heatmapPEI(GSE7621.clena.cluster, "hierarchical", "20")
heatmapCluster(GSE7621.clena.cluster, "kmeans", "3")
geneInCluster(GSE7621.clena.cluster, "sota", "3", "2")

heatmapPEI(GSE7621.SAM.clena.cluster, "som", "3")
tmp <- enrichment(GSE7621.SAM.clena.cluster, "som", "3")
