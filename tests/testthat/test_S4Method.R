#！ /usr/bin/Rscript

setwd("/home/zjia/data/Jobbing/PD/src")

#require(clena)
require(devtools)
load_all("~/clena")
load("../result/PD_clena_analysis.c2.cp.v4.0.symbols.gmt.RData")
#load("../result/PD_clena_analysis.c2.cp.kegg.v4.0.symbols.gmt.RData")

summary(GSE7621.clena.cluster)
heatmap.PEI(GSE7621.clena.cluster, "hierarchical", "12")
heatmap.PEI(GSE7621.clena.cluster, "som", "5")


summary(GSE7621.SAM.clena.cluster)
heatmap.PEI(GSE7621.SAM.clena.cluster, "sota", "12", orderMethod="max")
heatmap.PEI(GSE7621.SAM.clena.cluster, "diana", "20", CutoffNumGeneset=20)
heatmap.Cluster(GSE7621.SAM.clena.cluster, "sota", "12")

