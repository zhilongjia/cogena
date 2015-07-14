# cogena

[![Travis-CI Build Status](https://travis-ci.org/zhilongjia/cogena.png?branch=master)](https://travis-ci.org/zhilongjia/cogena)


To discovery smaller scale, but highly correlated cellular events that may be of great biological relevance, co-expressed gene set enrichment analysis, cogena, clusters gene expression profiles (coExp) and then make enrichment analysis for each clusters (clEnrich) based on hyper-geometric test. The heatmapCluster and heatmapPEI can visualise the results. See vignette for the detailed workflow.


![cogena_workflow](inst/figure/Cogena_workflow.png)

The workflow of cogena

![cogena_heatmapCluster](inst/figure/cogena_heatmapCluster.png)

The heatmap of co-expressed gene set


![cogena_heatmapPEI](inst/figure/cogena_heatmapPEI.png)

 The enrichment score for 10 clusters, together with I, II and All (Down-regulated, Up-regualted and All DE genes)



Installation:

	devtools::install_github("zhilongjia/cogena")

Acknowledgement:

cogena was originally based on the [clValid](http://cran.r-project.org/web/packages/clValid/index.html) package.

