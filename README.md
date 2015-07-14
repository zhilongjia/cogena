# cogena

[![Travis-CI Build Status](https://travis-ci.org/zhilongjia/cogena.png?branch=master)](https://travis-ci.org/zhilongjia/cogena)


To discovery smaller scale, but highly correlated cellular events that may be of great biological relevance, co-expressed gene set enrichment analysis, cogena, clusters gene expression profiles (coExp) and then makes enrichment analysis for each clusters (clEnrich) based on hyper-geometric test. The heatmapCluster and heatmapPEI can visualise the results. See vignette for the detailed workflow.


![cogena_workflow](inst/figure/Cogena_workflow.png)

The workflow of cogena

![cogena_heatmapCluster](inst/figure/heatmapCluster_Kmeans8.png)

The heatmap of co-expressed gene set based on K-means methods with 8 clusters.


![cogena_heatmapPEI](inst/figure/heatmapPEI_Kmeans8.png)

 The enrichment score for 8 clusters, together with I, II and All (Down-regulated, Up-regualted and All DE genes). The values shown is the -log2(FDR).



Installation:

	devtools::install_github("zhilongjia/cogena")

Acknowledgement:

cogena was originally based on the [clValid](http://cran.r-project.org/web/packages/clValid/index.html) package.

