#' Gene set enrichment for clusters (for one clustering method and a certain 
#'	number of clusters)
#' 
#' Gene set enrichment for clusters sourced from coExp function. the enrichment
#' score are based on -log(p) with p from hyper-geometric test.
#' 
#' Gene sets availiable (See vignette for more): 
#' \itemize{
#' \item c2.cp.kegg.v5.0.symbols.gmt (From Msigdb)
#' \item c2.cp.kegg.v4.0.symbols.gmt (From Msigdb)
#' \item c2.cp.reactome.v5.0.symbols.gmt (From Msigdb)
#' \item c5.bp.v5.0.symbols.gmt (From Msigdb)
#' \item c2.cp.biocarta.v5.0.symbols.gmt (From Msigdb)
#' \item c2.all.v5.0.symbols.gmt (From Msigdb)
#' \item c2.cp.v5.0.symbols.gmt (From Msigdb)
#' \item c5.mf.v5.0.symbols.gmt (From Msigdb)
#' \item Transcription_Factor_PPIs.gmt (From Enrichr)
#' \item TargetScan_microRNA.gmt (From Enrichr)
#' }
#' 
#' @param genecl_obj a genecl object
#' @param method as clMethods in genecl function
#' @param as nClust in cogena function
#' @param annofile gene set annotation file
#' @param sampleLabel sameple Label
#' @param TermFreq a value from [0,1) to filter low-frequence gene sets
#' @param ncore the number of cores used
#' 
#' @return a list containing the enrichment score for each clustering methods 
#' and cluster numbers included in the genecl_obj
#' @import parallel
#' @import foreach
#' @import doParallel
#' @examples 
#' 
#' #annotaion
#' annoGMT <- "c2.cp.kegg.v5.0.symbols.gmt"
#' annofile <- system.file("extdata", annoGMT, package="cogena")
#' 
#' data(PD)
#' clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich_one(genecl_result, "kmeans", "3", annofile=annofile, sampleLabel=sampleLabel)
#' 
#' @export
#' 
clEnrich_one <- function(genecl_obj, method=clusterMethods(genecl_obj), 
                          nCluster=nClusters(genecl_obj), annofile=NULL, sampleLabel=NULL, TermFreq=0, ncore=1){
    
    ############################################################################
    # Annotation data
    if (is.null(annofile)) {
        annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", 
                                package="cogena")
    }
    annotation <- gene2set(annofile, genecl_obj@labels, TermFreq=TermFreq)
    # the background gene gene-sets matrix
    AllGeneSymbols=NULL
    data(AllGeneSymbols, envir = environment())
    annotationGenesPop <- gene2set(annofile, AllGeneSymbols, TermFreq=TermFreq)
    annotationGenesPop <- annotationGenesPop[,colnames(annotation)]
    
    if (ncol(annotationGenesPop) ==0 || ncol(annotation)==0) {
        stop("Error in annotation as ncol equals zero. 
        Maybe lower the TermFreq value.")
    }
    
    ############################################################################
    # Enrichment score for Up-regulated, Down-regulated genes and All DE genes
    geneExp <- as.data.frame(genecl_obj@mat)
    geneExp$meanA <- apply(geneExp[,names(sampleLabel)[which((sampleLabel == names(table(sampleLabel))[1]))]], 1, mean)
    geneExp$meanB <- apply(geneExp[,names(sampleLabel)[which((sampleLabel == names(table(sampleLabel))[2]))]], 1, mean)
    geneExp$logFC <- ifelse( log2(geneExp$meanB/geneExp$meanA)>0, 1, -1)
    
    upGene <- rownames(geneExp[geneExp$logFC >0,])
    dnGene <- rownames(geneExp[geneExp$logFC <0,])
    Up <- cogena::PEI(upGene, annotation=annotation, annotationGenesPop=annotationGenesPop)
    Down <- cogena::PEI(dnGene, annotation=annotation, annotationGenesPop=annotationGenesPop)
    All <- PEI(genecl_obj@labels, annotation=annotation, annotationGenesPop=annotationGenesPop)
    upDnPei <- as.matrix(rbind(Up, Down))
    

    ############################################################################
    # Gene sets enrichment analysis for clusters
    clen <- list()
    cluster <- cogena::geneclusters(genecl_obj, method, nCluster)
    if (nCluster != length(unique(cluster))) {
        # warning (paste("Cluster", nc, "(aim) only have", length(unique(cluster)), "(result) clusters"))
        pei <- NA
    } else {
        pei <- matrix(NA, nrow=length(unique(cluster)), 
                      ncol=ncol(annotation))
        rownames(pei) <- sort(unique(cluster))
        colnames(pei) <- colnames(annotation)
        for (k in  sort(unique(cluster))) {
            genenames <- names(which(cluster==k))
            pei[as.character(k),] <- cogena::PEI(genenames, annotation=annotation, 
                                                 annotationGenesPop=annotationGenesPop)
        }
        pei <- rbind(pei, upDnPei, All)
        
        # negative log2 p value
        logAdjPEI <- function (pei) {
            # fdr based on pval
            pei.adjust <- matrix(p.adjust(pei, "fdr"), ncol=ncol(pei))
            dimnames(pei.adjust) <- dimnames(pei)
            pei.NeglogPval <- -log2(pei.adjust)
        }
        pei <- logAdjPEI(pei)
    }
    clen[[method]][[nCluster]] <- pei
    
    upDn <- list(upGene=upGene, dnGene=dnGene)
    
    ############################################################################
    ############################################################################
    res <- new("cogena", 
               mat=genecl_obj@mat, 
               measures=clen,
               upDn=upDn,
               Distmat=genecl_obj@Distmat, 
               clusterObjs=genecl_obj@clusterObjs, 
               clMethods=genecl_obj@clMethods, 
               labels=genecl_obj@labels, 
               annotation=annotation,
               sampleLabel=as.factor(sampleLabel),
               nClust=genecl_obj@nClust, 
               metric=genecl_obj@metric, 
               method=genecl_obj@method, 
               ncore=ncore,
               gmt=basename(annofile),
               call=match.call() )
    return (res)
}


