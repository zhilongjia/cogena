#' Gene set enrichment for clusters (for one clustering method and a certain 
#' number of clusters)
#' 
#' Gene set enrichment for clusters sourced from coExp function. the enrichment
#' score are based on -log(p) with p from hyper-geometric test.
#' 
#' Gene sets availiable (See vignette for more): 
#' \itemize{
#' \item c2.cp.kegg.v5.0.symbols.gmt.xz (From Msigdb)
#' \item c2.cp.reactome.v5.0.symbols.gmt.xz (From Msigdb)
#' \item c5.bp.v5.0.symbols.gmt.xz (From Msigdb)
#' \item c2.cp.biocarta.v5.0.symbols.gmt.xz (From Msigdb)
#' \item c2.all.v5.0.symbols.gmt.xz (From Msigdb)
#' \item c2.cp.v5.0.symbols.gmt.xz (From Msigdb)
#' \item c5.mf.v5.0.symbols.gmt.xz (From Msigdb)
#' }
#' 
#' @param genecl_obj a genecl or cogena object
#' @param method as clMethods in genecl function
#' @param nCluster as nClust in cogena function
#' @param annofile gene set annotation file
#' @param sampleLabel sameple Label. Do make the label of interest located after
#' the control label in the order of factor. See details. 
#' @param TermFreq a value from [0,1) to filter low-frequence gene sets
#' 
#' @return a list containing the enrichment score for each clustering methods 
#' and cluster numbers included in the genecl_obj
#' @source
#' Gene sets are from
#' 
#' 1. http://www.broadinstitute.org/gsea/msigdb/index.jsp
#' 
#' 2. http://amp.pharm.mssm.edu/Enrichr/
#' @import parallel
#' @import foreach
#' @import doParallel
#' @examples 
#' 
#' #annotaion
#' annoGMT <- "c2.cp.kegg.v5.0.symbols.gmt.xz"
#' annofile <- system.file("extdata", annoGMT, package="cogena")
#' 
#' data(Psoriasis)
#' clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich_one(genecl_result, "kmeans", "3", annofile=annofile, sampleLabel=sampleLabel)
#' clen_res1 <- clEnrich_one(clen_res, "hierarchical", "2", annofile=annofile, sampleLabel=sampleLabel)
#' 
#' @export
#' 
clEnrich_one <- function(genecl_obj, method, 
                          nCluster, annofile=NULL, sampleLabel=NULL, TermFreq=0){
    
    method <- match.arg(method, c("hierarchical","kmeans","diana","fanny",
                                        "som","model","sota","pam","clara","agnes", "apcluster"), several.ok=FALSE)
    if (any(as.numeric(nCluster)<2)) {
        stop("argument 'nCluster' must be a positive integer vector")
    }
    nCluster <- as.character(nCluster)
    ############################################################################
    # Annotation data
    if (is.null(annofile)) {
        annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt.xz", 
                                package="cogena")
    }
    
    if (is.null(names(sampleLabel))) {
        stop (paste("No name for parameter sampleLabel." ))
    }
    
    annotation <- gene2set(annofile, genecl_obj@labels, TermFreq=TermFreq)
    # the background gene gene-sets matrix
    AllGeneSymbols=NULL
    utils::data(AllGeneSymbols, envir = environment())
    annotationGenesPop <- gene2set(annofile, AllGeneSymbols, TermFreq=TermFreq)
    annotationGenesPop <- annotationGenesPop[,colnames(annotation)]
    
    if (ncol(annotationGenesPop) ==0 || ncol(annotation)==0) {
        stop("Error in annotation as ncol equals zero. 
        Maybe lower the TermFreq value.")
    }
    
    ############################################################################
    # Enrichment score for Up-regulated, Down-regulated genes and All DE genes
    geneExp <- logfc(as.data.frame(genecl_obj@mat), sampleLabel)
    upGene <- rownames(geneExp[geneExp$logFC >0,])
    dnGene <- rownames(geneExp[geneExp$logFC <0,])
    Up <- cogena::PEI(upGene, annotation=annotation, annotationGenesPop=annotationGenesPop)
    Down <- cogena::PEI(dnGene, annotation=annotation, annotationGenesPop=annotationGenesPop)
    All <- PEI(genecl_obj@labels, annotation=annotation, annotationGenesPop=annotationGenesPop)
    upDnPei <- as.matrix(rbind(Up, Down))
    

    ############################################################################
    # Gene sets enrichment analysis for clusters
    if ( class(genecl_obj)  == "genecl" ) {
        clen <- list()
    } else {
        clen <- genecl_obj@measures
    }
    
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
            pei.adjust <- matrix(stats::p.adjust(pei, "fdr"), ncol=ncol(pei))
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
               ncore=genecl_obj@ncore,
               gmt=basename(annofile),
               call=match.call() )
    return (res)
}


