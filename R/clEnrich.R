#' Gene set enrichment for clusters
#' 
#' Gene set enrichment for clusters sourced from coExp function. the enrichment
#' score are based on -log2(p) with p from hyper-geometric test.
#' 
#' sampleLable:
#' Use factor(c("Normal", "Cancer", "Normal"), levels=c("Normal", "Cancer")), 
#' instead of factor(c("Normal", "Cancer","Normal")). This parameter will affect
#' the direction of gene regulation in cogena.
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
#' @param genecl_obj a genecl object
#' @param annofile gene set annotation file
#' @param sampleLabel sameple Label. Do make the label of interest located after
#' the control label in the order of factor. See details. 
#' @param TermFreq a value from [0,1) to filter low-frequence gene sets
#' @param ncore the number of cores used
#' 
#' @return a list containing the enrichment score for each clustering methods 
#' and cluster numbers included in the genecl_obj
#' 
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
#' utils::data(Psoriasis)
#' clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#'     
#' @export
#' 
clEnrich <- function(genecl_obj, annofile=NULL, sampleLabel=NULL, TermFreq=0, ncore=1){
    
    ############################################################################
    # Annotation data
    if (is.null(annofile)) {
        annofile <- system.file("extdata", "c2.cp.kegg.v6.1.symbols.gmt.xz", 
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
    
    # check intersect input genes and AllGeneSymbols
    input_gene_len <- nrow(genecl_obj@mat)
    intersect_gene_len <- length(intersect(rownames(genecl_obj@mat), AllGeneSymbols))
    Biobase::note(paste(intersect_gene_len, "out of", input_gene_len, "exist in the genes population.") )
    
    
    if (ncol(annotationGenesPop) ==0 || ncol(annotation)==0) {
        stop("Error in annotation as ncol equals zero. 
        Maybe lower the TermFreq value.")
    }
    
    
    ############################################################################
    ############################################################################
    clen <- list()
    
    ############################################################################
    # Enrichment score for Up-regulated, Down-regulated genes and All DE genes
    geneExp <- logfc(as.data.frame(genecl_obj@mat), sampleLabel)
    upGene <- rownames(geneExp[geneExp$logFC >0,])
    dnGene <- rownames(geneExp[geneExp$logFC <0,])
    Up <- cogena::PEI(upGene, annotation=annotation, annotationGenesPop=annotationGenesPop)
    Down <- cogena::PEI(dnGene, annotation=annotation, annotationGenesPop=annotationGenesPop)
    All <- PEI(genecl_obj@labels, annotation=annotation, annotationGenesPop=annotationGenesPop)
    upDnPei <- as.matrix(rbind(Up, Down))
    
    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)
    nc <- NULL
    ############################################################################
    # Gene sets enrichment analysis for clusters
    for (i in clusterMethods(genecl_obj) ) {
        pei_tmp <- foreach::foreach (nc = as.character(nClusters(genecl_obj)) ) %dopar% {
            cluster <- cogena::geneclusters(genecl_obj, i, nc)
            if (nc != length(unique(cluster))) {
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
            return (pei)
        }
        names(pei_tmp) <- nClusters(genecl_obj)
        clen[[i]] <- pei_tmp
    }
    
    # stopImplicitCluster()
    parallel::stopCluster(cl)
    
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


