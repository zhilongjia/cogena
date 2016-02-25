#' Enrichment score or Similarity based on PAGE for clusters
#' 
#' Enrichment score (-log2(FDR)) or Similarity based on Parametric Analysis of 
#' Gene Set Enrichment (PAGE) for clusters.
#' 
#' The CMap rank data are used internally (top ranked genes are Up-regulated).
#' 
#' @inheritParams clEnrich_one
#' @param ncore the number of cores used
#' @param verbose verbose
#' @import parallel
#' @import foreach
#' @import doParallel
#' @export
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
#' cogena_res <- clPage(genecl_result, "kmeans", "3", sampleLabel=sampleLabel)
#' 
#' 
clPage <- function(genecl_obj, method, nCluster, sampleLabel=NULL, ncore=1, similarity=FALSE, verbose=FALSE){
    
    method <- match.arg(method, c("hierarchical","kmeans","diana","fanny",
                                  "som","model","sota","pam","clara","agnes", "apcluster"), several.ok=FALSE)
    if (any(as.numeric(nCluster)<2)) {
        stop("argument 'nCluster' must be a positive integer vector")
    }
    nCluster <- as.character(nCluster)
    
    data(CMap_rank)
    dat <- scale(CMap_rank, scale=FALSE)
    
    # genes in cluster
    cl_gene_list <- list()
    geneExp <- logfc(as.data.frame(genecl_obj@mat), sampleLabel)
    upGene <- rownames(geneExp[geneExp$logFC >0,])
    dnGene <- rownames(geneExp[geneExp$logFC <0,])
    for (i in seq_len(nCluster)) {
        cl_i_gene <- cogena::geneInCluster(genecl_obj, method, nCluster, as.character(i))
        cl_gene_list[[paste0(i, "#", length(cl_i_gene))]] <- intersect(cl_i_gene, rownames(dat) )
    }
    cl_gene_list[[paste0("Up#", length(upGene))]] <- intersect(upGene, rownames(dat) )
    cl_gene_list[[paste0("Down#", length(dnGene))]] <- intersect(dnGene, rownames(dat) )
    
    #PAGE
    sigma <- apply(dat, 2, sd, na.rm=TRUE)
    mu <- apply(dat, 2, mean, na.rm=TRUE)
    page_single <- function(dat, gene_name, geneExp) {
        gene_name <- intersect(gene_name, rownames(dat))
        sm <- apply(dat[gene_name,], 2, mean, na.rm=TRUE)
        zscore <- (sm-mu)*length(gene_name)^(1/2) / sigma
        # anti-correlation is negative
        if (sum(geneExp[gene_name,"logFC"])>0) {
            zscore= -zscore
        }
        zscore
    }
    
    # Zscore for each cluster
    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)
    if (verbose) {print(paste("getDoParWorkers:", foreach::getDoParWorkers()))}
    nc=NULL
    Zscore <- foreach(nc=1:length(cl_gene_list), .combine=cbind) %dopar% {
        cl_i_gene <- cl_gene_list[[names(cl_gene_list)[nc]]]
        page_single(dat, cl_i_gene, geneExp)
    }
    parallel::stopCluster(cl)
    
    Zscore <- cbind(Zscore, rowMeans(Zscore[,tail(colnames(Zscore), 2)]))
    colnames(Zscore) <- c(names(cl_gene_list), paste0("All#", nrow(genecl_obj@mat)) )
    
    # similarity
    page_max <- function(top){
        i <- 1 # for rank input, i is irrelative
        gene_name <- names(sort(dat[,i], decreasing=TRUE)[1:top])
        (mean(dat[gene_name,i])-mean(dat[,i]))*top^(1/2) / sd(dat[,i])
    }
    factor_value <- sapply(sapply(cl_gene_list, length), page_max)
    factor_value_all <- page_max(nrow(genecl_obj@mat))
    factor_value <- c(factor_value, factor_value_all)
    # Need debug
    sim <- Zscore %*% diag(1/factor_value)
    colnames(sim) <- colnames(Zscore)
    
    # enrichment score: -log2 (fdr)
    pval <- pnorm(Zscore)
    fdr <- p.adjust(pval, method = "fdr")
    fdr <- matrix(fdr, nrow=nrow(pval), ncol=ncol(pval), dimnames = dimnames(pval))
    score <- t(-log2(fdr) )
    
    rownames(score) <- sub("#.*", "", rownames(score) )
    
    if ( class(genecl_obj)  == "genecl" ) {
        clen <- list()
    } else {
        clen <- genecl_obj@measures
    }
    clen[[method]][[nCluster]] <- score
    
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
               annotation=sim,
               sampleLabel=as.factor(sampleLabel),
               nClust=genecl_obj@nClust, 
               metric=genecl_obj@metric, 
               method=genecl_obj@method, 
               ncore=genecl_obj@ncore,
               gmt="CMap_rank",
               call=match.call() )
    return (res)
}
