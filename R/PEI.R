#' Pathway enrichment index
#' 
#' Internal clena functions.

PEI <- function(genenames, annotation, annotationGenesPop, verbose=FALSE) {

    ######################################################################################
    # presteps for calculating p-value
    #pei <- vector("numeric", ncol(annotation))
    pei <- vector(mode="numeric")
    NumGenes <- nrow(annotationGenesPop)
    for (j in colnames(annotation) ){
        NumGenesWithinPathway <- sum(annotationGenesPop[,j])
        NumGenesWithinCluster <- length(genenames)
        NumGenesWithinCluster.Pathway <- sum(annotation[genenames,j]) -1 #minus one due to phyper
        if (NumGenesWithinPathway<2 | NumGenesWithinCluster<2 | NumGenesWithinCluster.Pathway<2){
            pei[j] = NA
        } else {
            pei[j] <- phyper(NumGenesWithinCluster.Pathway, NumGenesWithinPathway, NumGenes-NumGenesWithinPathway, NumGenesWithinCluster, lower.tail=FALSE)
        }
    }
    pei
}
