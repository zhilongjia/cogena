#' stop Implicit Cluster
#'
#' Internal functions.
#'
#' Internal functions for stop cluster. This is a patch for \code{doParallel}.
#' Call \code{stopImplicitCluster()} after foreach/doParallel ending
#' @return NA 
#' @import doParallel
#' @author Dan from fredhutch
#' @keywords internal
#' 
stopImplicitCluster <- function(){
    if(exists(".revoDoParCluster", where=doParallel:::.options) &&
       !is.null(doParallel:::.options[['.revoDoParCluster']]))
    {
        stopCluster(doParallel:::.options[['.revoDoParCluster']])
        remove('.revoDoParCluster', envir=doParallel:::.options)
    }
}
