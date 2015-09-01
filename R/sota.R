#' Self-organizing tree algorithm (SOTA)
#' 
#' Computes a Self-organizing Tree Algorithm (SOTA) clustering of a dataset 
#' returning a SOTA object. This function comes from 
#' \code{\link[clValid]{sota}} in the clValid package with litter change.
#' 
#' The Self-Organizing Tree Algorithm (SOTA) is an unsupervised neural network 
#' with a binary tree topology. It combines the advantages of both hierarchical 
#' clustering and Self-Organizing Maps (SOM). The algorithm picks a node with 
#' the largest Diversity and splits it into two nodes, called Cells. This 
#' process can be stopped at any level, assuring a fixed number of hard 
#' clusters. This behavior is achieved with setting the unrest.growth parameter
#'to TRUE. Growth of the tree can be stopped based on other criteria, 
#'like the allowed maximum Diversity within the cluster and so on. 
#'Further details regarding the inner workings of the algorithm can be 
#'found in the paper listed in the Reference section.
#' 
#' Please note the 'euclidean' is the default distance metric different from 
#' \code{\link[clValid]{sota}}
#' 
#' @param data data matrix or data frame. Cannot have a profile ID as the 
#' first column.
#' @param maxCycles integer value representing the maximum number of iterations 
#' allowed. The resulting number of clusters returned by sota is maxCycles+1 
#' unless unrest.growth is set to FALSE and the maxDiversity criteria is 
#' satisfied prior to reaching the maximum number of iterations
#' @param maxEpochs integer value indicating the maximum number of training 
#' epochs allowed per cycle. By default, maxEpochs is set to 1000.
#' @param distance character string used to represent the metric to be used 
#' for calculating dissimilarities between profiles. 'euclidean' is the default,
#' with 'correlation' being another option.
#' @param wcell alue specifying the winning cell migration weight. 
#' The default is 0.01.
#' @param pcell value specifying the parent cell migration weight. 
#' The default is 0.005.
#' @param scell value specifying the sister cell migration weight. 
#' The default is 0.001.
#' @param delta value specifying the minimum epoch error improvement. 
#' This value is used as a threshold for signaling the start of a new cycle. 
#' It is set to 1e-04 by default.
#' @param neighb.level integer value used to indicate which cells are 
#' candidates to accept new profiles. This number specifies the number of 
#' levels up the tree the algorithm moves in the search of candidate cells 
#' for the redistribution of profiles. The default is 0.
#' @param maxDiversity value representing a maximum variability allowed within 
#' a cluster. 0.9 is the default value.
#' @param unrest.growth logical flag: if TRUE then the algorithm will run 
#' maxCycles iterations regardless of whether the maxDiversity criteria is 
#' satisfied or not and maxCycles+1 clusters will be produced; if FALSE then 
#' the algorithm can potentially stop before reaching the maxCycles based on 
#' the current state of cluster diversities. A smaller than usual number of 
#' clusters will be obtained. The default value is TRUE.
#' @param ... Any other arguments.
#' 
#' @return A SOTA object.
#' \item{data}{data matrix used for clustering}
#' \item{c.tree}{complete tree in a matrix format. Node ID, its Ancestor, and 
#' whether it's a terminal node (cell) are listed in the first three columns. 
#' Node profiles are shown in the remaining columns.}
#' \item{tree}{incomplete tree in a matrix format listing only the terminal 
#' nodes (cells). Node ID, its Ancestor, and 1's for a cell indicator are 
#' listed in the first three columns. Node profiles are shown in the remaining
#' columns.}
#' \item{clust}{integer vector whose length is equal to the number of profiles 
#' in a data matrix indicating the cluster assingments for each profile in the 
#' original order.}
#' \item{totals}{integer vector specifying the cluster sizes.}
#' \item{dist}{character string indicating a distance function used in the 
#' clustering process.}
#' \item{diversity}{vector specifying final cluster diverisities.}
#' @author Vasyl Pihur, Guy Brock, Susmita Datta, Somnath Datta
#' @references  Herrero, J., Valencia, A, and Dopazo, J. (2005). A hierarchical 
#' unsupervised growing neural network for clustering gene expression patterns. 
#' Bioinformatics, 17, 126-136.
#' @examples 
#' #please ref the manual of sota function from clValid.
#' data(Psoriasis)
#' 
#' sotaCl <- sota(as.matrix(DEexprs), 4)
#' @export
#' @keywords cluster
#' @export
#' @rdname sota


sota <- function(data, maxCycles, maxEpochs=1000, distance="euclidean",
    wcell=.01, pcell=.005, scell=.001, delta=.0001, neighb.level=0,
    maxDiversity = .9, unrest.growth=TRUE, ...){

    tree <- sota.init(data)
    pr <- 4:ncol(tree)
    n <- 3
    genes<- 1:nrow(data)
    clust <- rep(1, length(genes))
    Node.Split <- 1

    Res.V <- getResource(data, tree, clust, distance, pr)
    if(distance=="correlation")
        Res.V <- 1-Res.V
    diversity <- Res.V

    for(k in 1:maxCycles){                                #loop for the Cycles
        trainNode <- Node.Split
        trainSamp <- genes[clust==trainNode]
        curr.err <- 1e10
        ep <- 1

        while(ep <= maxEpochs){                      #loop for the Epochs
            last.err <- 0
            left.ctr <- right.ctr <- 0
            left.d <- right.d <- 0
            for(i in trainSamp){
                cells <- tree[c(n-1,n),]
                dist <- rep(0, nrow(cells))
                for(j in 1:2)
                    dist[j] <- dist.fn(data[i,], cells[j,pr], distance=distance)

                or <- which.min(dist)
                if(or==1)
                    left.ctr <- left.ctr + 1
                else
                    right.ctr<- right.ctr + 1

                closest <- cells[or,1]
                sis <- ifelse(closest%%2==0,closest+1,closest-1)
                sis.is.cell <- ifelse(tree[sis,"cell"]== 1, 1, 0)

                ##   updating the cell and its neighbourhood
                if(sis.is.cell==1){
                    parent <- tree[closest, "anc"]
                    tree[closest, pr] <- 
                        tree[closest, pr]+wcell*(data[i,]-tree[closest, pr])
                    tree[sis, pr] <- 
                        tree[sis, pr]+scell*(data[i,]-tree[sis,pr])
                    tree[parent, pr] <- 
                        tree[parent, pr]+pcell*(data[i,]-tree[parent, pr])
                } else {
                    tree[closest, pr] <- 
                        tree[closest, pr]+wcell*(data[i,]-tree[closest, pr])
                }
            }
            cells <- tree[c(n-1,n),]
            for(i in trainSamp){
                for(j in 1:2)
                    dist[j] <- dist.fn(data[i,], cells[j,pr], distance=distance)
                last.err <- last.err+min(dist)}
            last.err <- last.err/length(trainSamp)
            
            if(ifelse(last.err==0, 0, abs((curr.err-last.err)/last.err)) < delta
                && left.ctr !=0 && right.ctr !=0)
                break
            ep <- ep + 1
            curr.err <- last.err
        }
        clust <- assignGenes(data, trainSamp, clust, tree, n, distance, pr, 
            neighb.level)
        Res.V <- getResource(data, tree, clust, distance, pr)
        if(distance=="correlation")
            Res.V <- 1-Res.V

        tempRes <- Res.V
        tempRes[tempRes == 0] <- diversity[tempRes==0]
        diversity <- tempRes
        
        if(k==maxCycles || (max(Res.V) < maxDiversity & unrest.growth==FALSE))
            break  ## do not split the cell
        newCells <- splitNode(Res.V, tree, n)
        tree <- newCells$tree
        n <- newCells$n
        Node.Split <- newCells$toSplit
    }

    tree <- trainLeaves(data, tree, clust, pr, wcell, distance, n, delta)
    Res.V <- getResource(data, tree, clust, distance, pr)
    Res.V <- Res.V[Res.V!=0]
    if(distance=="correlation")
        Res.V <- 1-Res.V
    
    diversity[(length(diversity)-length(Res.V)+1):length(diversity)] <- Res.V

    treel <- tree[tree[,"cell"]==1,]
    old.cl <- treel[,1]
    treel[,1] <- 1:nrow(treel)
    old.clust <- clust

    clust <- cl.ID(old.clust, old.cl, 1:nrow(treel))
    totals <- table(clust)
    #add the gene name to the clust
    names(clust) <- rownames(data)

    out <- list(data=data, c.tree=cbind(tree[1:n,],Diversity=diversity), 
        tree=treel, cluster=clust, totals=totals, dist=distance, 
        diversity=Res.V)

    class(out) <- "sota"
    return(out)
}

sota.init <- function(data){
    nodes <- matrix(0, nrow(data)*2, 3+ncol(data))
    if(is.null(colnames(data)))
        colnames(data) <- paste("V", 1:ncol(data))
    colnames(nodes) <- c("ID", "anc", "cell", colnames(data))
    nodes[,"ID"]=1:(nrow(data)*2)

    nodes[1,] <- c(1, 0, 0, apply(data,2, function(x) mean(x, na.rm=TRUE)))
    nodes[2,] <- c(2, 1, 1, nodes[1,][-c(1,2,3)])
    nodes[3,] <- c(3, 1, 1, nodes[1,][-c(1,2,3)])
    return(nodes)
}

dist.fn <- function(input, profile, distance){
    if(distance=="correlation")
        return(1-cor(input,profile, use="pairwise.complete.obs"))
    else
        return(sqrt(sum((input-profile)^2)))
}

cl.ID <- function(clust, old.cl, new.cl){
    for(i in 1:length(clust))
        clust[i] <- new.cl[which(old.cl==clust[i])]
    clust
}

getResource <- function(data, tree, clust, distance, pr){
    dist <- rep(0, length(clust))
    resource <- rep(0, max(clust))

    for(i in unique(clust)){
        temp <- data[clust==i,]
        if(is.vector(temp))
            temp <- matrix(temp, nrow=1, ncol=ncol(data))
        if(distance=="correlation")
            resource[i] <- mean(apply(temp, 1, dist.fn, profile=tree[i,pr],
                distance=distance))
        else
            resource[i] <- mean(apply(temp, 1, dist.fn, profile=tree[i,pr],
                distance=distance))}
    resource
}

getCells <- function(tree, neighb.level, n){
    or.n <- n
    cells <- c(n-1,n)
    for(i in 1:(neighb.level+1)){
        n  <- tree[n, "anc"]
        if(n==1)
            break
    }
    for(j in 2:(or.n-2)){
        z <- j
        if(tree[j,"cell"]!=1)
            next
        while(z > 0){
            z <- tree[z, "anc"]
            if(z==n){
                cells <- c(cells, j)
                break}
        }
    }
    return(tree[cells,])
}

trainLeaves <- function(data, tree, clust, pr, wcell, distance, n, delta){
    nc <- ncol(data)
    for(i in 1:n){
        if(!is.element(i, clust))
            next
        temp <- matrix(data[clust==i,], ncol=nc)
        converged <- FALSE
        init.err <- getCellResource(temp, tree[i,pr], distance)
        while(!converged){
            for(j in 1:nrow(temp))
                tree[i, pr] <- tree[i, pr]+wcell*(temp[j,]-tree[i, pr])
            
            last.err <- getCellResource(temp, tree[i,pr], distance)
            converged <- 
                ifelse(abs((last.err-init.err)/last.err) < delta, TRUE, FALSE)
            init.err <- last.err
        }
    }
    return(tree)
}


assignGenes <- function(data, Sample, clust, tree, n, distance, pr, 
    neighb.level){
    if(neighb.level==0)
        cells <- tree[c(n-1,n),]
    else
        cells <- getCells(tree, neighb.level, n)

    for(i in Sample){
        dist <- rep(0, nrow(cells))
        for(j in 1:nrow(cells))
            dist[j] <- dist.fn(data[i,], cells[j,pr], distance)
        or <- which.min(dist)
        closest <- cells[or,1]
        clust[i] <- closest
    }
    clust
}

splitNode <- function(Res.V, tree, n){
    maxheter <- which.max(Res.V)
    cl.to.split <- tree[maxheter,1]
    tree[n<-n+1,-1] <- tree[cl.to.split,-1]
    tree[n, "anc"] <- cl.to.split
    tree[n<-n+1,-1] <- tree[cl.to.split,-1]
    tree[n, "anc"] <- cl.to.split
    tree[cl.to.split, "cell"] <- 0
    return(list(tree=tree, n=n, toSplit=cl.to.split))
}



getCellResource <- function(temp, profile, distance){
    if(distance=="correlation")
        resource <- mean(apply(temp, 1, dist.fn, profile, distance=distance))
    else(distance=="euclidean")
    resource <- mean(apply(temp, 1, dist.fn, profile, distance=distance))
    resource
}

#' @param x an object of sota
#' @method print sota
#' @rdname sota
#' @export
print.sota <- function(x, ...){
    results <- as.matrix(cbind(as.numeric(names(x$totals)), 
        as.numeric(x$totals), x$diversity))  ## changed
    colnames(results) <- c("ID","Size", "Diversity")
    rownames(results) <- rep("", nrow(results))
    cat("\nClusters:\n")
    print(results, ...)
    cat("\nCentroids:\n")

    print(format(data.frame(x$tree[,-c(1:3)])), ...)
    cat("\n")
    cat(c("Distance: ", x$dist, "\n"))
    invisible(x)
}

#' @param cl cl specifies which cluster is to be plotted by setting it to the 
#' cluster ID. By default, cl is equal to 0 and the function plots all clusters 
#' side by side.
#' @method plot sota
#' @rdname sota
#' @export
plot.sota <- function(x, cl=0, ...){

    op <- par(no.readonly=TRUE)
    on.exit(par(op))
    if(cl!=0)
        par(mfrow=c(1,1)) else
        {
            pdim <- c(0,0)
            for(i in 1:100){
                j <- i
                if(length(x$totals) > i*j)
                    j <- j+1
                else{
                    pdim <- c(i,j)
                    break}
                if(length(x$totals) > i*j)
                    i <- i+1
                else{
                    pdim <- c(i,j)
                    break}
            }
            par(mfrow=pdim)
        }

    ylim = c(min(x$data), max(x$data))
    pr <- 4:ncol(x$tree)
    if(cl==0)
        cl.to.print <- 1:length(table(x$clust)) else  ## changed
            cl.to.print <- cl
    cl.id <- sort(unique(x$clust))  ## changed

    for(i in cl.to.print){
        plot(1:ncol(x$data), x$tree[i, pr], col="red", type="l",
            ylim=ylim, xlab=paste("Cluster ",i), ylab="Expr. Level", ...)
        legend("topleft", legend=paste(x$totals[i], " Genes"), cex=.7,
            text.col="navy", bty="n")
        cl <- x$data[x$clust==cl.id[i],]  ## changed
        if(is.vector(cl))
            cl <- matrix(cl, nrow=1)
        for(j in 1:x$totals[i])
            lines(1:ncol(x$data), cl[j,], col="grey")
        lines(1:ncol(x$data), x$tree[i, pr], col="red", ...)
        
    }
}
