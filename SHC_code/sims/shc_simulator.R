#' @title Data Simulator for SHC Testing
#'
#' This function simulates a dataset comprised of a mixture of t
#' or Gaussian distributions based on the specified parameters.
#' To specify random number of samples to each cluster, pull
#' randomly from 1:K for k*ni times
#' 
#' @param k integer specifying number of clusters (mixtures)
#' @param d integer specifying dimension of problem
#' @param ni integer or integer vector of length \code{k} specifying
#'        number of samples per cluster
#' @param mui numeric matrix of dimension k x d with row ki corresponding
#'        to the mean of cluster ki
#' @param vsd numeric
#' @param df integer specifying the degrees of freedom for each mixture
#'        t-distribution, with NULL corresponding to a Gaussian
#'        distribution (default = NULL)
#'
#' @return
#' a list containing \code{data} and \code{labels}. The data matrix is of
#' dimension (d x n), and the labels are a length n vector of integers, where
#' n is \code{k*ni} or \code{\sum(ni)}. 
#' 
#' @author Patrick Kimes
shc_simulator <- function(k, d, ni, mui, vsd, df = NULL) {

    ## parse cluster size input
    if (length(ni) == 1) {
        ## randomly sample from ni*k
        nVec <- as.numeric(table(sample(1:k, ni*k, replace=TRUE)))
    } else if (length(ni) == 2) {
        nVec <- rmultinom(1, ni[[1]]*k, ni[[2]])[, 1]
    } else {
        stop("ni must be length 1 or 2")
    }
    
    ## create mean matrix and label vector
    muMat <- c()
    labels <- c()
    for (ki in 1:k) {
        muMat <- cbind(muMat, replicate(nVec[ki], mui[ki,]))
        labels <- c(labels, rep(ki, nVec[ki]))
    }

    ## simulate from t-distribution if dof specified else Gaussian
    if (df < 0) {
        data <- muMat + (diag(vsd) %*% matrix(rnorm(d*sum(nVec)), d, sum(nVec)))
    } else {
        data <- muMat + (diag(vsd) %*% matrix(rt(d*sum(nVec), df), d, sum(nVec)))
    }
        
    ## return results as a list
    list(data=data, labels=labels)
}

