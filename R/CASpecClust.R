#' @useDynLib CASC
#' @importFrom Rcpp sourceCpp
NULL

#' Covariate-assisted Spectral Clustering
#'
#' \code{CASpecClust} returns the covariate-assisted spectral clustering assignments.
#' This is a function that implements the methods described in Binkiewicz et.al. (2016).
#' \code{CASpecClust} takes the adjacency matrix \code{A}, covariate matrix \code{X},
#' calculates a similarity matrix and peforms spectral clustering.
#'
#' The adjacency matrix and covariate matrices are first converted into sparse format
#' to calculate the graph laplacian:
#'  \deqn{L_\tau = D^{−1/2}_\tau A D^{−1/2}_\tau}
#'
#' spectral clusterging is performed on a n-by-n "similarity matrix", which is defined
#' according to method used:
#'
#' \enumerate{
#'   \item regularized spectral clustering:
#'       \deqn{L_\tau}
#'   \item assortative covariate-assisted spectral clustering:
#'       \deqn{L_\tau + \alpha XX^t}
#'   \item covariate-assisted spectral clustering:
#'       \deqn{L_\tau L_\tau + \alpha XX^t}
#' }
#'
#' The nodes are clustered into \code{K} clusters.
#'
#' @param A an adjacency matrix of the nodes, n by n
#' @param X covariate matrix, n by p.
#'        Categorical variables should be re-expressed with dummy variables using
#'        one-hot encoding, see vignettes for example.
#' @param K number of clusters
#' @param n.alpha number of alpha values in the valid range to search. Defaults
#'   to 100
#' @param center whether to center the covariates X. Defaults to FALSE
#' @param method clustering method. See details
#' @param randStarts number of random starts in K-means. Defaults to 20
#'
#' @return A list of three elememts:
#' \enumerate{
#'   \item alpha: values of tuning parameter in the valid range
#'   \item WCSS: within-cluster sum of squares corresponding to the alpha values
#'   \item cluster: the optimal clustering assignments
#' }
#'
#' @examples See vignettes
#'
#' @export

CASpecClust <- function(A, X = NULL, K = 2, n.alpha = 100, center = FALSE,
                        method = 3, randStarts = 20) {

  ## create similarity matrics
  A <- as(A, "dgCMatrix")

  # calculate L.tau (regualrized Laplacian)
  rSums = Matrix::rowSums(A)
  tau = mean(rSums)
  normMat = Matrix::Diagonal(length(rSums), 1 / sqrt(rSums + tau))
  L.tau = normMat %*% A %*% normMat

  if (method == "1") {

    res <- casc_alpha(A, X, L.tau, K, alpha=0, method, randStarts)
    alphaSeq <- NULL
    WCSSalpha <- res$WCSS
    bestAssignment <- res$cluster

  } else {

    # center and scale covariate matrix
    X <- scale(X, center = center, scale = sqrt(Matrix::colSums(X^2)))
    X <- as(X, "dgCMatrix")

    # setting tuning parameters here
    TuningRange <- getTuningRange(L.tau, X, K, method)
    alphaSeq <- seq(TuningRange$alpha.max,
                    TuningRange$alpha.min,
                    length.out = n.alpha)

    # search over the alpha values
    WCSSalpha <- vector("numeric",length(alphaSeq))
    bestAssignment <- vector("numeric",length(alphaSeq))

    for (i in 1:length(alphaSeq)) {

      res <- casc_alpha(A, X, L.tau, K, alphaSeq[i], method, randStarts)
      WCSSalpha[i] <- res$WCSS
      if (res$WCSS == min(WCSSalpha[1:i])) {
        bestAssignment <-  res$cluster
      }

    }

  }

  # return clustering result
  return(list(alpha = alphaSeq, WCSS = WCSSalpha, cluster = bestAssignment+1))
}



