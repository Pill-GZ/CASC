#' Get range of tuning weight alpha
#'
#' \code{getTuningRange} determins the valid range of tuning variable alpha to
#' search over
#'
#' @param A an adjacent matrix of the nodes, n-by-n
#' @param X covariate matrix, n-by-p
#' @param K number of clusters
#' @param method
#'
#' @return the max and min of possible alpha values
#'
#' @examples
#'
#' \dontrun{
#' sum("a")
#' }


getTuningRange <- function(L.tau, X, K, method) {

  n.cov <- ncol(X)

  # RSpectra::svds works with matrices with dim > 2 only
  if (n.cov > 2 && K < n.cov) {
    singValCov <- RSpectra::svds(X, k = min(K + 1, n.cov))$d
  } else {
    singValCov <- svd(as.matrix(X), nu = min(K + 1, n.cov))$d
  }

  # calculate permissible range of alpha
  if (method == "2") {
    eigenValGraph <- RSpectra::eigs_sym(L.tau, k = K + 1)$values
    if (n.cov > K) {
      alpha.max <- eigenValGraph[1] /
        (singValCov[K] ^ 2 - singValCov[K + 1] ^ 2)
    } else {
      alpha.max <- eigenValGraph[1] / singValCov[n.cov] ^ 2
    }
    alpha.min <- (eigenValGraph[K] - eigenValGraph[K + 1]) /
      singValCov[1] ^ 2
  } else if (method == "3") {
    eigenValGraph <- (RSpectra::eigs_sym(L.tau, k = K + 1)$values) ^ 2
    if (n.cov > K) {
      alpha.max <- eigenValGraph[1] /
        (singValCov[K] ^ 2 - singValCov[K + 1] ^ 2)
    } else {
      alpha.max <- eigenValGraph[1] / singValCov[n.cov] ^ 2
    }
    alpha.min <- (eigenValGraph[K] - eigenValGraph[K + 1]) /
      singValCov[1] ^ 2
  }

  # return
  list(alpha.max = alpha.max, alpha.min = alpha.min)
}
