#' Covariate-assisted spectural clustering using fixed weight alpha
#'
#' \code{casc_alpha} returns the covariate-assisted pectrual clustering
#' assignments with fixed alpha weight.
#'
#'
#' @param A an adjacent matrix of the nodes, n-by-n
#' @param X covariate matrix, n-by-p
#' @param alpha the weight of the XX' matrix to be added to regularized
#'   Laplacian. Default value 100
#' @param method xlustering method. See details
#'
#' @return The clustering assignments along with the within-cluster sum of
#'   squares.
#'
#' @examples
#'
#' \dontrun{
#' sum("a")
#' }

casc_alpha <- function(A, X, L.tau, K, alpha, method, randStarts) {

  # input matrix L.tau or L.alpha(assortive or not)
  if (method == "1"){

    input.matrix <- L.tau

  } else if (method == "2") {

    input.matrix <- L.tau + alpha * X %*% Matrix::t(X)

  } else if (method == "3") {

    input.matrix <- L.tau %*% L.tau + alpha * X %*% Matrix::t(X)

  } else {

    stop(paste("Error: method =", method, "Not valid"))
  }

  N <- ncol(A)

  # form matrix with eigenvectors of similarity matrix as columns
  Ustar <- RSpectra::eigs_sym(input.matrix, k=K)$vectors

  # normalize the rows
  Ustar <- t(scale(t(Ustar), center=FALSE, scale = sqrt(Matrix::rowSums(Ustar^2))))

  #run k-means here and assign.
  rcpp_kmeans(Ustar, randStarts)

}
