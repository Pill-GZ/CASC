# ---------------------------------------------------------------------
# gets a good range for the tuning parameter in CASC
# ---------------------------------------------------------------------

getTuningRangeRohe = function(graphMatrix, covariates, nBlocks,
                          assortative) {

  nCov = ncol(covariates)

  singValCov = svd(covariates, nu = min(nBlocks, nCov))$d

  if(assortative == T) {
    # eigenValGraph = eigs(graphMatrix, nBlocks + 2, which = "LR",
    #     opts = list(retvec = F))$values #eigenvalues only

    eigenValGraph = irlba::irlba(graphMatrix, nu = nBlocks + 1, nv = 0)$d

    if(nCov > nBlocks) {
      hmax = eigenValGraph[1]/(singValCov[nBlocks]^2 - singValCov[nBlocks+1]^2)
    } else {
      hmax = eigenValGraph[1]/singValCov[nCov]^2
    }
    hmin = (eigenValGraph[nBlocks] - eigenValGraph[nBlocks + 1])/singValCov[1]^2
  } else {
    # eigenValGraph = eigs(graphMatrix, nBlocks + 2,
    #     opts = list(retvec = F))$values #eigenvalues only
    # eigenValGraph = sort(eigenValGraph^2, decreasing=T)

    eigenValGraph = (irlba::irlba(graphMatrix, nu = nBlocks + 1, nv = 0)$d)^2

    if(nCov > nBlocks) {
      hmax = eigenValGraph[1]/(singValCov[nBlocks]^2 - singValCov[nBlocks+1]^2)
    } else {
      hmax = eigenValGraph[1]/singValCov[nCov]^2
    }
    hmin = (eigenValGraph[nBlocks] - eigenValGraph[nBlocks + 1])/singValCov[1]^2
  }

  # return
  list( hmax = hmax, hmin = hmin )
}
