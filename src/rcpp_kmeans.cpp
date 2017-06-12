#include <Rcpp.h>
#include <time.h>
using namespace Rcpp;

// [[Rcpp::export]]

List rcpp_kmeans(NumericMatrix x, int nstart) {
  //random seed
  srand (time(NULL));

  // x is a N by K matrix of orthogonal columns
  int N = x.nrow();
  int K = x.ncol();

  // create cluster assignment matrix for nstart
  NumericMatrix clust_assign(N, nstart);
  // sum of within cluster distances for each random initialization
  double WCSS[nstart];

  // cluster assignment to be returned
  NumericVector cluster_assignment(N);

  // repeat for nstart different initialization in kmeans
  for (int iter = 0; iter < nstart; iter++) {
    // create centers for clusters
    // each row is a k-dim vector of the coordinates
    NumericMatrix centers(K, K);
    // create storage for cluster sizes
    int clust_size[K];
    // and set it to zeros
    for (int k = 0; k < K; k++) {
      clust_size[k] = 0.0;
    }
    // create storage for cluster assignments
    int cluster[N];

    // initialize the centers by picking at random
    // K out of the N points
    std::set<int> init;
    while (init.size() != K) {
      init.insert(rand() % N);
    }
    std::vector<int> init2(init.begin(), init.end());

    // initialize the centroids coordinates
    for (int k = 0; k < K; k++) {
      for (int j = 0; j < K; j++) {
        centers(k,j) = x(init2[k],j);
      }
    }

    // initialize cluster assignments of the N points
    for (int i = 0; i < N; i++) { // for each point
      double dist[K];
      for (int j = 0; j < K; j++) { // for each centroid
        dist[j] = 0.0;
      }
      // calculate distance to centroids
      for (int j = 0; j < K; j++) {
        for (int k = 0; k < K; k++) {
          dist[j] += (x(i,k) - centers(j,k))
          * (x(i,k) - centers(j,k));
        }
      }
      // find nearest centroid
      int assign = 0;
      double m = dist[0];
      for (int j = 1; j < K; j++) {
        if (dist[j] < m) {
          assign = j;
          m = dist[j];
        }
      }
      // and assign the index of the nearest centroid
      cluster[i] = assign;
    } // end of initialize cluster assignments of the N points

    // create storage for updated centroids
    NumericMatrix update_centers(K, K);
    // update the coordinates of centroids
    for (int i = 0; i < N; i++) {
      clust_size[cluster[i]] += 1;
      for (int j = 0; j < K; j++) {
        update_centers(cluster[i],j) = ((clust_size[cluster[i]] - 1)
                                          / (double) clust_size[cluster[i]])
        * update_centers(cluster[i],j)
        + (x(i,j) / (double) clust_size[cluster[i]]);
      }
    }
    centers = update_centers; // end of centroid updates

    // `change' is a flag for assignment change
    int change = 1;
    // sum of within cluster distance
    double sum_dist;

    // repeat untill converge
    while (change > 0) {
      sum_dist = 0.0;
      change = 0;


      for (int i = 0; i < N; i++) { // for each point
        double dist[K];
        for (int j = 0; j < K; j++) { // for each centroid
          dist[j] = 0.0;
        }
        // calculate distance to centroids
        for (int j = 0; j < K; j++) {
          for (int k = 0; k < K; k++) {
            dist[j] += (x(i,k) - centers(j,k))
            * (x(i,k) - centers(j,k));
          }
        }
        // find nearest centroid
        int assign = 0;
        double m = dist[0];
        for (int j = 1; j < K; j++) {
          if (dist[j] < m) {
            assign = j;
            m = dist[j];
          }
        }
        // increment the sum of distances by
        // the distance to the nearest centroid
        sum_dist += m;
        // if cluster assignments change, update the centroid
        if (assign != cluster[i]) {
          change += 1;
          for (int j = 0; j < K; j++) { // for each coordinate
            // update the cluster that incremented
            update_centers(assign,j) = (clust_size[assign]
                                          / (double) (clust_size[assign] + 1))
            * update_centers(assign,j)
            + (x(i,j)
            / (double) (clust_size[assign] + 1));
            // update the cluster that decremented
            update_centers(cluster[i],j) =
            (clust_size[cluster[i]]
               / (double) (clust_size[cluster[i]] - 1))
              * update_centers(cluster[i],j)
              - (x(i,j)
              / (double) (clust_size[cluster[i]]
              - 1));
          } // end of centroid updates
          clust_size[assign] += 1;
          clust_size[cluster[i]] -= 1;
          cluster[i] = assign;
        } // end of if
      } // end of loop for each point

      centers = update_centers;

    } // end of while loop (repeat untill convergence)

    // store the within cluster sum of squares
    // for this iteration
    WCSS[iter] = sum_dist;
    for (int i = 0; i < N; i++) {
      clust_assign(i,iter) = cluster[i];
    }

  } // end of the nstart iterations

  // finding the assignemt that minimizes the WCSS
  int index = 0;
  double min_WCSS = WCSS[0];
  for (int iter = 1; iter < nstart; iter++) {
    if (WCSS[iter] < min_WCSS) {
      index = iter;
      min_WCSS = WCSS[iter];
    }
  }

  // return the best assigment
  cluster_assignment = clust_assign(_,index);
  List ret;
  ret["WCSS"] = min_WCSS; ret["cluster"] = cluster_assignment;
  return(ret);
}

/*
 * given Matrix615<double> x(N,K), nBlocks=K, nstart
 * Let centers K by K matrix : each row represents the centroid
 * cluster: clustering assignment N by 1 vector
 * Input Matrix615<double> x(N,K)=Ustar, nstart
 */



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
x <- matrix(c(1,2,4,5,1,1,3,4),4,2)
rcpp_kmeans(x,5)

*/
