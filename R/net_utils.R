
#' Returns spectral cluster memberships
#' @param x an adjacency matrix
#' @param n number of rows/columns in x
#' @param k number of clusters
#' @keywords internal
#' @noRd
dcspectral <- function(x, n, k){
  d <- rowSums(x)
  d <- d + mean(d)
  deg <- diag(1/sqrt(d))
  l <- deg %*% x %*% deg
  if(!isSymmetric(l)){
    l <- (l + t(l))/2
  }
  spectra <- eigen(l)
  specmat <- spectra$vectors[, 1:k]
  rownorm <- apply(specmat, 1, function(a){(sum(a^2))^0.5})
  rownorm <- ifelse(rownorm < 10^(-06), 10^(-06), rownorm)
  specnorm <- specmat/rownorm
  speck <- kmeans(specnorm, k, nstart=5)
  return(speck$cluster)
}

#' Returns probability matrix for SBM
#' @param n number of vertices
#' @param k number of clusters
#' @param p community probability matrix
#' @param clusters cluster membership vector
#' @keywords internal
#' @noRd
get_P_SBM <- function(n, k, p, clusters){
  Z <- matrix(0, n, k)
  for(i in 1:n){
    Z[i, clusters[i]] <- 1
  }
  P <- Z%*%p%*%t(Z)
  return(P)
}

#' Returns network sampled from a SBM
#' @param n number of vertices
#' @param P probability matrix
#' @keywords internal
#' @noRd
sample_from_SBM <- function(n, P){
  A <- matrix(0, n, n)
  for(l in 1:(n - 1)){
    for(j in (l + 1):n){
      A[l, j] <- rbinom(1,1, P[l, j])
      A[j, l] <- A[l, j]
    }
  }
  G <- igraph::graph_from_adjacency_matrix(A, mode= "undirected")
  return(G)
}

#' Returns probability matrix for DCSBM
#' @param n number vertices
#' @param k number of clusters
#' @param theta degree parameter vector
#' @param p community probability matrix
#' @param clusters cluster membership vector
#' @keywords internal
#' @noRd
get_P_DCSBM <- function(n, k, theta, p, clusters){
  Z <- matrix(0, n, k)
  for(i in 1:n){
    Z[i, clusters[i]] <- 1
  }
  P <- diag(theta)%*%Z%*%p%*%t(Z)%*%diag(theta)
  return(P)
}

#' Returns network sampled from a DCSBM
#' @param n number of vertices
#' @param P probability matrix
#' @keywords internal
#' @noRd
sample_from_DCSBM <- function(n, P){

  A <- matrix(0, n, n)

  for(l in 1:(n - 1)){
    for(j in (l + 1):n){
      A[l, j] <- ifelse(rpois(1, P[l, j]) >= 1, 1, 0)
      A[j, l] <- A[l, j]
    }
  }
  G <- igraph::graph_from_adjacency_matrix(A)
  return(G)
}

#' Returns mean distance in Chung-Lu network
#' @param g igraph network
#' @keywords internal
#' @noRd
cl_dist <- function(g){
  if(igraph::is_connected(g)){
    return(igraph::mean_distance(g))
  }
  dist.mat <- igraph::distances(g)
  #Get off diagonal elements
  dist <- dist.mat[row(dist.mat) != col(dist.mat)]
  #Replace inf distances with diamater
  dist[dist == Inf] <- igraph::diameter(g)
  return(mean(dist))
}

#' Returns network sampled from the CL model
#' @param P probability matrix
#' @keywords internal
#' @noRd
sample_from_CL <- function(P){
  n <- nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  A[col(A) > row(A)] <- runif(n*(n - 1)/2)
  A <- A + t(A)
  diag(A) <- runif(n)
  A <- (A < P) + 0
  G <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  return(G)
}

#' Prints output of sw_test
#' @param res result object
#' @param baseline.model model against which comparison was conducted
#' @keywords internal
#' @noRd
sw_print_res <- function(res, baseline.model){
  cat("\nTest against ", baseline.model, "\n")
  cat("Clustering coefficient: ", res[[3]], ", p-value for clustering coefficient: ", res[[1]], "\n")
  cat("Average path length: ", res[[4]], ", p-value for clustering coefficient: ", res[[2]], "\n")
}
