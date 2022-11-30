#' Conducts bootstrap test against ER
#' @param x adjacency matrix
#' @param b number of bootstrap iterations
#' @keywords internal
#' @noRd
small_world_ER <- function(x, b){
  G <- igraph::graph_from_adjacency_matrix(x, mode = "undirected")
  n <- igraph::gorder(G)
  L <- igraph::mean_distance(G)
  C <- igraph::transitivity(G, type = "global")
  phat <- igraph::gsize(G)/choose(n,2)

  tstar <- rep(NA, b)
  Lstar <- rep(NA, b)

  for (iter in 1:b){
    Gstar <- igraph::sample_gnp(n = n, p = phat)
    Lstar[iter] <- igraph::mean_distance(Gstar)
    tstar[iter] <- igraph::transitivity(Gstar, type = "global")
    if(iter%%50 == 0){
      cat("Iteration ", iter, "\n")
    }
  }

  pval1 <- sum(tstar > C)/b
  pval2 <- sum(Lstar < L)/b

  return(list(pval1 = pval1, pval2 = pval2,
              coef1 = C, coef2 = L))
}

#' Conducts bootstrap test against SBM
#' @param x adjacency matrix
#' @param b number of bootstrap iterations
#' @keywords internal
#' @noRd
small_world_SBM <- function(x, b){

  G <- igraph::graph_from_adjacency_matrix(x, mode = "undirected")
  n <- igraph::gorder(G)
  L <- igraph::mean_distance(G)
  C <- igraph::transitivity(G, type = "global")
  A <- igraph::as_adj(G, type = "both", sparse = FALSE)

  k <- length(igraph::cluster_louvain(G))
  clusters <- dcspectral(igraph::as_adj(G, type = "both", sparse = FALSE), n, k)

  ones <- matrix(1, n, n)
  diag(ones) <- 0
  phat <- matrix(0, k, k)
  for (i in 1:k){
    for ( j in i:k){
      phat[i, j] <- sum(A[clusters == i,clusters == j])/sum(ones[clusters == i, clusters == j])
      phat[j, i] <- phat[i,j]
    }
  }

  tstar <- rep(NA,b)
  Lstar <- rep(NA,b)

  P <- get_P_SBM(n, k, phat, clusters)

  for (iter in 1:b){
    Gstar <- sample_from_SBM(n, P)
    Lstar[iter] <- igraph::mean_distance(Gstar)
    tstar[iter] <- igraph::transitivity(Gstar, type = "global")
    if(iter%%50 == 0){
      cat("Iteration ",iter, "\n")
    }
  }

  pval1 = sum(tstar > C)/b
  pval2 = sum(Lstar < L)/b

  return(list(pval1 = pval1, pval2 = pval2,
              coef1 = C, coef2 = L))

}

#' Conducts bootstrap test against DCSBM
#' @param x adjacency matrix
#' @param b number of bootstrap iterations
#' @keywords internal
#' @noRd
small_world_DCSBM <- function(x, b){

  G <- igraph::graph_from_adjacency_matrix(x, mode = "undirected")
  n <- igraph::gorder(G)
  L <- igraph::mean_distance(G)
  C <- igraph::transitivity(G, type = "global")
  A <- as.matrix(igraph::as_adj(G, type = "both"))

  k <- length(igraph::cluster_louvain(G))
  clusters <- dcspectral(igraph::as_adj(G, type = "both", sparse = FALSE), n, k)

  thetahat <- rep(0, n)
  thetau <- rowSums(A)
  phat <- matrix(0, k, k)

  for (i in 1:k){
    sizek <- length(thetau[clusters == i])
    thetak <- sum(thetau[clusters == i])
    thetahat[clusters == i] <- sizek*(thetau[clusters == i])/thetak
  }

  for (i in 1:k){
    for ( j in i:k){
      phat[i, j] <- sum(A[clusters == i, clusters == j])/length(A[clusters == i, clusters == j])
      phat[j, i] <- phat[i, j]
    }
  }

  tstar <- rep(NA, b)
  Lstar <- rep(NA, b)

  P <- get_P_DCSBM(n, k, thetahat, phat, clusters)

  for (iter in 1:b){
    Gstar <- sample_from_DCSBM(n, P)
    Lstar[iter] <- igraph::mean_distance(Gstar)
    tstar[iter] <- igraph::transitivity(Gstar, type = "global")
    if(iter%%50 == 0){
      cat("Iteration ",iter, "\n")
    }
  }

  pval1 <- sum(tstar > C)/b
  pval2 <- sum(Lstar < L)/b

  return(list(pval1 = pval1,pval2 = pval2,
              coef1 = C, coef2 = L))

}

#' Conducts bootstrap test against CL
#' @param x adjacency matrix
#' @param b number of bootstrap iterations
#' @keywords internal
#' @noRd
small_world_CL <- function(G, b){

  G <- igraph::graph_from_adjacency_matrix(x, mode = "undirected")
  n <- igraph::gorder(G)
  d <- igraph::degree(G)
  P <- d %*%t (d)/sum(d)

  C <- igraph::transitivity(G, type = "global")
  L <- igraph::mean_distance(G)

  tstar <- rep(NA, b)
  Lstar <- rep(NA, b)
  tLstar <- rep(NA, b)

  for(iter in 1:b){
    Gstar <- sample_from_CL(P)
    Lstar[iter] <- cl_dist(Gstar)
    tstar[iter] <- igraph::transitivity(Gstar, type = "global")
    if(iter%%50 == 0){
      cat("Iteration ",iter, "\n")
    }
  }

  pval1 <- sum(tstar > C)/b
  pval2 <- sum(Lstar < L)/b

  return(list(pval1 = pval1, pval2 = pval2,
              coef1 = C, coef2 = L))
}
