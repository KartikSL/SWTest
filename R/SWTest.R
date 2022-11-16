#' @title sw_test
#'
#' @description Performs a hypothesis test for the network small-world property
#'
#' @param graph An input network
#' @param baseline.model Model to test against
#' @param num.iter Number of bootstrap iterations
#'
#' @return Prints result of the test
#' @examples
#' data(dolphins)
#' sw_test(dolphins, "ER", 100)
#' @export

sw_test <- function(graph, baseline.model, num.iter){

  if(startsWith(tolower(baseline.model), "e")){
    res <- small_world_ER(graph, num.iter)
  }
  else if(startsWith(tolower(baseline.model), "s")){
    res <- small_world_SBM(graph, num.iter)
  }
  else if(startsWith(tolower(baseline.model), "d")){
    res <- small_world_DCSBM(graph, num.iter)
  }
  else if(startsWith(tolower(baseline.model), "c")){
    res <- small_world_CL(graph, num.iter)
  }
  else if(startsWith(tolower(baseline.model), "a")){
    res.er <- small_world_ER(graph, num.iter)
    res.sbm <- small_world_SBM(graph, num.iter)
    res.dcsbm <- small_world_DCSBM(graph, num.iter)
    res.cl <- small_world_CL(graph, num.iter)
    res.out <- list(res.er, res.sbm, res.dcsbm, res.cl)
    model.list <- list("ER", "SBM", "DCSBM", "Chung-Lu")
    for(i in 1:4){
      sw_print_res(res.out[[i]], model.list[[i]])
    }
  }
  else{
    cat("Invalid input/n")
  }
}


# Add messages in functions
# Add documentation for all functions
