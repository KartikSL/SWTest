#' @title sw_test
#'
#' @description Provides an overview table for the time and scope conditions of
#'     a data set
#'
#' @param dat A data set object
#' @param id Scope (e.g., country codes or individual IDs)
#' @param time Time (e.g., time periods are given by years, months, ...)
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' @export

sw_test <- function(graph, baseline.model, num.iter){

  if(startsWith(baseline.model, "E")){
    res <- small_world_ER(graph, num.iter)
  }
  else if(startsWith(baseline.model, "S")){
    res <- small_world_SBM(graph, num.iter)
  }
  else if(startsWith(baseline.model, "D")){
    res <- small_world_DCSBM(graph, num.iter)
  }
  else if(startsWith(baseline.model, "C")){
    res <- small_world_CL(graph, num.iter)
  }
  else if(startsWith(baseline.model, "A")){
    res.er <- small_world_ER(graph, num.iter)
    res.sbm <- small_world_SBM(graph, num.iter)
    res.dcsbm <- small_world_DCSBM(graph, num.iter)
    res.cl <- small_world_CL(graph, num.iter)
  }
  else{

  }
}

# To do:
# Decide output format
# Add documentation for internal functions
# Add example for main function
