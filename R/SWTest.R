#' @title sw_test
#'
#' @description Performs a hypothesis test for the network small-world property
#'
#' @param x A square adjacency matrix.
#' @param model.list Vector of baseline models to test against. Possible values are ER, SBM, DCSBM and CL.
#' @param num.iter Vector of number of bootstrap iterations for each baseline model.
#'
#' @return Prints result of the test
#' @examples
#' data(dolphins)
#' sw_test(dolphins, c("ER", "SBM"), c(100, 50))
#' @export

sw_test <- function(x, model.list, num.iter){
  res.list <- vector(mode = "list", length(model.list))
  for(iter in 1:length(model.list)){
    if(startsWith(tolower(model.list[iter]), "e")){
      cat("\nTesting against ER\n")
      res.list[[iter]] <- small_world_ER(x, num.iter[iter])
    }
    else if(startsWith(tolower(model.list[iter]), "s")){
      cat("\nTesting against SBM\n")
      res.list[[iter]] <- small_world_SBM(x, num.iter[iter])
    }
    else if(startsWith(tolower(model.list[iter]), "d")){
      cat("\nTesting against DCSBM\n")
      res.list[[iter]] <- small_world_DCSBM(x, num.iter[iter])
    }
    else if(startsWith(tolower(model.list[iter]), "c")){
      cat("\nTesting against CL\n")
      res.list[[iter]] <- small_world_CL(x, num.iter[iter])
    }
    else{
      message("\n", model.list[iter], " is an invalid input\n")
      model.list <- model.list[-iter]
      num.iter <- num.iter[-iter]
    }
  }
  cat("Alternative hypothesis: Beta is not equal to 0 or 1\n")
  for(iter in 1:length(model.list)){
    sw_print_res(res.list[[iter]], model.list[iter])
  }
}

# Exact ER test
# Remove rproj file
# Change name?
