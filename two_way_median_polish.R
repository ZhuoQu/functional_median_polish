source("two_way_median_polish_ufd.R")
source("two_way_median_polish_mfd.R")

two_way_median_polish <- function(Y, partition_row_index, partition_col_index) {
  ### 1. when Y is a matrix, it represents univariate functional observations 
  ## with the number of rows = n (sample size) and the number of cols = the grid length
  ### when Y is an array, it represents multivariate functional observation, and [, , i]
  ## represents its ith marginal observations
  ### 2. partition_index is a list of the disjoint indexes and their union is 1:n.
  ### It shows the possible categories of the effect and the length of partition_index >=2
  if (class(Y) != "matrix" && class(Y) != "array") {
    stop("Y should be functional data")
  } else if (length(class(Y)) == 2 && 
             class(Y)[1] == "matrix" && 
             class(Y)[2] == "array") {
    nbvar <- 1
    sample_size <- nrow(Y)
    grid_size <- ncol(Y)
  } else if (length(class(Y)) == 1 && class(Y) != "matrix") {
    nbvar <- dim(Y)[3]
    sample_size <- dim(Y)[1]
    grid_size <- dim(Y)[2]
  }
  
  if (length(partition_row_index) < 2 || length(partition_col_index) < 2 ) {
    stop("partition_index is not a qualified partition")
  }
  
  if (nbvar == 1) {
    result <- two_way_median_polish_ufd(Y, partition_row_index, partition_col_index)
  } else { ### implementation for the multivariate functional data
    result <- two_way_median_polish_mfd(Y, partition_row_index, partition_col_index)
  }
  return (result)
}
  



