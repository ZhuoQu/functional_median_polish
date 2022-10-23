packages <- c("fdaoutlier", "rlist", "mrfDepth")

## Now load or install&load all
package_check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
source("multivariate_functional_depth.R")


two_way_median_polish_ufd <- function(Y, partition_row_index, partition_col_index) {
  ### 1. Y is a matrix, it represents univariate functional observations 
  ## with the number of rows = n (sample size) and the number of cols = the grid length
  ### 2. partition_index is a list of the disjoint indexes and their union is 1:n.
  ### It shows the possible categories of the effect and the length of partition_index >=2
  rematch_Y <- Y
  grand_effect <- rep(0, length.out = grid_size)
  row_effect <- matrix(0, nrow = length(partition_row_index), ncol = grid_size)
  col_effect <- matrix(0, nrow = length(partition_col_index), ncol = grid_size)
  add_row_effect <- row_effect + 1
  Iteration <- 1
  while (sum(add_row_effect) != 0) { ### criterion of the stop of the iteration
    cat(paste(Iteration, "\n"))
    partition_row_Y <- lapply(1:length(partition_row_index), function(k) {
    rematch_Y[partition_row_index[[k]], ]})
  
    median_row_effect <- lapply(partition_row_Y, function(l) {
      if (class(l) == "numeric") {
        median <- l ### allow only one observation in the subgroup
      } else {
      fast_mbd <- modified_band_depth(l) ## we use the fast modified band depth to obtain the median
      median <- l[order(fast_mbd, decreasing = TRUE)[1], ]
      }
      return (median)
    })
  
    if (length(median_row_effect) <= nbvar + 1) {
      median_median_row_effect <- colMeans(list.rbind(median_row_effect))
    } else {
      matrix_effect <- t(list.cbind(median_row_effect))
      fast_mbd_effect <- modified_band_depth(matrix_effect)
      median_median_row_effect <- matrix_effect[order(fast_mbd_effect, decreasing = TRUE)[1], ]
    }
  
    grand_effect <- grand_effect + median_median_row_effect ### the update of grand effect
  
    add_row_effect <- t(sapply(median_row_effect, function(k) {k - median_median_row_effect}))
    row_effect <- row_effect + add_row_effect ### the update of the row effect
  
    partition_row_Y <- lapply(1:length(partition_row_index), function(k) {
      t(apply(partition_row_Y[[k]], 1, function(rownum) {
        rownum - median_row_effect[[k]]
      })) 
    }) 
  
    rematch_Y[unlist(partition_row_index), ] <- list.rbind(partition_row_Y)
  
    partition_col_Y <- lapply(1:length(partition_col_index), function(k) {
      rematch_Y[partition_col_index[[k]], ]}) 
  
    median_col_effect <- lapply(partition_col_Y, function(l) {
      if (class(l) == "numeric") {
        median <- l ### allow only one observation in the subgroup
      } else {
      fast_mbd <- modified_band_depth(l) ## we use the fast modified band depth to obtain the median
      median <- l[order(fast_mbd, decreasing = TRUE)[1], ]
      }
      return (median)
    })
  
    if (length(median_col_effect) <= nbvar + 1) {
      median_median_col_effect <- colMeans(list.rbind(median_col_effect))
    } else {
      matrix_effect <- t(list.cbind(median_col_effect))
      fast_mbd_effect <- modified_band_depth(matrix_effect)
      median_median_col_effect <- matrix_effect[order(fast_mbd_effect, decreasing = TRUE)[1], ]
    }
  
    grand_effect <- grand_effect + median_median_col_effect ### the update of grand effect
  
    add_col_effect <- t(sapply(median_col_effect, function(k) {k - median_median_col_effect}))
    col_effect <- col_effect + add_col_effect ### the update of the column effect
  
    partition_col_Y <- lapply(1:length(partition_col_index), function(k) {
      t(apply(partition_col_Y[[k]], 1, function(rownum) {
        rownum - median_col_effect[[k]]
      })) 
    }) 

    rematch_Y[unlist(partition_col_index), ] <- list.rbind(partition_col_Y)
    Iteration <- Iteration + 1
  }
  return(list(grand_effect = grand_effect,
              row_effect = row_effect,
              col_effect = col_effect))
}