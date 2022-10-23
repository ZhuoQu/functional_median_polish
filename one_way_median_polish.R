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
source("multi_functional_depth.R")

one_way_median_polish <- function(Y, partition_index) {
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
  
  if (length(partition_index) < 2) {
    stop("partition_index is not a qualified partition")
  }
  
  if (nbvar == 1) {  ### implementation for the univariate functional data
    partition_Y <- lapply(1:length(partition_index), function(k) {
      Y[partition_index[[k]], ]})
    
    
    iterative_median <- lapply(partition_Y, function(l) {
      if (class(l) == "numeric") {
        median <- l ### allow only one observation in the subgroup
      } else {
      fast_mbd <- modified_band_depth(l) ## allow repetatives in subgroups
      ###we use the fast modified band depth to obtain the median in each subggroup
      median <- l[order(fast_mbd, decreasing = TRUE)[1], ]
      } 
      return (median)
    })
    
    if (length(iterative_median) <= nbvar + 1) {
      grand_effect <- colMeans(list.rbind(iterative_median))
    } else {
      matrix_effect <- t(list.cbind(iterative_median))
      fast_mbd_effect <- modified_band_depth(matrix_effect) 
      ## use fast modified band depth if the length of the effect is larger than nbvar + 1
      grand_effect <- matrix_effect[order(fast_mbd_effect, decreasing = TRUE)[1], ]
    }
    ### the update of grand effect
    
    row_effect <- t(sapply(iterative_median, function(k) {k - grand_effect}))
     ### the update of the row effect
      
    # partition_Y <- lapply(1:length(partition_Y), function(k) {
    #   t(apply(partition_Y[[k]], 1, function(rownum) {rownum - iterative_median[[k]]}))
    # })
    
  } else { ### implementation for the multivariate functional data
    partition_Y <- lapply(1:length(partition_index), function(k) {
      Y[partition_index[[k]], , ]})
    
    iterative_median <- lapply(partition_Y, function(l) {
      if (length(class(l)) == 2 && class(l)[1] == "matrix") {
        median <- l
      } else  {    
        list_data <- lapply(1:grid_size, function(time_point) {
        l[, time_point, ]
      })
      
      fd_process <- multi_functional_depth(x = list_data, 
                                           depth_name = "MFHD", ##multivariate functional halfspace depth
                                           weight_type = "volume") ### we use multivariate functional halfspace depth
      median <- l[order(fd_process, decreasing = TRUE)[1], , ]
      }
    })
    
    
    if (length(iterative_median) <= nbvar + 1) {
      grand_effect <- sapply(1:nbvar, function(k) {
        rowMeans(sapply(iterative_median, function(subset) {
          subset[, k]
        }))
      })  ### iterative_median is a matrix with grid_size rows and nbvar cols
    } else {
      list_data <- lapply(1:grid_size, function(time_point) {
        t(sapply(iterative_median, function(sub) {
          sub[time_point, ]
        }))
      })
      fd_process <- multi_functional_depth(list_data, 
                                           depth_name = "MFHD", ##multivariate functional halfspace depth
                                           weight_type = "volume") ### we use multivariate functional halfspace depth
      med_index <- order(fd_process, decreasing = TRUE)[1]
      grand_effect <- iterative_median[[med_index]]
    } ### obtain grand effect
    
    row_effect <- lapply(iterative_median, function(k) {k - grand_effect})
    ### obtain row effect
  }
  return(list(grand_effect = grand_effect,
              row_effect = row_effect))
}


