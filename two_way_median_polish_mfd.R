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

two_way_median_polish_mfd <- function(Y, partition_row_index, partition_col_index) {
  ### 1. when Y is an array, it represents multivariate functional observation, and [, , i]
  ## represents its ith marginal observations
  ### 2. partition_index is a list of the disjoint indexes and their union is 1:n.
  ### It shows the possible categories of the effect and the length of partition_index >=2
  rematch_Y <- Y
  grand_effect <- matrix(0, nrow = grid_size, ncol = nbvar)
  row_effect <- lapply(1:length(partition_row_index), function(k) {
    matrix(0, nrow = grid_size, ncol = nbvar)
  })
  col_effect <- lapply(1:length(partition_col_index), function(k) {
    matrix(0, nrow = grid_size, ncol = nbvar)
  })

  add_row_effect <- lapply(row_effect, function(k) {k + 1})
  Iteration <- 1


####### Iteration of the two-way median polish
  while (sum(add_row_effect[[1]])!= 0) {
    cat(paste(Iteration, "\n"))
    partition_row_Y <- lapply(1:length(partition_row_index), function(k) {
      rematch_Y[partition_row_index[[k]], , ]})
  
    median_row_effect <- lapply(partition_row_Y, function(l) {
      
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
  
    if (length(median_row_effect) <= nbvar + 1) {
      median_median_row_effect <- sapply(1:nbvar, function(k) {
      rowMeans(sapply(median_row_effect, function(subset) {
        subset[, k]
      }))
      })  ### median_median_row_effect is a matrix with grid_size rows and nbvar cols
    } else {
      list_data <- lapply(1:grid_size, function(time_point) {
        t(sapply(median_row_effect, function(sub) {
        sub[time_point, ]
        }))
      })
    fd_process <- multi_functional_depth(list_data, 
                                         depth_name = "MFHD", ##multivariate functional halfspace depth
                                         weight_type = "volume") ### we use multivariate functional halfspace depth
    med_index <- order(fd_process, decreasing = TRUE)[1]
    median_median_row_effect <- median_row_effect[[med_index]]
    }
  
  
  grand_effect <- grand_effect + median_median_row_effect ### the update of grand effect
  
  add_row_effect <- lapply(median_row_effect, function(k) {k - median_median_row_effect})
  ## effect is a list of the number of effect categories. 
  ## In each sublist, it is a matrix of grid_size rows and nbvar columns
  
  for (k in 1:length(partition_row_index)) {
    row_effect[[k]] <- row_effect[[k]] + add_row_effect[[k]] ### the update of row effect
  }
  
  #  plot(1:grid_size, grand_effect[, 1])
  #  lines(1:grid_size, grand_effect[, 1] + row_effect[[1]][, 1], col = "yellow")
  #  lines(1:grid_size, grand_effect[, 1] + row_effect[[2]][, 1], col = "yellow")
  # for (k in 1:length(partition_row_index)) {
  #   for (var in 1:nbvar) {
  #      temp <-  apply(partition_row_Y[[k]][, , var], 1, function(rownum) {
  #         rownum - median_row_effect[[k]][, var]})
  #      partition_row_Y[[k]][, , var] <- temp
  #     }
  #   }
  
  rematch_Y[unlist(partition_row_index), , ] <- abind(partition_row_Y, along = 1)
  
  # par(mfrow = c(1, 3))
  # for (var in 1:nbvar) {
  #   plot(1:grid_size, rematch_Y[1, , var], ylim = range(rematch_Y[, , var]))
  #   sapply(1:sample_size, function(k) {
  #     lines(1:grid_size, rematch_Y[k, , var])
  #   })
  # }
  
  partition_col_Y <- lapply(1:length(partition_col_index), function(k) {
    rematch_Y[partition_col_index[[k]], , ]}) 
  
  median_col_effect <- lapply(partition_col_Y, function(l) {
    
    if (length(class(l)) == 2 && class(l)[1] == "matrix") {
      median <- l
    } else  {    
    list_data <- lapply(1:grid_size, function(time_point) {
      l[, time_point, ]
    })
    fd_process <- multi_functional_depth(list_data, 
                                         depth_name = "MFHD", ##multivariate functional halfspace depth
                                         weight_type = "volume") ### we use multivariate functional halfspace depth
    median <- l[order(fd_process, decreasing = TRUE)[1], , ]
    }
  })
  
  if (length(median_col_effect) <= nbvar + 1) {
    median_median_col_effect <- sapply(1:nbvar, function(k) {
      rowMeans(sapply(median_col_effect, function(subset) {
        subset[, k]
      }))
    })  ### median_median_col_effect is a matrix with grid_size rows and nbvar cols
  } else {
    list_data <- lapply(1:grid_size, function(time_point) {
      t(sapply(median_col_effect, function(sub) {
        sub[time_point, ]
      }))
    })
    fd_process <- multi_functional_depth(x = list_data, 
                                         depth_name = "MFHD", ##multivariate functional halfspace depth
                                         weight_type = "volume") ### we use multivariate functional halfspace depth
    med_index <- order(fd_process, decreasing = TRUE)[1]
    median_median_col_effect <- median_col_effect[[med_index]]
  }
  
  grand_effect <- grand_effect + median_median_col_effect ### the update of grand effect
  
  add_col_effect <- lapply(median_col_effect, function(k) {
    k - median_median_col_effect})
  
  for (k in 1:length(partition_col_index)) {
    col_effect[[k]] + add_col_effect[[k]]
  } ### the update of the column effect
  
  for (k in 1:length(partition_col_index)) {
    for (var in 1:nbvar) {
      temp <- t(apply(partition_col_Y[[k]][, , var], 1, function(v) {
        v - median_col_effect[[k]][, var]
      }))
      partition_col_Y[[k]][, , var] <- temp
    }
  }
  rematch_Y[unlist(partition_col_index), , ] <- abind(partition_col_Y, along = 1)
  Iteration <- Iteration + 1
  }
  return (list(grand_effect = grand_effect,
               row_effect = row_effect,
               col_effect = col_effect))
}