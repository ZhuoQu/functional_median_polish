packages <- c("mrfDepth", "ddalpha", "geometry")

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

hdepth <- mrfDepth::hdepth
options <- list(type = "Rotation",
                ndir = 750,
                approx = FALSE)

multi_depth <- function(mat, depth_name) {
   ### check whether all elements in a column in mat are zero
  nbvar <- ncol(mat)
   degeneracy_variable <- which(apply(mat, 2, var) == 0)
   if (length(degeneracy_variable) > 0 && ### efficient variables are at least 2
       length(degeneracy_variable) < nbvar - 1) {
     mat <- mat[, -degeneracy_variable]
   } else if (length(degeneracy_variable) >= nbvar - 1)
   {
     stop ("the efficient dimension of mat should be at least 2")
   }
   
   if (nrow(unique.matrix(mat)) <= nbvar + 1) { ### efficient observations are less than nbvar + 1
     func <- rep(0, nrow(mat)) 
   } else {
    if (depth_name == "MFHD") {
      options <- list(type = "Rotation", approx = TRUE)
      if (ncol(mat) <= 3 ||nrow(mat) < 120) {
        func <- hdepth(mat)$depthZ
      } else {
        func <- hdepth(mat, options = options)$depthZ
      }
    }
    if (depth_name == "MFPD") {
      func <- projdepth(mat)$depthZ
    }
    if (depth_name == "MFSPD") {
      func <- sprojdepth(mat)$depthZ
    }
    if (depth_name == "MFDPD") {
      func <- dprojdepth(mat)$depthZ
    } else if (depth_name == "MSBD") {
      func <- depth.simplicial(mat, mat, k = 0.05, exact = F)
    }
   }
  return (func)
}

multi_functional_depth <- function(x, depth_name, weight_type) { 
  #c("MFHD","MFPD","MFDPD","MFSPD","MSBD")
  ## x is a list of grid_size. In each list, it is a matrix of obs_nb * var_nb
  grid_size <- length(x)
  sample_size <- nrow(x[[1]])
  dp <- lapply(x, function(x_list) {
    multi_depth(x_list, depth_name)
    }) ###dp <- matrix(NA, nrow = t, ncol = n)
  
  depth_thres <- quantile(unlist(dp), 0.25)
  if (weight_type == "volume") {
    weight_bin <- as.numeric(sapply(1:grid_size, function(l) {
      data <- x[[l]]
      match_point <- data[which(dp[[l]] > depth_thres), ]
      return (try(convhulln(match_point, "FA")$vol, silent = TRUE))
    }))
    weight_bin[is.na(weight_bin)] <- 0
    weight_bin <- weight_bin / sum(weight_bin, na.rm = TRUE)
  } else {
    weight_bin <- rep(1 / grid_size, grid_size)
  }
  result <- sapply(1:sample_size, function(samp) {
      dp_col <- unlist(lapply(dp, function(k) {k[samp]}))
      weight_bin %*% dp_col   ###### dont understand what happened 
      ###### the length of dp is 365 but the length of dp_col is only 357
    })
  return (result)
}



