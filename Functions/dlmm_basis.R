# Distributed lag interaction (mixture) model basis function

# x = exposures
# basis_all = basis generating info provided to all exposures
# basis_each = list of basis information with the same length as x supplied individually (uniquely) to each exposure
# lag = list containing the number of lag days for each exposure
# interactions = type of interactions to examine - "none", "noself", or "all"

dlmm_basis <- function(x = list(),
                       basis_all = list(),
                       basis_each = list(),
                       lag = list(),
                       interactions = "noself") {
  
  if (!is.list(x) | !is.list(basis_all) | !is.list(basis_each) | !is.list(lag))
    stop("all inputs must be of type list")
  
  M <- length(x)
  x_names <- names(x)
  
  if (is.null(x_names)) {
    names(x) <- paste0("x", 1:M)
    x_names <- names(x)
  }
  
  if (!all(sapply(x, is.matrix)))
    stop("all items in list `x` must be of type `matrix`")
  if (!all(sapply(x, nrow) == nrow(x[[1]])))
    stop("rows of all matrices in list `x` must be equal")
  
  N <- nrow(x[[1]])
  
  # Check basis
  if (length(basis_all) != 0) {
    basis_each <- lapply(1:M, function(i) basis_all)
    names(basis_each) <- x_names
  } else if (length(basis_each) != M) {
    stop("if not using basis_all, must specify basis for each x in list")
  } else {
    if (!is.null(colnames(basis_each)) && colnames(basis_each) != x_names)
      stop("names in basis_each must match names in x")
  }
  
  names(basis_each) <- x_names
  
  # Check lag
  if (length(lag) == 0) {
    lag <- lapply(x, function(i) c(0, ncol(i) - 1))
  } else if (length(lag) != M) {
    stop("length of list lag must equal length of list x")
  }
  
  names(lag) <- x_names

  # create basis functions
  seqlag <- function(range) seq(range[1], range[2], by = 1)
  
  b <- list()
  w <- list()
  
  for (i in 1:M) {
    n <- x_names[i]
    basis_each[[n]]$x <- seqlag(lag[[n]])
    basis_each[[n]]$intercept <- TRUE
    b[[n]] <- do.call(onebasis, basis_each[[n]])
    w[[n]] <- x[[n]] %*% b[[n]]
    colnames(w[[n]]) <- paste0("e", i, ".", 1:ncol(b[[n]]))
  }
  
  # Create crossbasis
  if (length(x) > 1 & interactions != "none") {
    
    cb <- list()
    
    for (i in 1:ifelse(interactions == "noself", M-1, M)) {
      for (j in ifelse(interactions == "noself", i+1, i):M) {
        n1 <- x_names[i]
        n2 <- x_names[j]
        n <- paste0(n1, "-", n2)
        cb[[n]] <- b[[n1]] %x% b[[n2]]
        w[[n]] <- t(sapply(1:N, function(k) {
          crossprod(cb[[n]], (x[[n1]][k,] %x% x[[n2]][k,])) }))
        colnames(w[[n]]) <-
          paste0(paste0("m", i, ".", rep(1:ncol(b[[n1]]), each = ncol(b[[n2]])), "."),
                 paste0("m", j, ".", rep(1:ncol(b[[n2]]), ncol(b[[n1]]))))
      }
    }
    
  }

  w <- do.call(cbind, w)

  # Save attributes
  attributes(w)$basis <- basis_each
  attributes(w)$lag <- lag
  attributes(w)$x_names <- x_names
  attributes(w)$interactions <- interactions

  class(w) <- c("dlmm_basis", "matrix")
  return(w)
  
}
