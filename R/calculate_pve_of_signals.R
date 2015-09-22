#' Calculates the percentage of variance explained for signals.
#' 
#' Calculates the percentage of variance explained for signals from mixtures
#' of signals.
#' 
#' @param X Mixtures matrix containing mixtures of signals in columns. Must have
#' the column means removed (centred)
#' @param A Mixing matrix
#' @param S Separated signals in columns of matrix
#' @param vcm vector containing column means of X 
#' @return List containing the percentage of variance explained for each signal.
#' @section Warning: There is currently no known formal proof 
#' for this calculation.
#' @seealso fastICA
#' @examples
#' (see vignette for worked example) 
calculate_pve_of_signals = function(X, A, S, vcm){
  original_var <- var(matrix(X))
  pve_list <- list()
  for (i in 1:ncol(S)){
    a <- matrix(A[i, ], nrow = 1, ncol = ncol(A))
    s <- matrix(S[ ,i], nrow = nrow(S), ncol = 1)
    var_s <- var(matrix((s %*% a))  + (vcm[i] / ncol(S)))
    pve_s <- var_s / original_var * 100
    pve_list <- c(pve_list, pve_s)
  }
  pve.list <- as.numeric(pve_list)
  return(pve.list)
}
