#' Calculates the percentage of variance explained for signals.
#' 
#' Calculates the percentage of variance explained for signals 
#' from mixtures of signals.
#' 
#' @param X Mixtures matrix containing mixtures of signals in 
#' columns, with column means removed (centred)
#' @param A Mixing matrix
#' @param S Separated signals in columns of matrix
#' @return List containing the percentage of variance explained
#' for each signal.
#' @seealso fastICA
#' @examples
#' ## Calculates PVE for two separated signals (see vignette)
#' n <- 100
#' s1 <- sin(0:(n - 1) / 2)
#' s2 <- -(15 - abs(0:(n - 1) %% (2 * 15) - 15)) / 5
#' S <- cbind(s1, s2)
#' A <- matrix(c(-0.141, -0.293, -0.301, 0.603), nrow=2, ncol=2)
#' X <- S %*% A
#' vcm <- matrix(NA, nrow=1, ncol=ncol(X))
#' for (i in 1:ncol(X)){
#'   vcm[i] <- mean(X[, i])
#'   X[, i] <- X[, i] - vcm[i]
#' }
#' ica_results <- fastICA::fastICA(X, n.comp=ncol(X), 
#'   alg.typ="parallel", method="C", fun="logcosh", verbose=TRUE, 
#'   maxit=200, tol=1e-04, alpha=1)
#' signal_pves <- calculate_pve_of_signals(X, ica_results$A, 
#'   ica_results$S, vcm)
#' print(signal_pves)
#' par(mfcol=c(2, 3))
#' for (i in 1:ncol(S)){
#'   plot(S[, i], type="l", main=paste("Artificial Signal", i), 
#'     xlab="", ylab="")
#' }
#' for (i in 1:ncol(X)){
#'   plot(X[, i], type="l", main=paste("Mixture", i), xlab="", 
#'     ylab="")
#' }
#' for (i in 1:ncol(ica_results$S)){
#'   plot(ica_results$S[, i], type="l", 
#'     main=paste("Separated Signal ", i, " PVE=", 
#'     round(signal_pves[i], 2), "%",  sep=""), xlab="", ylab="")
#' }
calculate_pve_of_signals <- function(X, A, S){
  original_var <- var(matrix(X))
  pve_list <- list()
  for (i in 1:ncol(S)){
    a <- matrix(A[i, ], nrow = 1, ncol = ncol(A))
    s <- matrix(S[, i], nrow = nrow(S), ncol = 1)
    var_s <- var(matrix( (s %*% a)))
    pve_s <- var_s / original_var * 100
    pve_list <- c(pve_list, pve_s)
  }
  pve.list <- as.numeric(pve_list)
  return(pve.list)
}
