#' Truncated Gaussian knockoff
#'
#' Construct Gaussian knockoff by truncating the diagonal variables
#'
#' @param x original data
#' @param Sigma covariance matrix of the original X
#' @param amp The amplitude of the gaussian copula, default value is set to 1, upper limit is 2.
#'
#' @return knockoff generated
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @keywords truncated, knockoff
#'
#' @examples
#'
#' library(mvtnorm)
#' p = 20
#' n = 100
#' u = eigen(cov(matrix(rnorm(p * n), nrow = n)))$vectors
#' d = runif(p, 1, 10)
#' Sigma = u %*% diag(d) %*% t(u)
#' x = rmvnorm(n, rep(0, p), Sigma)
#' truncatedGaussianKnockoff(x, Sigma, 1)
#' @export


truncatedGaussianKnockoff <- function(x, Sigma, amp = 1){
    p = nrow(Sigma)
    n = nrow(x)
    s = eigen(Sigma, symmetric = TRUE)$values[p] * amp
    A = 2 * diag(s, nrow = p) - diag(s, nrow = p) %*% solve(Sigma) %*% diag(s, nrow = p)
    C = chol(A)
    # U = eigen(svd(x)$u %*% t(svd(x)$u))$vectors[, (p + 1) : (2 * p)]
    U = svd((diag(n) - svd(x)$u %*% t(svd(x)$u)) %*% matrix(rnorm(nrow(x) * nrow(x)), nrow(x)))$u[, 1 : p]
    xHat = x %*% (diag(p) - solve(Sigma) %*% diag(s, nrow = p)) + U %*% C
    return(xHat)
}

