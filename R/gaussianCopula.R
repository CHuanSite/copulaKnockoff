#' Gaussian Copula Generation
#'
#' Generate a Gaussian copula with specific covariance structures
#'
#' @param p Dimension of the covariance
#' @param n Number of samples
#' @param type Structure of the covariance matrix, has option "diag", "toeplitz" and "unstructured"
#'
#' @importFrom clusterGeneration rcorrmatrix
#' @importFrom MASS mvrnorm
#'
#' @return A list of correlation anddata
#'
#' @keywords Gaussian, Copula
#'
#' @examples out = gaussian_copula(1000, 100, "diag")
#'
#'
#' @export

gaussian_copula <- function(p, n, type = "diag"){
    ## Random correlation matrix
    if(type == "diag"){
        R = diag(p)
    } else if (type == "toeplitz") {
        R = toeplitz(0.8^(0:(p-1)))
    } else{
        R = rcorrmatrix(p)
    }
    A = t(chol(R))
    X = c()
    for(i in 1 : n){
        X = rbind(X, mvrnorm(1, mu = rep(0, p), Sigma = R))
    }
    U = pnorm(X)

    return(list(U = U, R = R, X = X))
}
