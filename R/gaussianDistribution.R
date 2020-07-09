#' Gaussian Distribution Generation
#'
#' Generate a Gaussian Distribution with specific covariance
#'
#' @param p Dimension of the covariance
#' @param n Number of samples
#' @param diag_cov Structure of the covariance matrix
#'
#' @importFrom clusterGeneration rcorrmatrix
#' @importFrom MASS mvrnorm
#'
#' @return A list contains correlation, simulated data, multivariate mean, multivariate standard deviation, covariance
#'
#' @keywords covariance, gaussian
#'
#' @examples
#' out = gaussian_distribution(100, 50, "unstructure)
#' out = gaussian_distribution(100, 50, "diagonal)
#' out = gaussian_distribution(100, 50, "toeplitz)
#'
#' @export

gaussian_distribution <- function(p, n, diag_cov = "unstructure"){
    if(diag_cov == "diagonal"){
        R = diag(p)
    }else if(diag_cov == "toeplitz"){
        R = toeplitz(0.8^(0:(p-1)))
    }else if(diag_cov == "unstructure"){
        R = rcorrmatrix(p)
    }else{
        R = rcorrmatrix(p)
    }
    mean_list = rnorm(p)
    sd_list = rgamma(n = p, shape = 2)
    X = mvrnorm(n, mu = mean_list, Sigma = diag(sd_list) %*% R %*% diag(sd_list))
    # X = t(chol(diag(sd_list) %*% R %*% diag(sd_list))) %*% matrix(rnorm(p * n), nrow = n) + mean_list
    return(list(R = R, X = X, mean_list = mean_list, sd_list = sd_list, Sigma = diag(sd_list) %*% R %*% diag(sd_list)))
}
