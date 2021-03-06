#' Gamma Distribution Generation with Copula
#'
#' Generate gamma distribution with provided copula of p features and n samples
#'
#' @param copula Copula used to generate gamma distribution
#' @param p Number of features
#' @param n Number of samples
#'
#' @return A list of data, multivariate mean and multivariate standard deviation
#'
#' @keywords copula, gamma distribution
#'
#' @examples
#'
#' out_copula = gaussianCopula(p = 100, n = 500, type = "diag")
#' out_data = gammaDistributionGeneration(copula = out_copula, p = 100, n = 500)
#'
#' @export

gammaDistributionGeneration <- function(copula, p, n){
    ## Gamma distribution
    shape_list = list()
    rate_list = list()
    for(i in 1 : p){
        rate_list[[i]] = sample(2:10, 1)
        shape_list[[i]] = sample(2:10, 1)
    }
    data = c()
    for(i in 1 : p){
        data = cbind(data, unlist(lapply(copula$U[, i], function(x){qgamma(x, shape = shape_list[[i]], rate = rate_list[[i]])})))
    }

    return(list(data = data, rate_list = rate_list, shape_list = shape_list))
}
