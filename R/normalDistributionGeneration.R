#' Normal Distribution Generation with Copula
#'
#' Generate normal distribution with provided copula of p features and n samples
#'
#' @param copula Copula used to generate normal distribution
#' @param p Number of features
#' @param n Number of samples
#'
#' @return A list of data, multivariate mean and multivariate standard deviation
#'
#' @keywords copula, normal distribution
#'
#' @examples
#'
#' out_copula = gaussian_copula(p = 100, n = 500, type = "diag")
#' out_data = normal_distribution_generate(copula = out_copula, p = 100, n = 500)
#'
#' @export

normal_distribution_generate <- function(copula, p, n){
    ## Normal distribution
    mean_list= list()
    sd_list = list()
    for(i in 1 : p){
        mean_list[[i]] = rnorm(1)
        sd_list[[i]] = rgamma(n = 1, shape = 2)
    }
    data.normal = c()
    for(i in 1 : p){
        data.normal = cbind(data.normal, unlist(lapply(copula$U[, i], function(x){qnorm(x, mean = mean_list[[i]], sd = sd_list[[i]])})))
    }
    return(list(data.normal = data.normal, mean_list = mean_list, sd_list = sd_list))
}
