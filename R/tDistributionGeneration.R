#' t Distribution Generation with Copula
#'
#' Generate t distribution with provided copula of p features and n samples
#'
#' @param copula Copula used to generate t distribution
#' @param p Number of features
#' @param n Number of samples
#'
#' @return A list of data, multivariate mean and multivariate standard deviation
#'
#' @keywords copula, t distribution
#'
#' @examples
#'
#' out_copula = gaussianCopula(p = 100, n = 500, type = "diag")
#' out_data = tDistributionGeneration(copula = out_copula, p = 100, n = 500)
#'
#' @export

tDistributionGeneration <- function(copula, p, n){
    ## t distribution
    df_list = list()
    for(i in 1 : p){
        df_list[[i]] = sample(5 : 20, 1)
    }
    data.t = c()
    for(i in 1 : p){
        data.t = cbind(data.t, unlist(lapply(copula$U[, i], function(x){qt(x, df_list[[i]])})))
    }
    return(list(data.t = data.t, df_list = df_list))
}
