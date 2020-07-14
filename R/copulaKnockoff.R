#' Knockoff Construction through Copula
#'
#' Construct knockoffs by estimating Gaussian Copulas first
#'
#' @param data.x input data
#' @param R correlation matrix, can be provided, default value is NULL, where the algorithm will estimate it
#' @param marginal_list marginal list, containing mean and standard deviation
#'
#' @import stats
#' @importFrom knockoff create.gaussian
#'
#' @return A list of correlation and constructed knockoff
#'
#' @keywords copula, knockoff
#'
#' @examples
#' out_copula = gaussianCopula(p = 100, n = 500, type = "diag")
#' out_data = gaussianDistributionGeneration(copula = out_copula, p = 100, n = 500)
#' out = copulaKnockoff(out_data$data)
#'
#' @export

copulaKnockoff <- function(data.x, R = NULL, marginal_list = NULL){
    ## Estimate the empirical cumulative distribution
    ecdf.list = list()
    for(i in 1 : ncol(data.x)){
        ecdf.list[[i]] = ecdf(data.x[, i])
    }

    ## Estimate the Gaussian Copula
    u = matrix(0, nrow = nrow(data.x), ncol = ncol(data.x))
    x.knockoff.copy = matrix(0, nrow = nrow(data.x), ncol = ncol(data.x))

    if(is.null(marginal_list)){
        ## Latent variable
        for(j in 1 : ncol(data.x)){
            u[, j] = ifelse(is.infinite(qnorm(ecdf.list[[j]](data.x[, j]))), nrow(data.x) / (nrow(data.x) + 1), nrow(data.x) / (nrow(data.x) + 1) * qnorm(ecdf.list[[j]](data.x[, j])))
        }
        ## correlation of R
        if(is.null(R)){
            R = cor(u)
        }
        ## Generate the Gaussian knockoff
        u.knockoff.copy = create.gaussian(u, mu = 0, Sigma = R)
        prob.knockoff.copy = pnorm(u.knockoff.copy)
        for(j in 1 : ncol(data.x)){
            x.knockoff.copy[, j] = unname(quantile(data.x[, j], probs = prob.knockoff.copy[, j]))
        }
    }else{
        ## Latent variable
        for(j in 1 : ncol(data.x)){
            u[, j] = qnorm(pnorm(data.x[, j], mean = marginal_list$mean_list[[j]], sd = marginal_list$sd_list[[j]]))
        }
        ## correlation of R
        if(is.null(R)){
            R = cor(u)
        }
        ## Generate the Gaussian knockoff
        u.knockoff.copy = create.gaussian(u, mu = 0, Sigma = R)
        for(j in 1 : ncol(data.x)){
            x.knockoff.copy[, j] = marginal_list$sd_list[[j]] * u.knockoff.copy[, j] + marginal_list$mean_list[[j]]
        }
    }
    return(list(R = R, x.knockoff.copy = x.knockoff.copy))
}
