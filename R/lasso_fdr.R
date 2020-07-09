#' Knockoff Variable Selection through Lasso while Controlling the FDR
#'
#' Select variable based on knockoff constructed from LASSO
#'
#' @param x Original data
#' @param x.knockoff knockoff copy of the original data
#' @param y Response variable
#' @param mask The variables that are masked from the purpose
#' @param family what kind of
#' @param fdr_rate Controlling the False Discovery Rate of the knockoff
#'
#' @keywords Lasso, knockoff,FDR
#'
#' @examples
#' out_copula = gaussian_copula(p = 100, n = 500, type = "diag")
#' out_data = normal_distribution_generate(copula = out_copula, p = 100, n = 500)
#' out = copula_knockoff(out_data)
#' lasso_out = Lasso_generation_model(data_signal = out_data, p_freq = 0.3)
#' lasso_fdr(x = out_data, x.knockoff = out$x.knockoff.copy, y = lasso_out$y)
#'
#' @export




lasso_fdr <- function(x, x.knockoff, y, mask = NULL, family = NULL, fdr_rate = 0.1){
    p = ncol(x)

    # swap = rbinom(ncol(x), 1, 0.5)
    # swap.M = matrix(swap, nrow = nrow(x), ncol = length(swap), byrow = TRUE)
    # x.swap = x * (1 - swap.M) + x.knockoff * swap.M
    # Xk.swap = x * swap.M + x.knockoff * (1 - swap.M)
    # cvfit.knockoff = cv.glmnet(normc(cbind(x.swap, Xk.swap)), y, family = family, alpha = 0.5)
    # Z = coef(cvfit.knockoff, s = "lambda.1se")[2 : (2 * p + 1)]
    # orig = 1 : p
    # W = abs(Z[orig]) - abs(Z[orig + p])
    # W = W * (1 - 2 * swap)

    X = cbind(x, x.knockoff)
    X = scale(X)
    # cvfit.knockoff = cv.glmnet(X, y, family = family, alpha = 1)
    cvfit.knockoff = linearRidge(y ~ X)
    # W = abs(coef(cvfit.knockoff, s = "lambda.1se")[2 : (p + 1)]) - abs(coef(cvfit.knockoff, s = "lambda.1se")[(p + 2) : (2 * p + 1)])
    W = abs(coef(cvfit.knockoff)[2 : (p + 1)]) - abs(coef(cvfit.knockoff)[(p + 2) : (2 * p + 1)])
    knockoff_res = knockoff.threshold(W, fdr = fdr_rate, offset = 1)
    selected_x = which(W >= knockoff_res)
    if(is.null(mask)){
        ## If no mask, then just return the selected covariate
        return(list(selected_x = selected_x, W = W, knockoff_res = knockoff_res))
    }else{
        ## If there exists mask, then erturn the power, FDR
        original_x = which(mask == 0)
        FDR = ifelse(length(selected_x) == 0, 0, 1 - sum(selected_x %in% original_x) / max(1, length(selected_x)))
        power = sum(selected_x %in% original_x) / max(1, length(original_x))
        return(list(selected_x = selected_x, original_x = original_x, FDR = FDR, power = power, W = W, knockoff_res = knockoff_res))
    }
}
