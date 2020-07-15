#' Data Generation with Lasso Model
#'
#' Generate data with the Lasso Model
#'
#' @param data_signal The covariate used to generate the data
#' @param p_freq The percentage of covariates that are set to be 0
#' @param signal_low_significant The lower bound on the signal significance, default value is set to be 3
#' @param signal_high_significant The higher bound on the signal significance, default value is set to be 5
#'
#' @import stats
#'
#' @return A list containing the data, covariant, the mask and p_freq
#'
#' @keywords lasso, generation
#'
#' @examples
#' out_X = gaussianDistribution(100, 50, "unstructure")
#' out_Y_Lasso = LassoGenerationModel(out_X$X, 0.5)
#'
#' @export

LassoGenerationModel <- function(data_signal, p_freq, signal_low_significant = 3, signal_high_significant = 5) {
    p = ncol(data_signal)
    p_covariate = runif(ncol(data_signal), signal_low_significant, signal_high_significant) * (rbinom(ncol(data_signal), 1, prob = 0.5) - 0.5) * 2
    mask = sample(0 : 1, p, replace = TRUE, prob = c(p_freq, 1 - p_freq))
    p_covariate[which(mask == 1)] = 0
    y = data_signal %*% p_covariate + matrix(rnorm(nrow(data_signal)), ncol = 1)
    return(list(y = y, beta = p_covariate, mask = mask, p_freq = p_freq))
}
