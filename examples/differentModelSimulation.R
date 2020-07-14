cov_type = c("diag", "toeplitz", "general")
freq_list = list(freq_p = c(0.4, 0.6, 0.8), name_p  = c("04", "06", "08"))

p = 50
n = 200

out_table = c()

for(k in 1 : 3){
    for(j in 1 : 3){
        copula = gaussianCopula(p = p, n = n, type = cov_type[j])
        gaussian.list = gaussianDistributionGeneration(copula, p, n)
        gaussian.lasso = LassoGenerationModel(gaussian.list$data, freq_list$freq_p[k])
        X = gaussian.list$data
        y = gaussian.lasso$y

        set.seed(1001)
        ## Form 1: Knockoff with unknown Gaussian distribution
        knockoff.1 <- list(x.knockoff.copy = create.gaussian(X, mu = apply(X, 2, mean), Sigma = cov(X)))
        lasso.gaussian.1 <- lassoFDR(X, knockoff.1$x.knockoff.copy, y = y, mask = gaussian.lasso$mask, family = "gaussian")
        ridge.gaussian.1 <- ridgeFDR(X, knockoff.1$x.knockoff.copy, y = y, mask = gaussian.lasso$mask, family = "gaussian")
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "unknown Gaussian distribution", lasso.gaussian.1$power, lasso.gaussian.1$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "unknown Gaussian distribution", ridge.gaussian.1$power, ridge.gaussian.1$FDR, "Ridge"))

        set.seed(1001)
        ## Form 2: Knockoff with known Gaussian distribution
        knockoff.2 <- list(x.knockoff.copy = create.gaussian(X, mu = unlist(gaussian.list$mean_list), Sigma = diag(unlist(gaussian.list$sd_list)) %*% copula$R %*% diag(unlist(gaussian.list$sd_list))))
        lasso.gaussian.2 <- lassoFDR(X, knockoff.2$x.knockoff.copy, y = y, mask = gaussian.lasso$mask, family = "gaussian")
        ridge.gaussian.2 <- ridgeFDR(X, knockoff.2$x.knockoff.copy, y = y, mask = gaussian.lasso$mask, family = "gaussian")
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "known Gaussian distribution", lasso.gaussian.2$power, lasso.gaussian.2$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "known Gaussian distribution", ridge.gaussian.2$power, ridge.gaussian.2$FDR, "Ridge"))

        set.seed(1001)
        ## Form 3: Knockoff with known Gaussian copula and known marginal
        knockoff.3 <- copulaKnockoff(X, R = copula$R, marginal_list = gaussian.list)
        lasso.gaussian.3 <- lassoFDR(X, knockoff.3$x.knockoff.copy, y, gaussian.lasso$mask, family = "gaussian")
        ridge.gaussian.3 <- ridgeFDR(X, knockoff.3$x.knockoff.copy, y = y, mask = gaussian.lasso$mask, family = "gaussian")
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "known Gaussian copula and known marginal", lasso.gaussian.3$power, lasso.gaussian.3$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "known Gaussian copula and known marginal", ridge.gaussian.3$power, ridge.gaussian.3$FDR, "Ridge"))

        set.seed(1001)
        ## Form 4: Knockoff with known Gaussian copula and unknown marginal
        knockoff.4 <- copulaKnockoff(X, R = copula$R, marginal_list = NULL)
        lasso.gaussian.4 <- lassoFDR(X, knockoff.4$x.knockoff.copy, y, gaussian.lasso$mask, family = "gaussian")
        ridge.gaussian.4 <- ridgeFDR(X, knockoff.4$x.knockoff.copy, y = y, mask = gaussian.lasso$mask, family = "gaussian")
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "known Gaussian copula and unknown marginal", lasso.gaussian.4$power, lasso.gaussian.4$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "known Gaussian copula and unknown marginal", ridge.gaussian.4$power, ridge.gaussian.4$FDR, "Ridge"))

        set.seed(1001)
        ## Form 5: Knockoff with unknown Gaussian copula and known marginal
        knockoff.5 <- copulaKnockoff(X, R = NULL, marginal_list = gaussian.list)
        lasso.gaussian.5 <- lassoFDR(X, knockoff.5$x.knockoff.copy, y, gaussian.lasso$mask, family = "gaussian")
        ridge.gaussian.5 <- ridgeFDR(X, knockoff.5$x.knockoff.copy, y = y, mask = gaussian.lasso$mask, family = "gaussian")
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "unknown Gaussian copula and known marginal", lasso.gaussian.5$power, lasso.gaussian.5$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "unknown Gaussian copula and known marginal", ridge.gaussian.5$power, ridge.gaussian.5$FDR, "Ridge"))

        set.seed(1001)
        ## Form 6: Knockoff with estimated Gaussian copula
        knockoff.6 <- copulaKnockoff(X, R = NULL, marginal_list = NULL)
        lasso.gaussian.6 <- lassoFDR(X, knockoff.6$x.knockoff.copy, y, gaussian.lasso$mask, family = "gaussian")
        ridge.gaussian.6 <- ridgeFDR(X, knockoff.6$x.knockoff.copy, y = y, mask = gaussian.lasso$mask, family = "gaussian")
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "estimated Gaussian copula", lasso.gaussian.6$power, lasso.gaussian.6$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "estimated Gaussian copula", ridge.gaussian.6$power, ridge.gaussian.6$FDR, "Ridge"))

    }
}


out_table = print(as.data.frame(out_table))
colnames(out_table) = c("Covariance", "Sparsity", "Form", "Power", "FDR", "Model")
print(out_table)

