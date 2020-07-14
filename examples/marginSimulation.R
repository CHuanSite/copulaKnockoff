## Generate copulas and normal, t, gamma distributions for the data
## With Gaussian Copula to estimate the FDR
p = 50
n = 200

cov_type = c("diag", "toeplitz", "general")
freq_list = list(freq_p = c(0.4, 0.6, 0.8), name_p  = c("04", "06", "08"))

out_table = c()

for(k in 1 : 3){
    for(j in 1 : 3){
        freq_p = freq_list$freq_p[k]
        copula = gaussianCopula(p = p, n = n, type = cov_type[j])
        gaussian.list = gaussianDistributionGeneration(copula, p, n)
        t.list = tDistributionGeneration(copula, p, n)
        gamma.list = gammaDistributionGeneration(copula, p, n)

        set.seed(1001)
        ## Generate the data for lasso regression
        gaussian.lasso = LassoGenerationModel(gaussian.list$data.normal, freq_p)
        t.lasso = LassoGenerationModel(t.list$data.t, freq_p)
        gamma.lasso = LassoGenerationModel(gamma.list$data, freq_p)

        set.seed(1001)
        ## generate knockoff copies
        gaussian.knockoff = copulaKnockoff(gaussian.list$data)
        t.knockoff = copulaKnockoff(t.list$data.t)
        gamma.knockoff = copulaKnockoff(gamma.list$data)

        set.seed(1001)
        ## Fit lasso model with knockoff variables for gaussian, t and gamma
        lasso.gaussian.FDR = lassoFDR(gaussian.list$data, gaussian.knockoff$x.knockoff.copy, gaussian.lasso$y, gaussian.lasso$mask, family = "gaussian")
        lasso.t.FDR = lassoFDR(t.list$data, t.knockoff$x.knockoff.copy, t.lasso$y, t.lasso$mask, family = "gaussian")
        lasso.gamma.FDR = lassoFDR(gamma.list$data, gamma.knockoff$x.knockoff.copy, gamma.lasso$y, gamma.lasso$mask, family = "gaussian")

        ## Fit ridge model with knockoff variables for gaussian, t and gamma
        ridge.gaussian.FDR = lassoFDR(gaussian.list$data, gaussian.knockoff$x.knockoff.copy, gaussian.lasso$y, gaussian.lasso$mask, family = "gaussian")
        ridge.t.FDR = lassoFDR(t.list$data, t.knockoff$x.knockoff.copy, t.lasso$y, t.lasso$mask, family = "gaussian")
        ridge.gamma.FDR = lassoFDR(gamma.list$data, gamma.knockoff$x.knockoff.copy, gamma.lasso$y, gamma.lasso$mask, family = "gaussian")

        ## Output the table
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "Gaussian", lasso.gaussian.FDR$power, lasso.gaussian.FDR$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "t", lasso.t.FDR$power, lasso.t.FDR$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "Gamma", lasso.gamma.FDR$power, lasso.gamma.FDR$FDR, "Lasso"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "Gaussian", ridge.gaussian.FDR$power, ridge.gaussian.FDR$FDR, "Ridge"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "t", ridge.t.FDR$power, ridge.t.FDR$FDR, "Ridge"))
        out_table = rbind(out_table, c(cov_type[j], freq_list$freq_p[k], "Gamma", ridge.gamma.FDR$power, ridge.gamma.FDR$FDR, "Ridge"))


    }
}

out_table = print(as.data.frame(out_table))
colnames(out_table) = c("Covariance", "Sparsity", "Marginal Family", "Power", "FDR", "Model")
print(out_table)

## cov_type freq_list margin power FDR model
