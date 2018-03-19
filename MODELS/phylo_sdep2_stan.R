# Marina Golivets
# Last modofied: 10/25/2017

meta.stan <- function (data, X, VCV_nb, VCV_tg, sigma2_d, mm_nb, mm_tg, niter, warmup, thin, chains) {

    modelString <- "
data {
    int <lower=0> N; // num obs
    int <lower=0> J; // num studies
    int <lower=1, upper=J> jj[N]; // study index
    int <lower=1> K; // num sampling groups
    int<lower=1, upper=K> kk[N]; // sampling group index
    int <lower=1> M_nb; // num neighbor species
    int <lower=1> M_tg; // num target species
    int <lower=1, upper=M_nb> mm_nb[N]; // neighbor species index
    int <lower=1, upper=M_tg> mm_tg[N]; // target species index
    vector[M_nb] I_nb; // vector of 1s for neighbor species
    vector[M_tg] I_tg; // vector of 1s for target species
    corr_matrix[M_nb] VCV_nb; // phylogenetic correlation matrix for neighbor species
    corr_matrix[M_tg] VCV_tg; // phylogenetic correlation matrix for target species
    int <lower=1> p; // num vars in design matrix
    matrix [N, p] X; // design matrix
    vector [N] d; // Hedges' d estimates
    vector[N] sigma_d;
    real <lower=0> sigma2_d; // total sampling variance 
}
transformed data {
    cholesky_factor_corr[M_nb] VCV_nb_chol;
    cholesky_factor_corr[M_tg] VCV_tg_chol;
    VCV_nb_chol = cholesky_decompose(VCV_nb); // Cholesky decomposed phylogenetic correlation matrix for neighbor species
    VCV_tg_chol = cholesky_decompose(VCV_tg); // Cholesky decomposed phylogenetic correlation matrix for target species
}
parameters {
    vector [p] beta;
    real <lower=0> sigma_sdep;
    real <lower=0> sigma_study; 
    real <lower=0> sigma_obs; 
    real <lower=0> sigma_nb;
    real <lower=0> sigma_tg;
    real <lower=0> sigma_nb_phylo; 
    real <lower=0> sigma_tg_phylo;
    vector[K] eta_sdep;
    vector[J] eta_study; 
    vector[N] eta_obs;
    vector[M_nb] eta_nb;
    vector[M_tg] eta_tg;
    vector[M_nb] eta_nb_phylo; 
    vector[M_tg] eta_tg_phylo;
}
transformed parameters {
    vector[N] mu;
    vector[K] eff_sdep; // sapling dependency group effect
    vector[J] eff_study; // study effect
    vector[M_nb] eff_nb; // neighbor species effect
    vector[M_tg] eff_tg; // target species effect
    vector[M_nb] sigma_nb_phylo_v;
    vector[M_tg] sigma_tg_phylo_v;
    vector[M_nb] eff_nb_phylo; // neighbor species phylogenetic effect
    vector[M_tg] eff_tg_phylo; // target species phylogenetic effect
    vector[N] theta;
    //real <lower=0> prec_nb_phylo;
    //real <lower=0> prec_tg_phylo;
    mu = X * beta;
    eff_sdep = sigma_sdep * eta_sdep;
    eff_study = sigma_study * eta_study;
    eff_nb = sigma_nb * eta_nb;
    eff_tg = sigma_tg * eta_tg;
    sigma_nb_phylo_v = sigma_nb_phylo * I_nb;
    sigma_tg_phylo_v = sigma_tg_phylo * I_tg;
    eff_nb_phylo = sigma_nb_phylo_v .* (VCV_nb_chol * eta_nb_phylo);
    eff_tg_phylo = sigma_tg_phylo_v .* (VCV_tg_chol * eta_tg_phylo);
    for (i in 1:N) theta[i] = mu[i] + eff_sdep[kk[i]] + eff_study[jj[i]] + eff_nb[mm_nb[i]] + eff_tg[mm_tg[i]] + eff_nb_phylo[mm_nb[i]] + eff_tg_phylo[mm_tg[i]] + sigma_obs * eta_obs[i];
    //prec_nb_phylo = 1/sigma_nb_phylo;
    //prec_tg_phylo = 1/sigma_tg_phylo;
}
model {
    d ~ normal(theta, sigma_d);
    // priors
    beta ~ normal(0, 5);
    sigma_sdep ~ normal(0, 5);
    sigma_study ~ normal(0, 5);
    sigma_obs ~ normal(0, 5);
    sigma_nb ~ normal(0, 5);
    sigma_tg ~ normal(0, 5);
    sigma_nb_phylo ~ normal(0, 5); //gamma(2, 2);
    sigma_tg_phylo ~ normal(0, 5); //gamma(2, 2);
    //prec_nb_phylo ~ gamma(2, 4);
    //prec_tg_phylo ~ gamma(2, 4);
    eta_sdep ~ normal(0, 1);
    eta_study ~ normal(0, 1);
    eta_obs ~ normal(0, 1);
    eta_nb ~ normal(0, 1);
    eta_tg ~ normal(0, 1);
    eta_nb_phylo ~ normal(0, 1);
    eta_tg_phylo ~ normal(0, 1);
}
generated quantities {
    vector[N] log_lik; // log likelihood
    //real <lower=0> p_g12; // probability of BOG1 < BOG2
    //real <lower=0> p_g13;
    //real <lower=0> p_g14;
    //real <lower=0> p_g23;
    //real <lower=0> p_g24;
    //real <lower=0> p_g34;
    real <lower=0> sigma2_tot; // total heterogeneity
    real <lower=0, upper=1> I2_sdep;
    real <lower=0, upper=1> I2_study; 
    real <lower=0, upper=1> I2_nb; 
    real <lower=0, upper=1> I2_tg; 
    real <lower=0, upper=1> I2_nb_phylo; 
    real <lower=0, upper=1> I2_tg_phylo;
    real <lower=0, upper=1> I2_obs; // proportion of heterogeneity explained by observation
    real <lower=0, upper=1> H2_nb; // phylogenetic signal for neighbor species
    real <lower=0, upper=1> H2_tg; // phylogenetic signal for target species
    real <lower=0, upper=1> I2_total;
    for (i in 1:N) log_lik[i] = normal_lpdf(d[i] | theta[i], sigma_d[i]);
    sigma2_tot = sigma_sdep^2 + sigma_study^2 + sigma_nb_phylo^2 + sigma_tg_phylo^2 + sigma_nb^2 + sigma_tg^2 + sigma_obs^2 + sigma2_d;
    I2_sdep = sigma_sdep^2 / sigma2_tot; 
    I2_study = sigma_study^2 / sigma2_tot;
    I2_nb = sigma_nb^2 / sigma2_tot;
    I2_tg = sigma_tg^2 / sigma2_tot;
    I2_nb_phylo = sigma_nb_phylo^2 / sigma2_tot;
    I2_tg_phylo = sigma_tg_phylo^2 / sigma2_tot;
    I2_obs = sigma_obs^2 / sigma2_tot;
    H2_nb = sigma_nb_phylo^2 / (sigma2_tot - sigma2_d);
    H2_tg = sigma_tg_phylo^2 / (sigma2_tot - sigma2_d);
    I2_total = I2_sdep + I2_study + I2_nb + I2_tg + I2_nb_phylo + I2_tg_phylo + I2_obs;
    //p_g12 = (beta[1] >= beta[2]);
    //p_g13 = (beta[1] >= beta[3]);
    //p_g14 = (beta[1] >= beta[4]);
    //p_g23 = (beta[2] >= beta[3]);
    //p_g24 = (beta[2] >= beta[4]);
    //p_g34 = (beta[3] >= beta[4]);
}
"

# compile the model
writeLines(modelString, con = "TEMPmodelStan.txt")
stan_code <- readChar("TEMPmodelStan.txt", file.info("TEMPmodelStan.txt")$size)

# data for Stan
dataList <- list(N = nrow(data), 
                 J = length(unique(data$pb.id)), 
                 jj = as.numeric(as.factor(data$pb.id)),
                 K = length(unique(data$s.dep)), 
                 kk = as.numeric(as.factor(data$s.dep)),
                 M_nb = length(unique(data$neighbor.name)),
                 M_tg = length(unique(data$target.name)),
                 mm_nb = mm_nb,
                 mm_tg = mm_tg,
                 VCV_nb = VCV_nb,
                 VCV_tg = VCV_tg,
                 I_nb = rep(1, ncol(VCV_nb)),
                 I_tg = rep(1, ncol(VCV_tg)),
                 d = data$yi, 
                 sigma_d = sqrt(data$vi),
                 sigma2_d = sigma2_d,
                 X = X,
                 p = ncol(X))


# sample from posterior
mrem.stan <- stan(model_code = stan_code, data = dataList, 
                  chains = chains, iter = niter, warmup = warmup, 
                  thin = thin, seed = 500, control = list(stepsize = .5, adapt_delta = .999))
                  
return(mrem.stan)

}