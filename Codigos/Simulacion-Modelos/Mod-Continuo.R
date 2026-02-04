library(cmdstanr)

stan_code_pool <- '
data {
  int<lower=1> K;  // regions
  int<lower=0> N_obs; // no. of observed y_{ij}
  int<lower=0> N_mis; // no. of missing  y_{ij}
  array[N_obs] int<lower=1, upper=K> group; // region of observed
  array[N_mis] int<lower=1, upper=K> groupmis; // region of missing
  int<lower=1> p;  // predictors
  vector[N_obs] y_obs;
  matrix[N_obs, p] X_obs;
  matrix[N_mis, p] X_mis;
}
parameters {
  //vector<lower=0>[K] lambda;
  vector<lower=0>[K] lambda_tilde;
  //real mu_lambda;
  //real<lower=0> sigma_lambda;
  real<lower=0> sigma;           // sd for all region
  vector[p] beta;                // Regression coefficients
  //vector[K] mu;                // intercept for each region
  //vector[K] mu_tilde;            // intercept for each region
  real mu_tilde;                 // intercept for each region
  real mu0;                    // global mean for mu
  real<lower=0> sigma_mu;        // sd for mu
  vector[N_mis] y_mis;           // missing y_{ij}
  vector<lower=0>[K] z_tilde;    // Latent variable for each group
}
transformed parameters{
  //vector<lower=0>[K] lambda = mu_lambda + sigma_lambda * lambda_tilde;
  vector<lower=0, upper=1>[K] rho = lambda ./ sqrt(1 + square(lambda));
  vector<lower=0>[K] z = z_tilde .* sigma ./ sqrt(1 - 2/pi() * square(rho));
  //vector[K] mu = mu_mu + sigma_mu * mu_tilde;
  real mu = mu_mu + sigma_mu * mu_tilde;
}
model {
  // target += sum( 0.5 * log1p(square(rho)) - log1m(square(rho)) );
  beta ~ normal(0, 10);
  mu ~ normal(0, 10);
  sigma ~ inv_gamma(0.001, 0.001);
  //lambda ~ normal(mu_lambda, sigma_lambda);
  lambda_tilde ~ normal(0, 1);
  mu_lambda ~ normal(0, 1);
  sigma_lambda ~ normal(0, 1);
  mu_tilde ~ normal(0, 1);
  mu_mu ~ normal(0, 5);
  sigma_mu ~ normal(0, 1);
  z ~ normal(0, 1); 
  vector[N_obs] sig = sigma * sqrt(1 - square(rho[group])) ./ sqrt(1 - 2/pi() * square(rho[group]));
  vector[N_obs] intercept = rho[group].* z[group] - sigma*rho[group] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[group]));
  target +=normal_id_glm_lupdf(y_obs | X_obs, mu + intercept, beta, sig);
  vector[N_mis] sigmis = sigma * sqrt(1 - square(rho[groupmis]))./ sqrt(1 - 2/pi() * square(rho[groupmis]));
  vector[N_mis] interceptmis = rho[groupmis] .* z[groupmis]- sigma*rho[groupmis] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[groupmis]));
  target +=normal_id_glm_lupdf(y_mis | X_mis, mu + interceptmis, beta, sigmis);
}'

stan_model_pool <- write_stan_file(stan_code_pool) # write .stan file
mod_pool <- cmdstan_model(stan_model_pool, compile=T)
