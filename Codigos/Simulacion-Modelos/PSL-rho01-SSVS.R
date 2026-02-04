##______________________________________________________________________________
# Código STAN para la estimación en áreas pequeñas, modelo probit asimétrico
# Respuesta:
#   y \in {0, 1}
# Variables latentes:
#   Z: directas en el modelo
#   W: directas en el modelo
# Estados:
#   Únicamente áreas pequeñas en un estado
# Parámetro de forma:
#   Correlación en (0, 1)
# Parametrización:
#   Centrada:
#   E(Z)   = \mu + \sigma\rho\sqrt(2/\pi)
#   Var(Z) = \sigma^2 (1 - 2\rho^2/\pi)^2
#   Y = mu + sigma^2 (Y-E(Z))/SD(Z)
# Vale la pena notar que y_{ij} depende del valor de z_{ij} (var. latente),
#   además, z_{ij} se genera a partir de truncamiento oculto, es decir que
#   z_{ij} se genera a través de [U_{i}, W_{i}]^T
#   por lo tanto, también existe dependencia en z_{ij} y en w_{i}
##______________________________________________________________________________

codigo_stan <- '
data {
  int<lower=1> K;                         // Number of groups (municipios)
  int<lower=0> N_obs;                     // Number of observed data
  int<lower=0> N_mis;                     // Number of missing data
  int<lower=1> p;                         // Number of explanatory variables
  array[N_obs] int<lower=0, upper=1> y_obs;   // Observed dependent variable for all groups
  matrix[N_obs, p] X_obs;                     // Matrix of observed explanatory variables for all groups
  matrix[N_mis, p] X_mis;                 // Matrix of observed explanatory variables for all groups with missing y
  int<lower=1, upper=K> group[N_obs];     // Group assignment for each observation
  int<lower=1, upper=K> groupmis[N_mis];  // Group assignment for each missing observation
  // Restrictions
  int<lower=0> N0;                        // Number of 0 observations
  int<lower=0> N1;                        // Number of 1 observations
  array[N0] int<lower=1, upper=N_obs> n0; // Position of 0 observations
  array[N1] int<lower=1, upper=N_obs> n1; // Position of 1 observations
  // SSVS
  real<lower=0> tau2;                     // Estimation of zero
  real<lower=1> c2;                       // Estimation of a wide range of values
  real epsilon;                           // Threshold for significant coefficient
}

transformed data {
  real<lower=0> sigma=1.0;                // standard error parameter for latent variable z
}

parameters {
  vector<lower=0, upper=1>[K] rho;       // Correlation coefficient for each group
  vector<lower=0>[K] w;                  // hidden truncated variable w_{i}
  vector[K] mu;
  vector[p] b;                           // Regression coefficients for explanatory variables
  vector<upper=0>[N0] z0;                // No. of z_{ij} < 0, N0 < N_obs, restricted
  vector<lower=0>[N1] z1;                // No. of z_{ij} > 0, N1 < N_obs, restricted
  vector[N_mis] z_mis;                    // Latent variable for missing y_{ij}, no restriccions
  vector<lower=0, upper=1>[p] pr;        // Inclusion probability of each beta_i
}

transformed parameters {
  // Restricciones en el modelo probit
  vector[N_obs] z;  
  for (n in 1:N1) {
    z[n1[n]] = z1[n];
  }
  for (n in 1:N0) {
    z[n0[n]] = z0[n];
  }
}

model {
  // a priori de referencia para rho
  target += sum(0.5*log1p(rho.^2) - log1m(rho.^2));
  // a priori plana para mu
  mu ~ normal(0, 10);
  // a priori para pr (probabilidad de inclusión)
  pr ~ beta(0.5, 0.5);
  // a priori SSVS para beta (selección de variables)
  for(i in 1:p){
    target += log_mix(pr[i],
    normal_lpdf(b[i] | 0, tau2*c2),
    normal_lpdf(b[i] | 0, tau2));
  }
  
  // Para centrar parametros
  vector[K]       sd_sn        = sqrt(1 - 2*square(rho)/pi());
  vector[N_obs] sd_group_sn    = sqrt(1 - 2*square(rho[group])/pi());
  vector[N_mis] sd_groupmis_sn = sqrt(1 - 2*square(rho[groupmis])/pi());
  
  // Prior for w (hidden truncated variable) for each group
  w ~ normal(0, sigma ./ sd_sn);  // Truncated normal prior for latent variables
  
  // Likelihood for observed y given z and b for each group
  // `sig` means the square root of variance parameter in the centered skew normal model
  
  vector[N_obs] sig = sigma * sqrt(1 - square(rho[group])) ./ sd_group_sn;
  vector[N_obs] intercept = (rho[group] .* w[group]) - (sigma * rho[group] .* sqrt(2/pi()) ./ sd_group_sn);
  vector[N_obs] eta = mu[group] + X_obs*b + intercept;

  vector[N_mis] sigmis = sigma * sqrt(1 - square(rho[groupmis]))./ sd_groupmis_sn;
  vector[N_mis] interceptmis = (rho[groupmis] .* w[groupmis]) - (sigma * rho[groupmis] .* sqrt(2/pi()) ./ sd_groupmis_sn);
  vector[N_mis] etamis = mu[groupmis] + X_mis*b + interceptmis;
  
  // Latent variable for observed y_{ij}
  z ~ normal(eta, sig);
  
  // Latent variable for missing y_{ij}
  z_mis ~ normal(etamis, sigmis);
}

generated quantities {
  // SSVS
  vector<lower=0, upper=1>[p] m_ind;
  for (j in 1:p) {
      if (fabs(b[j]) > epsilon) {
          m_ind[j] = 1;
      } else {
          m_ind[j] = 0;
      }
  }
  // All fitted values
  //vector[N_obs+N_mis] y_all;
  //for(i in 1:(N_obs+N_mis)){
  //  j = region_all[i];
  //  y_all[i] = 
  //}
}
'
#_______________________________________________________________________________
