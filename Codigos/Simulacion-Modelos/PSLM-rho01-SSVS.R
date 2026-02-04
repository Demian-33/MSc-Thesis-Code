##______________________________________________________________________________
# Código STAN para la estimación áreas pequeñas bajo el modelo probit asimétrico
# v1 - 'modelo latente adaptado de Stan User's Guide (SUG)'

# Respuesta:
#   y \in {0, 1, 2, ..., L}
# Variables latentes:
#   Z: directas en el modelo
#   W: directas en el modelo
# Estados:
#   Único estado
# Parámetro de forma:
#   Correlación en (0, 1)
# Parametrización:
#   Centrada:
#   E(Z)   = \mu + \sigma\rho\sqrt(2/\pi)
#   Var(Z) = \sigma^2 (1 - 2\rho^2/\pi)^2
# Vale la pena notar que y_{ij} depende del valor de z_{ij} (var. latente),
#   además z_{ij} se genera a partir de truncamiento oculto, es decir que
#   z_{ij} se genera a través de [U_{i}, W_{i}]^T
#   es decir, también existe dependencia en z_{ij} y en w_{i}
##______________________________________________________________________________

codigo_stan <- '
data {
  int<lower=1> K;                         // Number of groups (municipalities)
  int<lower=0> N_obs;                     // Number of observed data
  int<lower=0> N_mis;                     // Number of missing data
  int<lower=1> p;                         // Number of explanatory variables
  array[N_obs] int<lower=0, upper=2> y;   // Observed dependent variable for all groups
  matrix[N_obs, p] X;                     // Matrix of observed explanatory variables for all groups
  matrix[N_mis, p] X_mis;                 // Matrix of observed explanatory variables for all groups with missing y
  int<lower=1, upper=K> group[N_obs];     // Group assignment for each observation
  int<lower=1, upper=K> groupmis[N_mis];  // Group assignment for each missing observation
  // Stan Users Guide Restrictions
  int<lower=0> N0;                        // Number of 0 observations
  int<lower=0> N1;                        // Number of 1 observations
  int<lower=0> N2;                        // Number of 2 observations
  array[N0] int<lower=1, upper=N_obs> n0; // Position of 0 observations
  array[N1] int<lower=1, upper=N_obs> n1; // Position of 1 observations
  array[N2] int<lower=1, upper=N_obs> n2; // Position of 2 observations
  // SSVS
  real<lower=0> tau2;                     // Estimation of zero
  real<lower=1> c2;                       // Estimation of a wide range of values
  real epsilon;  // Threshold for significant coefficient
}

transformed data {
  real<lower=0> sigma=1.0;                  // standard error parameter for latent variable z
  real delta0=0.0;                          // baseline threshold
}

parameters {
  vector<lower=0, upper=1>[K] rho;       // Correlation coefficient for each group
  vector[p] b;                           // Regression coefficients for explanatory variables
  // Restrictions for multilevel probit
  real<lower=0> delta1;                       // first threshold to be estimed
  vector<lower=0, upper=1>[p] pr;        // Inclusion probability of each beta_i
  vector<upper=delta0>[N0] z0;                //          z_{ij} < delta0, N0 < N_obs
  vector<lower=delta0, upper=delta1>[N1] z1;  // delta0 < z_{ij} < delta1, N1 < N_obs
  vector<lower=delta1>[N2] z2;                //          z_{ij} > delta1, N2 < N_obs
  vector[N_mis] zmis;                    // Latent variable for missing y_{ij}
  vector<lower=0>[K] w;                  // hidden truncated variable w_{i}
}

transformed parameters {
 // why not use log lambda? or log rho?
 //vector<lower=0,upper=1>[K] rho=lambda./sqrt(1+lambda.^2); // Correlation coef. for each group
 //vector<lower=-1,upper=1>[K] rho=lambda./sqrt(1+lambda.^2); // Correlation coef. for each group
 
  vector[N_obs] z;
  for (n in 1:N0) {
    z[n0[n]] = z0[n];
  }
  for (n in 1:N1) {
    z[n1[n]] = z1[n];
  }
  for (n in 1:N2) {
    z[n2[n]] = z2[n];
  }
}

model {
  // Here we consider centered parameterization
  
  // Reference prior for lambda
  // target += sum(0.5*log(1+2*lambda.^2) - 2*log(1+lambda.^2));
  
  // Reference prior for rho
  target += sum(0.5*log1p(rho.^2) - log1m(rho.^2));
  
  // Wide prior for the threshold
  delta1 ~ normal(0, 100);
  
  // Priors for regression coefficients b
  // b ~ normal(0, 100);  // Wide priors for regression coefficients
  
  // Prior for SSVS: seleccion of b and his inclusion probability pr.
  pr ~ beta(0.5, 0.5);
  for(i in 1:p){
    //pr[i] ~ beta(0.5, 0.5);
    target += log_mix(1-pr[i],
    normal_lpdf(b[i] | 0, tau2),
    normal_lpdf(b[i] | 0, tau2*c2));
  }
  
  // Prior for w (hidden truncated variable) for each group
  w ~ normal(0, sigma ./ sqrt(1 - 2/pi() * square(rho)) );  // Truncated normal prior for latent variables
  
  // Likelihood for observed y given z and b for each group
  // `sig` means the square root of variance parameter in the centered skew normal model
  
  vector[N_obs] sig = sigma * sqrt(1 - square(rho[group]))./ sqrt(1 - 2/pi() * square(rho[group]));
  vector[N_obs] intercept = (rho[group] .* w[group]) - sigma * rho[group] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[group]));
  vector[N_obs] eta = X*b + intercept;
  // vector[N_obs] eta = (X*b + intercept)./sig;
  
  vector[N_mis] sigmis = sigma * sqrt(1 - square(rho[groupmis]))./ sqrt(1 - 2/pi() * square(rho[groupmis]));
  vector[N_mis] interceptmis = (rho[groupmis] .* w[groupmis]) - sigma * rho[groupmis] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[groupmis]));
  vector[N_mis] etamis = X_mis*b + interceptmis;
  // vector[N_mis] etamis = (X_mis*b + interceptmis)./sigmis;
  
  // Latent variable for observed y_{ij}
  z ~ normal(eta, sig);
  
  // Latent variable for missing y_{ij}
  zmis ~ normal(etamis, sigmis);
}

generated quantities {
    vector<lower=0, upper=1>[p] m_ind;
    for (j in 1:p) {
        if (fabs(b[j]) > epsilon) {
            m_ind[j] = 1;
        } else {
            m_ind[j] = 0;
        }
    }
}
'
#_______________________________________________________________________________
