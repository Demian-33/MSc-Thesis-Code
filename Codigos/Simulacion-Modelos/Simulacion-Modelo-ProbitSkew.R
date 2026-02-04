
# Simulación para el modelo continuo ====
rm(list=ls())
getwd()

# (0) Cargar librerias, funciones y modelo ====
source('funciones.R')
library(cmdstanr)

out <- posterior::default_summary_measures()
myquant <- function(a){
  b <- quantile(a, c(0.025, 0.975))
  return(b)
}
out[5] <- 'myquant'

new_dir <- 'probit-sesgado-sim'
if(!dir.exists(new_dir)){
  dir.create(new_dir)
}


# En general, parece que VB no puede manejar varios mu_i para cada region
# con un solo mu, obtiene resultados buenos
# Parece que no es defecto del modelo, ya que HMC si puede manejar varios mu_i

# (1) Definir interceptos y porcentaje de muestreo ====
interceptos <- c('varios','unico','no') # no estoy seguro de incluir uno o muchos
porcentaje <- c(5, 25) # creo que 95 si se tarda mucho con HMC
int <- 'unico'
per <- 25

for(int in interceptos){
  
  print(paste0('La fecha actual es: ', Sys.time()))
  
  cat(int, 'interceptos.\n')
  
  if(int == 'varios'){
    mod_sn <- cmdstan_model(stan_file='ProbitSN-rho01-SSVS-genq.stan', compile=T)
  } else{
  if(int == 'unico'){
    mod_sn <- cmdstan_model(stan_file='ProbitSN-rho01-SSVS-genq-oneInt.stan', compile=T) # puede tener futuro
  } else{
  if(int == 'no'){
    mod_sn <- cmdstan_model(stan_file='ProbitSN-rho01-SSVS-genq-design.stan', compile=T) # puede no incluir mu[i]
    }
  }
}
  
for(per in porcentaje){
  
  cat('Usando el porcentaje', per, '%\n')
  
  # (2) Obtener datos y Simular y_{ij} ====
  set.seed(1)
  ifelse(int %in% c('varios', 'unico'), bool <- F, bool <- T)
  datos  <- get_sim(prop=per/100, type='Probit', center = F, noInt=bool, restrict_rho = T)
  datos$y_all <- get_levels(0, datos$y_all)
  datos$y_obs <- get_levels(0, datos$y_obs)
  datos$y_mis <- get_levels(0, datos$y_mis)
  ifelse(int=='varios', Actual <- c(datos$rho, datos$mu, datos$b),
  ifelse(int=='unico',  Actual <- c(datos$rho, paste0(datos$mu, collapse='|'), datos$b),
  Actual <- c(datos$rho, datos$b)))
  
  # datos$X_all <- cbind(generate_design(datos$N), datos$X_all)
  # datos$X_obs <- cbind(generate_design(datos$N)[-datos$idx_mis,], datos$X_obs)
  # datos$X_mis <- cbind(generate_design(datos$N)[datos$idx_mis,], datos$X_mis)
  # datos$p <- datos$p + datos$K

  # (3) Preparar datos y valores iniciales ====
  dat <- list(
    K = datos$K,
    N_obs = datos$N_obs,
    N_mis = datos$N_mis,
    N0 = length(which(datos$y_obs==0)),
    N1 = length(which(datos$y_obs==1)),
    n0 = which(datos$y_obs==0),
    n1 = which(datos$y_obs==1),
    p = datos$p,
    y_obs = datos$y_obs,
    X_obs = datos$X_obs,
    X_mis = datos$X_mis,
    X_all = datos$X_all,
    group = datos$group,
    group_mis = datos$groupmis,
    group_all = datos$group_all,
    tau=1/30,
    c=sqrt(10 * 30**2),
    epsilon=0.1
  )
  
  tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs)
  # r0 <- glm(y_obs ~ . + 0  , data=tmp, family=binomial(link='probit'))
  r0 <- glm(y_obs ~ . , data=tmp, family=binomial(link='probit'))
  mu0 <- coef(r0)[1]
  
  set.seed(1)
  ifelse(int=='varios', mu0 <- rep(0, datos$K), mu0 <- 0)
  init_list <- list(
    mu=mu0,
    b=coef(r0)[-1],
    rho=rep(0.5, times=datos$K),
    # z0=-abs(predict(r0, type='link'))[dat$n0],
    # z1=abs(predict(r0, type='link'))[dat$n1],
    # w_std=abs(rnorm(datos$K, mean=0.5)),
    pr=rep(0.5, times=datos$p)
  )
  
  # (4) Aproximación VB ----
  cat('Usando VB', '\n\n')
  print(paste0('La fecha actual es: ', Sys.time()))
  t0 <- proc.time()
  muestra_VB <- mod_sn$variational(
    data = dat,
    init = list(init_list),
    grad_samples=1,
    elbo_samples=100,
    eval_elbo=100,
    adapt_iter=1000,
    adapt_engaged=T,
    eta=0.05,
    tol_rel_obj=1e-5,
    algorith='meanfield',
    draws=1000,
    iter=75000
  )
  invisible(ifelse(!bool, out_var <- c('rho', 'mu', 'b'), out_var <- c('rho', 'b')))
  muestra_vb <- muestra_VB$summary(variables = out_var, out)
  tproc_vb <- proc.time()-t0 # Tiempo de proceso
  cat("\n","Tiempo de ejecución:",tproc_vb[3],"(elapsed)")
  print(muestra_vb, n=18)
  
  ## (4.0) Pathfinder ====
  # t0 <- proc.time()
  # muestra_Path <- mod_sn$pathfinder(
  #   seed=123,
  #   max_lbfgs_iters=2000,
  #   num_elbo_draws = 50,
  #   data = dat,
  #   num_paths = 8,
  #   init = list(init_list, init_list, init_list, init_list,
  #               init_list, init_list, init_list, init_list)
  # )
  # muestra_path <- muestra_Path$summary() # ajustar al 2.5% - 97.5%
  # tproc_path <- proc.time()-t0 # Tiempo de proceso
  # cat("\n","Tiempo de ejecución:",tproc_path[3],"(elapsed)")
  # print(muestra_path, n=18)
  
  ## (4.1) Fijar rutas ====
  path0 <- paste0('./probit-sesgado-sim/muestra_vb_', per, '_', int, '.csv')
  path1 <- paste0('./probit-sesgado-sim/predict_vb_', per, '_', int, '.csv')
  path2 <- paste0('./probit-sesgado-sim/SSVS_vb_', per, '_', int, '.csv')
  path3 <- paste0('./probit-sesgado-sim/mar_freq_vb_', per, '_', int, '.csv')
  ## (4.2) Guardar estimaciones ====
  df_out <- rbind(cbind(Actual, muestra_vb[, c(2, 4, 6, 7)]), tproc_vb[3])
  write.csv(df_out, path0)
  ## (4.3) Guardar predicciones ====
  yhat_SN <- Step(apply(muestra_VB$draws('z_all'), 2, mean)[datos$idx_mis])
  write.csv(cbind(yhat_SN, datos$y_mis), path1)
  ## (4.4) Guardar probabilidad posterior ====
  write.csv(posterior_prob(muestra_VB), path2)
  write.csv(marginal_pi(muestra_VB), path3)
  
  # (5) Usando HMC ----
  cat('Usando HMC', '\n\n')
  print(paste0('La fecha actual es: ', Sys.time()))
  t0 <- proc.time()
  muestra_HMC <- mod_sn$sample(
    data = dat,
    init = list(init_list),
    chains=1,
    iter_sampling = 1000,
    iter_warmup = 5000
  )
  muestra_hmc <- muestra_HMC$summary(variables = out_var, out)
  tproc_hmc <- proc.time()-t0 # Tiempo de proceso
  cat("\n","Tiempo de ejecución:",tproc_hmc[3],"(elapsed)")
  print(muestra_hmc, n=18)
  
  ## (5.1) Fijar rutas ====
  path0 <- paste0('./probit-sesgado-sim/muestra_hmc_', per, '_', int,  '.csv')
  path1 <- paste0('./probit-sesgado-sim/predict_hmc_', per, '_', int,  '.csv')
  path2 <- paste0('./probit-sesgado-sim/SSVS_hmc_', per, '_', int,  '.csv')
  path3 <- paste0('./probit-sesgado-sim/mar_freq_hmc_', per, '_', int,  '.csv')
  ## (5.2) Guardar estimaciones ====
  df_out <- rbind(cbind(Actual, muestra_hmc[, c(2, 4, 6, 7)]), tproc_hmc[3])
  write.csv(df_out, path0)
  ## (5.3) Guardar predicciones ====
  yhat_SN <- Step(apply(drop(muestra_HMC$draws('z_all')), 2, mean))[datos$idx_mis]
  write.csv(cbind(yhat_SN, datos$y_mis), path1)
  ## (5.4) Guardar probabilidad posterior ====
  write.csv(posterior_prob(muestra_HMC), path2)
  write.csv(marginal_pi(muestra_HMC), path3)
  }
}
# 
# print(muestra_vb, n=18)
# print(muestra_hmc, n=18)
# 
# b0 <- Step(apply(muestra_VB$draws('z_all'), 2, mean))
# metrics_cat(b0[datos$idx_mis], datos$y_mis)
# table(b0[datos$idx_mis], datos$y_mis)
# 
# b1 <- Step(apply(drop(muestra_HMC$draws('z_all')[, 1, ]), 2, mean))
# metrics_cat(b1[datos$idx_mis], datos$y_mis)
# table(b1[datos$idx_mis], datos$y_mis)
# 
# b2 <- Step(predict(r0, newdata = data.frame(datos$X_mis), type='response'), delta=0.5)
# table(b2, datos$y_mis)
# metrics_cat(b2, datos$y_mis)
# 
# b3 <- Step(apply(muestra_Path$draws('z_all'), 2, mean))
# metrics_cat(b3[datos$idx_mis], datos$y_mis)
