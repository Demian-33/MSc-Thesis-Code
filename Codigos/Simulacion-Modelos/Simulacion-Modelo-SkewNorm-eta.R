
# Simulación para el modelo continuo ====
rm(list=ls())

# (0) Cargar librerias y funciones ====
source('funciones.R')
library(cmdstanr)

new_dir <- 'log-normal-sesgado-sim-eta'
if(!dir.exists(new_dir)){
  dir.create(new_dir)
}

out <- posterior::default_summary_measures()
myquant <- function(a){
  b <- quantile(a, c(0.025, 0.975))
  return(b)
}
out[5] <- 'myquant'


# (1) Definir eta, interceptos y porcentaje de muestreo ====
interceptos <- c('varios') # ¿como se justifica no incluir mu's?
porcentaje <- c(25)
eta <- seq(0.0, 1.5, by=0.05)
eta[1] <- 0.01
length(eta)

for(step_size in eta){
  
  cat('Usando eta =', step_size, '\n')

for(int in interceptos){
  
  print(paste0('La fecha actual es: ', Sys.time()))
  
  cat(int, 'interceptos.\n')
  
  if(int == 'varios'){
    mod_sn <- cmdstan_model(stan_file='LogSN-rho01-SSVS-genq.stan', compile=T)
    # mod_sn <- cmdstan_model(stan_file='LogSN-rho11-SSVS-genq.stan', compile=T)
  } else{
    if(int == 'unico'){
      mod_sn <- cmdstan_model(stan_file='LogSN-rho01-SSVS-genq-oneInt.stan', compile=T) # puede tener futuro
    } else{
      if(int == 'no'){
        mod_sn <- cmdstan_model(stan_file='LogSN-rho01-SSVS-genq-design.stan', compile=T) # puede no incluir mu[i]
      }
    }
  }
  
  for(per in porcentaje){
    
    cat('Usando el porcentaje', per, '%\n')
    
    # (2) Obtener datos y Simular y_{ij} ====
    set.seed(1)
    ifelse(int %in% c('varios', 'unico'), bool <- F, bool <- T)
    datos  <- get_sim(prop=per/100, noInt = bool, restrict_rho = T)
    ifelse(int=='varios', Actual <- c(datos$rho, datos$sigma, datos$mu, datos$b),
    ifelse(int=='unico',  Actual <- c(datos$rho, datos$sigma, paste0(datos$mu, collapse='-'), datos$b),
    Actual <- c(datos$rho, datos$sigma, datos$b)))
    
    # (3) Preparar datos y valores iniciales ====
    dat <- list(
      K = datos$K,
      N_obs = datos$N_obs,
      N_mis = datos$N_mis,
      p = datos$p,
      y_obs = datos$y_obs,
      X_obs = datos$X_obs,
      X_mis = datos$X_mis,
      X_all = datos$X_all,
      group = datos$group,
      groupmis = datos$groupmis,
      group_all = datos$group_all,
      tau=1/30,
      c=sqrt(10 * 30**2),
      epsilon=0.1
    )
    
    tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs)
    init_list <- get_inits(tmp, restrict_rho=T)
    init_list$pr <- rep(0.5, times=length(datos$b))
    ifelse(int=='varios', function(){}, init_list$mu <- mean(init_list$mu))
    
    # (4) Aproximación VB ----
    cat('Usando VB', '\n\n')
    print(paste0('La fecha actual es: ', Sys.time()))
    t0 <- proc.time()
    muestra_VB <- mod_sn$variational(
      data = dat,
      init = list(init_list),
      grad_samples=1,
      elbo_samples=100,
      eval_elbo=1000,
      adapt_iter=1000,
      adapt_engaged=FALSE,
      eta=step_size,
      tol_rel_obj=1e-5,
      algorith='meanfield',
      draws=1000,
      iter=75000
    )
    invisible(ifelse(!bool, out_var <- c('rho', 'sigma', 'mu', 'b'), out_var <- c('rho', 'sigma', 'b')))
    muestra_vb <- muestra_VB$summary(variables = out_var, out)
    tproc_vb <- proc.time()-t0 # Tiempo de proceso
    cat("\n","Tiempo de ejecución:",tproc_vb[3],"(elapsed)")
    print(muestra_vb, n=18)
    # r1 <- muestra_VB$draws('rho')[, 1]
    # r4 <- muestra_VB$draws('rho')[, 4]
    # x11(); hist(r4)
    
    ## (4.1) Fijar rutas ====
    path0 <- paste0('./log-normal-sesgado-sim-eta/muestra_vb_', per, '_', int, '_eta', step_size, '.csv')
    path1 <- paste0('./log-normal-sesgado-sim-eta/predict_vb_', per, '_', int, '_eta', step_size, '.csv')
    path2 <- paste0('./log-normal-sesgado-sim-eta/SSVS_vb_', per, '_', int, '_eta', step_size, '.csv')
    path3 <- paste0('./log-normal-sesgado-sim-eta/mar_freq_vb_', per, '_', int, '_eta', step_size, '.csv')
    ## (4.2) Guardar estimaciones ====
    df_out <- rbind(cbind(Actual, muestra_vb[, c(2, 4, 6, 7)]), tproc_vb[3])
    write.csv(df_out, path0)
    ## (4.3) Guardar predicciones ====
    yhat_SN <- apply(muestra_VB$draws('y_all'), 2, mean)[datos$idx_mis]
    write.csv(cbind(yhat_SN, datos$y_mis), path1)
    ## (4.4) Guardar probabilidad posterior ====
    write.csv(posterior_prob(muestra_VB), path2)
    write.csv(marginal_pi(muestra_VB), path3)
    }
  }
}

# a <- apply(muestra_VB$draws('y_all'), 2, mean)
# metrics(a[datos$idx_mis], datos$y_mis)
# x11();plot(a[datos$idx_mis], datos$y_mis);grid()
# 
# b0 <- apply(drop(muestra_HMC$draws('y_all')[, 1, ]), 2, mean)
# b1 <- apply(drop(muestra_HMC$draws('y_all')[, 2, ]), 2, mean)
# 
# metrics(b0[datos$idx_mis], datos$y_mis)
# metrics(b1[datos$idx_mis], datos$y_mis)
# points(b0[datos$idx_mis], datos$y_mis, col=3)
# points(b1[datos$idx_mis], datos$y_mis, col=4)

