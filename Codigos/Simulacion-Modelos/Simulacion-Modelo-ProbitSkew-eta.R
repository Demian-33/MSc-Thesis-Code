
# Simulación para el modelo continuo ====
rm(list=ls())
getwd()

# (0) Cargar librerias, funciones y modelo ====
source('funciones.R')
library(cmdstanr)
library(dplyr)

out <- posterior::default_summary_measures()
myquant <- function(a){
  b <- quantile(a, c(0.025, 0.975))
  return(b)
}
out[5] <- 'myquant'

new_dir <- 'probit-sesgado-sim-eta'
if(!dir.exists(new_dir)){
  dir.create(new_dir)
}


# En general, parece que VB no puede manejar varios mu_i para cada region
# con un solo mu, obtiene resultados buenos
# Parece que no es defecto del modelo, ya que HMC si puede manejar varios mu_i

# (1) Definir eta, interceptos y porcentaje de muestreo ====
interceptos <- c('varios', 'unico', 'no') # no estoy seguro de incluir uno o muchos
porcentaje <- c(25) # creo que 95 si se tarda mucho con HMC
eta <- seq(0.01, 1.5, by=0.1) # explorar menores a 0.01)
eta[1] <- 0.01
eta <- seq(0.001, 0.01, by=0.001)
eta <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0)
gradiente <- c(1, 5, 10)

for(grad in gradiente){

for(step_size in eta){
  
for(int in interceptos){
  
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
    print(paste0('La fecha actual es: ', Sys.time()))
    
    # (2) Obtener datos y Simular y_{ij} ====
    set.seed(1)
    ifelse(int %in% c('varios', 'unico'), bool <- F, bool <- T)
    datos  <- get_sim(prop=per/100, type='Probit', center = F, noInt=bool)
    datos$y_all <- get_levels(0, datos$y_all)
    datos$y_obs <- get_levels(0, datos$y_obs)
    datos$y_mis <- get_levels(0, datos$y_mis)
    ifelse(int=='varios', Actual <- c(datos$rho, datos$mu, datos$b),
    ifelse(int=='unico',  Actual <- c(datos$rho, paste0(datos$mu, collapse='-'), datos$b),
    Actual <- c(datos$rho, datos$b)))

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
      grad_samples=grad,
      elbo_samples=100,
      eval_elbo=1000,
      adapt_iter=1000,
      adapt_engaged=F,
      eta=step_size,
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
    
    
    ## (4.1) Fijar rutas ====
    path0 <- paste0('./probit-sesgado-sim-eta/muestra_vb_', '_', per, '_', int, '_eta', step_size, '_grad', grad, '.csv')
    path1 <- paste0('./probit-sesgado-sim-eta/predict_vb_', '_', per, '_', int, '_eta', step_size, '_grad', grad,   '.csv')
    path2 <- paste0('./probit-sesgado-sim-eta/SSVS_vb_', '_', per, '_', int, '_eta', step_size, '_grad', grad,   '.csv')
    path3 <- paste0('./probit-sesgado-sim-eta/mar_freq_vb_', '_', per, '_', int, '_eta', step_size, '_grad', grad,   '.csv')
    ## (4.2) Guardar estimaciones ====
    df_out <- rbind(cbind(Actual, muestra_vb[, c(2, 4, 6, 7)]), tproc_vb[3])
    write.csv(df_out, path0)
    ## (4.3) Guardar predicciones ====
    yhat_SN <- get_levels(0, apply(muestra_VB$draws('z_all'), 2, mean)[datos$idx_mis])
    write.csv(cbind(yhat_SN, datos$y_mis), path1)
    ## (4.4) Guardar probabilidad posterior ====
    write.csv(posterior_prob(muestra_VB), path2)
    write.csv(marginal_pi(muestra_VB), path3)
      }
    }
  }
}
id <- 1:3
a <- c('Level 1', 'Level 2', 'Level 3')
b <- c(0.0, 1.0, 2.5)
d <- c(10, 20, 80)
r <- rnorm(n=3)
temp <- data.frame(id, a, b, d, r)

library(reshape2)
temp <- melt(temp, id.vars = id)

g0 <- ggplot(data=temp, aes(x=a, y=b, fill=r, linewidth = e)) + 
  geom_tile( color='pink') + 
  facet_grid(~d)
x11(); print(g0)



