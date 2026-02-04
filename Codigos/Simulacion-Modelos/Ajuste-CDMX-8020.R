# =============== Analisis a nivel hogares del Ing. cor. total per capita =================

# este ajuste se realiza con solo la información de la ENIGH 2024
# se divide la población en un conjunto de 80% entrenamiento y 20% prueba

# manipulación y visualización
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(kableExtra)
# ajuste del modelo con Stan
library(cmdstanr)

# (0) Datos para el modelo ====
rm(list=ls())
gc()
load('./Xyregion-CDMX.Rdata')
source('./funciones-final_.R')

out_stan <- posterior::default_summary_measures()
myquant <- function(a){
  b <- quantile(a, c(0.025, 0.975))
  return(b)
}
out_stan[5] <- 'myquant'

muni <- c(
  'Azcapotzalco', 'Coyoacán', 'Cuajimalpa de Morelos', 'Gustavo A. Madero',
  'Iztacalco', 'Iztapalapa', 'La Magdalena Contreras', 'Milpa Alta',
  'Álvaro Obregón', 'Tláhuac', 'Tlalpan', 'Xochimilco', 
  'Benito Juárez','Cuauhtémoc', 'Miguel Hidalgo', 'Venustiano Carranza') %>% factor

LPI_urb <- 4564.97
LPI_rur <- 3296.92
LPEI_urb <- 2354.65
LPEI_rur <- 1800.55

## (0.1) Definir areas pequeñas ====
group <- as.numeric(region)
K <- length(unique(group))

## (0.2) Covariables categoricas (no ordenadas) a matriz diseño ====
## no se consideran en el analisis
# for(i in 1:ncol(Xcat)){
#   Xcat[, i] <- as.factor(as.integer(Xcat[, i]))
# }
# tail(tibble(Xcat)) # ya todas son factor
# Xcat <- model.matrix( ~ 0 + . , data=Xcat) # generar a matriz diseño

## (0.3) Remover combinaciones lineales y correlaciones ====
Xtmp <- cbind(Xcon, Xbin)
a <- caret::findLinearCombos(Xtmp)
Xtmp <- Xtmp[, -c(a$remove)]
b <- caret::findCorrelation(cov(Xtmp), verbose=T, exact=T, cutoff = 0.8)
Xtmp <- Xtmp[, -c(b)]

dim(Xtmp) # 110 covariables

## (0.4) Componentes principales ====
# Eigen <- eigen(cov(Xcon))
# cumsum(Eigen$values) / ncol(Xcon) # 26: 95%
# Xpc <- Xcon %*% Eigen$vectors[, 1:26]
# Xtmp <- cbind(Xpc, Xbin)

tproc_hmc <- c()

# (1) Modelo log-normal ====
file <- file.path('LogSN-rho01-SSVS-genq-oneInt.stan') # estimacion, un int
file <- file.path('LogSN-rho01-SSVS-genq.stan') # estimacion
mod_sn <- cmdstan_model(file, compile=T) 

## (1.1) Seleccionar datos ====
## los argumentos ya están por defecto, vea get_obs
set.seed(1)
datos <- get_obs(censo=F, prop=0.2)

dat <- list(
  K = K,
  N_obs = length(datos$y_obs),
  N_mis = length(datos$y_mis),
  p = ncol(datos$X_obs),
  y_obs = datos$y_obs,
  X_obs = datos$X_obs,
  X_mis = datos$X_mis,
  X_all = Xtmp[datos$master, ],
  group = datos$group[datos$observed],
  group_mis = datos$group[datos$missing],
  group_all = region[datos$master],
  # ssvs. Stan parametriza la normal con la desv. estándar
  tau=1/300,
  c=3000, # sqrt(10 * 300**2)
  epsilon=0.01
)
## (1.2) Valores iniciales ====
tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs)
r0 <- lm(y_obs ~ ., data=tmp)
r1 <- coef(sn::selm(y_obs ~ 1 , data=tmp), param.type = 'DP')
mu0 <- r1[1]
sigma0 <- r1[2]
rho0 <- lambda2rho(r1[3])
out_var <- c('rho', 'sigma', 'mu', 'b')

## datos para Stan
set.seed(2)
init_list <- list(
  b = get_coefs(r0, F),
  pr=rep(0.5, ncol(datos$X_mis)),
  sigma=sigma0,
  mu=rep(mu0, times=K),
  # mu=mu0,
  rho=rep(rho0, times=K),
  z_std=rnorm(K, 0.5, sd=0.1)
)

## (1.3) Correr el modelo con HMC ====
t0 <- proc.time()
muestra_HMC <- mod_sn$sample(
  output_dir = r'(C:\Users\Lesau\Documents\EntornosR\log-normal\HMC-8020)',
  data = dat,
  seed = 123,
  init = list(init_list),
  chains=1,
  iter_warmup = 10000,
  thin=2,
  iter_sampling = 2000
)
muestra_hmc <- muestra_HMC$summary(variables = out_var, out_stan)
tproc_hmc <- c(tproc_hmc, (proc.time()-t0)[3])
save(file=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-log-normal-sesgado-8020.RData)',
     muestra_HMC)

# (2) Modelo probit sesgado latente ====
# file <- file.path('ProbitSN-rho01-SSVS-genq-design.stan') # estimacion
file <- file.path('ProbitSN-rho01-SSVS-genq-oneInt.stan') # estimacion
file <- file.path('ProbitSN-rho01-SSVS-genq.stan') # estimacion
mod_sn <- cmdstan_model(file, compile=T)
## (2.1) Seleccionar datos ====
## 0: urbano, 1: rural
## ordinal: 0: ambas umbrales, 1: un umbral, 2: arriba de los umbrales
tmp <- tibble(ym, contexto=Xide$rururb, region=muni[region])
tmp <- tmp %>% mutate(
  ym_ord = case_when(
    is.na(ym) ~ NA,
    contexto==0 & ym < log(LPEI_urb) ~ 0,
    contexto==1 & ym < log(LPEI_rur) ~ 0,
    contexto==0 & ym > log(LPEI_urb) & ym < log(LPI_urb)  ~ 1,
    contexto==1 & ym > log(LPEI_rur) & ym < log(LPI_rur)  ~ 1,
    TRUE  ~ 2),
  ym_bin = case_when(
    ym_ord <= 1 ~ 1, # si carencias
    is.na(ym_ord) ~ NA,
    TRUE ~ 0) # no carencias
)
table(tmp$ym_bin)
table(tmp$ym_ord)

set.seed(1)
datos <- get_obs(censo=F, y=tmp$ym_bin, prop=0.2)
# no invertir 0's y 1's
# datos$y_obs <- ifelse(datos$y_obs==1, 0, 1)
dat <- list(
  K = K,
  N_obs = length(datos$y_obs),
  N_mis = length(datos$y_mis),
  p = ncol(datos$X_obs),
  # y_obs = datos$y_obs,
  y_obs = datos$y_obs,
  X_obs = datos$X_obs,
  X_mis = datos$X_mis,
  X_all = Xtmp[datos$master, ],
  group = datos$group[datos$observed],
  group_mis = datos$group[datos$missing],
  group_all = region[datos$master],
  # probits
  N0 = length(which(datos$y_obs==0)),
  N1 = length(which(datos$y_obs==1)),
  N2 = length(which(datos$y_obs==2)),
  n0 = which(datos$y_obs==0),
  n1 = which(datos$y_obs==1),
  n2 = which(datos$y_obs==2),
  # ssvs. Stan parametriza la normal con la desv. estándar
  tau=1/300,
  c=3000, #sqrt(10 * 300**2)
  epsilon=0.1
)

## (2.2) Valores iniciales ====
tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs)
r0 <- glm(y_obs ~ ., data=tmp, family=binomial(link='probit'))
mu0 <- 0
rho0 <- 0.5
out_var <- c('rho', 'mu', 'b')
set.seed(2)
# z <- get_z(datos$y_obs) # ?
init_list <- list(
  b = get_coefs(r0, F),
  pr=rep(0.5, ncol(datos$X_mis)),
  mu=rep(mu0, times=K),
  # mu=mu0, # un intercepto
  # z0 = z[z<0],
  # z1 = z[z>0],
  rho=rep(rho0, times=K),
  w_tilde=rnorm(K, 0.5, sd=0.1)
)

## (2.3) Correr el modelo con HMC ====
t0 <- proc.time()
muestra_probitHMC <- mod_sn$sample(
  output_dir = r'(C:\Users\Lesau\Documents\EntornosR\probit\HMC-8020)',
  data = dat,
  seed = 123,
  init = list(init_list),
  chains=1,
  iter_warmup = 10000,
  thin=2,
  iter_sampling = 2000
)
muestra_hmc <- muestra_probitHMC$summary(variables = out_var, out_stan)
tproc_hmc <- c(tproc_hmc, (proc.time()-t0)[3])
save(file=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-probit-sesgado-8020.RData)',
     muestra_probitHMC)

# (3) Modelo probit-ordenado para estimación ====    
file <- file.path('OrdProbitSN-rho01-SSVS-genq-oneInt.stan') # only observed
file <- file.path('OrdProbitSN-rho01-SSVS-genq.stan') # only observed
mod_sn <- cmdstan_model(file, compile=T) 
## (3.1) Seleccionar datos ====
## 0: urbano, 1: rural
tmp <- tibble(ym, contexto=Xide$rururb)
tmp <- tmp %>% mutate(
  ym_ord = case_when(
    is.na(ym) ~ NA,
    contexto==0 & ym < log(LPEI_urb) ~ 0,
    contexto==1 & ym < log(LPEI_rur) ~ 0,
    contexto==0 & ym > log(LPEI_urb) & ym < log(LPI_urb)  ~ 1,
    contexto==1 & ym > log(LPEI_rur) & ym < log(LPI_rur)  ~ 1,
    TRUE  ~ 2),
  ym_bin = case_when(
    ym_ord <= 1 ~ 1, # si carencias
    is.na(ym_ord) ~ NA,
    TRUE ~ 0)# no carencias
)
table(tmp$ym_ord)
# no invertir
# tmp$ym_ord <- ifelse(tmp$ym_ord == 0, 2, ifelse(tmp$ym_ord == 2, 0, 1))
set.seed(1)
datos <- get_obs(censo=F, y=tmp$ym_ord, prop=0.2)
dat <- list(
  K = K,
  N_obs = length(datos$y_obs),
  N_mis = length(datos$y_mis),
  p = ncol(datos$X_obs),
  # y_obs = datos$y_obs,
  y_obs = datos$y_obs,
  X_obs = datos$X_obs,
  X_mis = datos$X_mis,
  X_all = Xtmp[datos$master, ],
  group = datos$group[datos$observed],
  group_mis = datos$group[datos$missing],
  group_all = region[datos$master],
  # probits
  N0 = length(which(datos$y_obs==0)),
  N1 = length(which(datos$y_obs==1)),
  N2 = length(which(datos$y_obs==2)),
  n0 = which(datos$y_obs==0),
  n1 = which(datos$y_obs==1),
  n2 = which(datos$y_obs==2),
  # ssvs. Stan parametriza la normal con la desv. estándar
  tau=1/300,
  c=3000, # sqrt(10 * 300**2)
  epsilon=0.1
)

tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs)
r0 <- lm(y_obs ~ ., data = tmp)
rho0 <- 0.5
delta10 <- 0.25
mu0  <- 0
out_var <- c('rho', 'delta1', 'mu', 'b')
## (3.2) Valores inciales ====
set.seed(2)
init_list <- list(
  b = get_coefs(r0, F),
  pr=rep(0.5, ncol(datos$X_mis)),
  # mu=mu0,
  mu=rep(mu0, times=K),
  rho=rep(rho0, times=K),
  delta=delta10,
  w_tilde=rnorm(K, 0.5, sd=0.1)
)
## (3.3) Correr el modelo con HMC ====
t0 <- proc.time()
muestra_probitordHMC <- mod_sn$sample(
  output_dir = r'(C:\Users\Lesau\Documents\EntornosR\probit-ordenado\HMC-8020)',
  data = dat,
  seed = 123,
  init = list(init_list),
  chains=1,
  iter_warmup = 10000,
  thin=2,
  iter_sampling = 2000
)
muestra_hmc <- muestra_probitordHMC$summary(variables = out_var, out_stan)
tproc_hmc <- c(tproc_hmc, (proc.time()-t0)[3])
save(file=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-probit-ordenado-sesgado-8020.RData)',
     muestra_probitordHMC)

writeLines(paste0(tproc_hmc, collapse = '~'),
           con=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-tiempos-8020.txt)')

cat('Done!\n')

# (4) cargar datos ====
# rm(list=ls())
# load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR-8020\muestraHMC-log-normal-sesgado-estimacion.RData)')
# out_var <- c('rho','sigma', 'mu', 'b')
# muestra_hmc <- muestra_HMC$summary(variables = out_var, out_stan)
# print(muestra_hmc, n=40)
# muestra_HMC$output_files() # puede vivir en documents sin problema :)
