
# =============== Analisis a nivel hogares del ictpc =================

# parece que este no esta completo...
# creo que le vamos a dar matarile

rm(list=ls())
gc()

# (0) Cargar librerías ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(cmdstanr)
library(kableExtra)
# library(showtext)
# showtext_auto()
# showtext_opts(dpi = 300)
# rebuild_cmdstan()

out <- posterior::default_summary_measures()
myquant <- function(a){
  b <- quantile(a, c(0.025, 0.975))
  return(b)
}
out[5] <- 'myquant'

# (1) Datos para el modelo ====
load('./Xyregion-CDMX.Rdata')
source('./funciones-final.R')

muni <- c(
  'Azcapotzalco', 'Coyoacán', 'Cuajimalpa de Morelos', 'Gustavo A. Madero',
  'Iztacalco', 'Iztapalapa', 'La Magdalena Contreras', 'Milpa Alta',
  'Álvaro Obregón', 'Tláhuac', 'Tlalpan', 'Xochimilco', 
  'Benito Juárez','Cuauhtémoc', 'Miguel Hidalgo', 'Venustiano Carranza') %>% factor

LPI_urb <- 4564.97
LPI_rur <- 3296.92
LPEI_urb <- 2354.65
LPEI_rur <- 1800.55

## (1.1) Definir areas pequeñas ====
group <- as.numeric(region)
K <- length(unique(as.numeric(region)))

## (1.2) Covariables categoricas (no ordenadas) a matriz diseño ====
for(i in 1:ncol(Xcat)){
  Xcat[, i] <- as.factor(as.integer(Xcat[, i]))
}
tail(tibble(Xcat)) # ya todas son factor
Xcat <- model.matrix( ~ 0 + . , data=Xcat) # generar a matriz diseño

## (1.3) Remover combinaciones lineales y correlaciones ====
# Xtmp <- cbind(Xcon, Xbin)
# a <- caret::findLinearCombos(Xtmp)
# b <- caret::findCorrelation(cov(Xtmp), verbose=T, exact=T, cutoff = 0.8)
# Xtmp <- Xtmp[, -c(a$remove, b)]

## (1.4) Componentes principales ====
Eigen <- eigen(cov(Xcon))
cumsum(Eigen$values) / ncol(Xcon) # 19: 90%, 26: 95%
Xpc <- Xcon %*% Eigen$vectors[, 1:26]
a <- caret::findLinearCombos(Xbin)$remove
Xtmp <- cbind(Xpc, Xbin[, -a])

# (2) Modelo para estimación ====
train_test <- c(TRUE, FALSE)
modelos <- c('log-normal', 'probit', 'probit-ordenado')
tt <- T
m <- modelos[1]
m <- modelos[2]
m <- modelos[3]

for(tt in train_test){
  
  for(m in modelos){
    
if(m == 'log-normal'){
  file <- file.path('LogSN-rho01-SSVS-genq.stan') # only observed
}
if(m == 'probit'){
  file <- file.path('ProbitSN-rho01-SSVS-genq.stan') # only observed
}
if(m == 'probit-ordenado'){
  file <- file.path('OrdProbitSN-rho01-SSVS-genq.stan') # only observed
  file <- file.path('OrdProbitSN-rho01-SSVS-genq-design.stan') # only observed
}
    
mod_sn <- cmdstan_model(file, compile=T) 

## (2.1) Seleccionar datos ====
# 0: urbano, 1: rural
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
    TRUE ~ 0 # no carencias
  ) )

if(m == 'log-normal'){
  datos <- get_obs(none=tt, y=ym, region=region, X=Xtmp, prop = 0.2)     # 80 train 20 test
}
if(m == 'probit'){
  datos <- get_obs(none=tt, y=tmp$ym_bin, region=region, X=Xtmp, prop=0.2)     # 80 train 20 test
}
if(m == 'probit-ordenado'){
  datos <- get_obs(none=tt, y=tmp$ym_ord, region=region, X=Xtmp, prop=0.2)     # 80 train 20 test
}

if(tt){
  Xtmp_tt <- Xtmp
  region_tt <- group
} else{
  Xtmp_tt <- Xtmp[datos$master, ]
  region_tt <- group[datos$master]
}
  
dat <- list(
  K = K,
  N_obs = length(datos$y_obs),
  N_mis = length(datos$y_mis),
  p = ncol(datos$X_obs),
  y_obs = datos$y_obs,
  X_obs = datos$X_obs,
  X_mis = datos$X_mis,
  X_all = Xtmp_tt, # [datos$master, ]
  group = datos$group[datos$observed],
  group_mis = datos$group[datos$missing],
  group_all = region_tt, # [datos$master]
  # probits
  N0 = length(which(datos$y_obs==0)),
  N1 = length(which(datos$y_obs==1)),
  N2 = length(which(datos$y_obs==2)),
  n0 = which(datos$y_obs==0),
  n1 = which(datos$y_obs==1),
  n2 = which(datos$y_obs==2),
  # ssvs. Stan parametriza la normal con la desv. estándar
  tau=1/300,  # El 'pico' es N(0, sd=0.01)
  c=300*5,     # La 'losa' es N(0, sd=100*0.01) = N(0, sd=10)
  epsilon=0.01 # aceptamos |beta_{i}| > 0.05; o 1/2
)

# delta(dat$tau**2, dat$c**2)
# delta(dat$tau**2, 300*10)

## (2.2) Valores iniciales ====

if(m == 'log-normal'){
  tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs) # Xpc[datos$master, ]
  r1 <- lm(y_obs ~ ., data=tmp)
  r0 <- sn::selm(y_obs ~ 1 , data=tmp)
  mu0 <- coef(r0, param.type = 'DP')[1]
  rho0 <- lambda2rho(coef(r0, param.type = 'DP')[3])
  sigma0 <- coef(r0, param.type = 'DP')[2]
  out_var <- c('rho', 'sigma', 'mu', 'b')
}
if(m == 'probit'){
  tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs) # Xpc[datos$master, ]
  r1 <- glm(y_obs ~ ., data = tmp, family = binomial(link='probit'))
  rho0 <- 0.5
  sigma0 <- 1
  mu0  <- 0
  out_var <- c('rho', 'mu', 'b')
}
if(m == 'probit-ordenado'){
  tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs) # Xpc[datos$master, ]
  r1 <- lm(y_obs ~ ., data = tmp)
  sigma0 <- 1
  rho0 <- 0.5
  delta1 <- 0.2
  mu0  <- 0
  out_var <- c('rho', 'delta1', 'mu', 'b')
}

set.seed(2)
init_list <- list(
  b = get_coefs(r1, F),
  pr=rep(0.5, ncol(datos$X_mis)),
  sigma=sigma0,
  mu=rep(mu0, times=K),
  rho=rep(rho0, times=K),
  z_std=rnorm(K, 0.5, sd=0.1)
)

## (2.3) Correr el modelo con VB ====
t0 <- proc.time()
fit_vb <- mod_sn$variational(
  data = dat,
  seed = 123,
  init = list(init_list),
  draws = 500,
  adapt_iter=1000,
  adapt_engaged=T, # se encontró óptimo eta=0.1
  eta=0.1, # 0.05
  iter=75000,
  tol_rel_obj=1e-6,
  algorith='meanfield', # fullrank
)
fit_hmc <- mod_sn$sample(
  data = dat,
  seed = 123,
  init = list(init_list),
  chains=1,
  iter_warmup = 1000,
  iter_sampling = 1000
)
tproc_hmc <- proc.time()-t0 # Tiempo de proceso
fit_hmc$summary()
muestra_hmc <- fit_hmc$summary(variables = out_var, out)
print(muestra_hmc, n=35)

# podemos hacer esto mas artesanal, es decir
# sin ser automatizado
getwd()
nombre <- paste0('C:/Maestría_CP/Tesis/Documento-Tesis-SAOM/CDMX/muestra-VB-', m, '-', tt, '.RData')
# save(fit_vb, file=nombre)
  }
}

## (2.4) Precisión y tiempo (ajustados) ====
cat("\n","Tiempo de ejecución:",tproc_vb[3],"(elapsed)\n")
out <- which(!is.na(ym))
fitted_vb <- apply(fit_vb$draws('z_all'), 2, mean)
metrics(fitted_vb[out], datos$y_test)
metrics(fitted(r1), datos$y_test)

out <- datos$idx
fitted_vb <- findInterval(fitted_vb, 0)
table(fitted_vb[out], datos$y_obs)

print(muestra_vb, n=30)
out <- datos$idx
fitted_vb <- findInterval(fitted_vb, c(0, 1.44))
table(fitted_vb[datos$master][out], datos$y_obs)

## (2.5) Gráfico para las estimaciones de rho, sigma, mu ====
library(bayesplot)
g_rho <- mcmc_intervals(fit_vb$draws('rho')) + 
  scale_y_discrete(labels = parse(text = paste0("rho[", 1:16, "]")))
g_sigma <- mcmc_hist(fit_vb$draws('sigma'), alpha = 0.2) + 
  labs(x = parse(text = 'sigma'))
g_mu <- mcmc_intervals(fit_vb$draws('mu')) +
  scale_y_discrete(labels = parse(text = paste0("mu[", 1:16, "]")))

x11(); print(g_rho)
x11(); print(g_sigma)
x11(); print(g_mu)

ggsave(filename='int_rho.png', plot=g_rho, width=15, heigh=15, units='cm')
ggsave(filename='hist_sigma.png', plot=g_sigma, width=15, heigh=15, units='cm')
ggsave(filename='int_mu.png', plot=g_mu, width=15, heigh=15, units='cm')

# (3) Estimación a posteriori ====

## (3.1) Estimación de beta ====
beta_ <- fit_vb$draws("b")
Media <- apply(beta_, 2, mean)
Error <- apply(beta_, 2, sd)
Lo025 <- apply(beta_, 2, quantile, 0.025)
Up975 <- apply(beta_, 2, quantile, 0.975)
Frecuencia <- 100*apply(fit_vb$draws("m_ind"), 2, mean)
# Frecuencia <- 100*apply(fit_vb$draws("pr"), 2, mean)
Parametro <- paste0("$\\beta_{", 1:ncol(datos$X_obs), "}$")

beta_out <- tibble(Parametro, Media, `Error est.`=Error, `0.025 \\%`=Lo025, `0.975 \\%`=Up975, `Frecuencia \\%`=Frecuencia)
print(beta_out, n=Inf)
beta_out %>% filter(`Frecuencia \\%` > 75) %>%
  print(n=Inf)

table(beta_out$`Frecuencia \\%`)

beta_out %>% filter(`Frecuencia \\%` > 75) %>%
  kable(
    format='latex', escape=F, digits = 4,
    label='tab:estimaciones-beta',
    caption = 'Elaboración propia basada en la muestra \\underline{a posteriori}.'
  )


posterior_prob(fit_vb)[1:3, 1:2]
sort(marginal_pi(fit_vb))

## (3.2) Covariables seleccionadas ====
beta_out2 <- beta_out %>%
  rename(Frecuencia=`Frecuencia \\%`) %>%
  mutate(Numero = 1:nrow(beta_out)) %>%
  select(Numero, Frecuencia) %>% 
  arrange(Frecuencia)
print(beta_out2, n=15)

bottom <- beta_out2$Numero[1:7] # menores del 50%
top <- beta_out2 %>% filter(Frecuencia>99) %>% select(Numero)
top <- top$Numero
colnames(Xtmp)[bottom]
colnames(Xtmp)[top]

## (3.3) Estimacion de rho, sigma y mu ====
rho_ <- fit_vb$draws("rho")
Media_rho <- apply(rho_, 2, mean)
Error_rho <- apply(rho_, 2, sd)
Lo025_rho <- apply(rho_, 2, quantile, 0.025)
Up975_rho <- apply(rho_, 2, quantile, 0.975)

sigma_ <- fit_vb$draws("sigma")
Media_sigma <- apply(sigma_, 2, mean)
Error_sigma <- apply(sigma_, 2, sd)
Lo025_sigma <- apply(sigma_, 2, quantile, 0.025)
Up975_sigma <- apply(sigma_, 2, quantile, 0.975)

mu_ <- fit_vb$draws("mu")
Media_mu <- apply(mu_, 2, mean)
Error_mu <- apply(mu_, 2, sd)
Lo025_mu <- apply(mu_, 2, quantile, 0.025)
Up975_mu <- apply(mu_, 2, quantile, 0.975)

# Cuando se usa pooling
# mu0_ <- fit_vb$draws("mu_0")
# Media_mu0 <- apply(mu0_, 2, mean)
# Error_mu0 <- apply(mu0_, 2, sd)
# Lo025_mu0 <- apply(mu0_, 2, quantile, 0.025)
# Up975_mu0 <- apply(mu0_, 2, quantile, 0.975)
# 
# sigma_mu_ <- fit_vb$draws("sigma_mu")
# Media_sigma_mu <- apply(sigma_mu_, 2, mean)
# Error_sigma_mu <- apply(sigma_mu_, 2, sd)
# Lo025_sigma_mu <- apply(sigma_mu_, 2, quantile, 0.025)
# Up975_sigma_mu <- apply(sigma_mu_, 2, quantile, 0.975)

Parametro <- c(
  paste0("$\\rho_{", 1:16, "}$"),
  paste0("$\\sigma$"),
  paste0("$\\mu$")
  # paste0("$\\mu_{", 1:16, "}$")
  #paste0("$\\mu_{0}$"),      # pooling
  #paste0("$\\sigma_{\\mu}$") # pooling
)

param_out <- tibble(
  Parametro=Parametro,
  Media=c(Media_rho, Media_sigma, Media_mu), # , Media_mu0, Media_sigma_mu
  `Error est.`=c(Error_rho, Error_sigma, Error_mu), # , Error_mu0, Error_sigma_mu
  `0.025 \\%`=c(Lo025_rho, Lo025_sigma, Lo025_mu), # , Lo025_mu0, Lo025_sigma_mu
  `0.975 \\%`=c(Up975_rho, Up975_sigma, Up975_mu), # , Up975_mu0, Up975_sigma_mu
)

print(param_out, n=Inf)
param_out %>% kable(
  format='latex', escape = F, digits = 4,
  label='tab:estimaciones',
  caption = 'Elaboración propia basada en la muestra \\underline{a posteriori}.'
)

K <- 16
# par1 <- Parametro[1:(K+1)]
# par2 <- c(NA, Parametro[(K+1):length(Parametro)])
Media_Error <- c()
Lo_Up <- c()
for(k in 1:(1+2*K)){
  a0 <- sprintf('%.4f', param_out$Media[k])
  a1 <- sprintf('%.4f', param_out$`Error est.`[k])
  Media_Error[k] <- paste('\\makecell[l]{', a0, '\\\\', '(', a1, ')', '}', sep='')
  b0 <- sprintf('%.4f', param_out$`0.025 \\%`[k])
  b1 <- sprintf('%.4f', param_out$`0.975 \\%`[k])
  Lo_Up[k] <- paste('(', b0, ', ', b1, ')', sep='')
}

param_out2 <- tibble(
  Parámetro=Parametro,
  `\\makecell[l]{Media\\\\(Error est.)}`= Media_Error,
  `(2.5 \\%, 97.5 \\%)`= Lo_Up
)

kable(
  list(param_out2[1:(1+K),], param_out2[(2+K):(2*K), ]),
  format='latex', escape = F,
  label='tab:estimaciones',
  caption = 'Elaboración propia basada en la muestra \\underline{a posteriori}.'
)


## (3.4) Medidas de desigualdad ====

# 0: ámbito urbano
# 1: ámbito rural
library(ineq)
library(dineq)

## (3.5) Gini-Lorenz a posteriori ====
out <- which(is.na(ym))
fitted_vb <- exp(apply(fit_vb$draws('y_all'), 2, mean))
# fitted_vb <- exp(fitted_vb)

dftmp <- data.frame(
  Ajustados=fitted_vb[out],
  Contexto=Xide$rururb[out],
  Factor=Xide$factor[out]
)

lc <- Lc(dftmp$Ajustados, dftmp$Factor)
lc_urb <- Lc(dftmp$Ajustados[dftmp$Contexto==0], dftmp$Factor[dftmp$Contexto==0])
lc_rur <- Lc(dftmp$Ajustados[dftmp$Contexto==1], dftmp$Factor[dftmp$Contexto==1])

round(Gini(dftmp$Ajustados),4)
round(gini.wtd((dftmp$Ajustados), dftmp$Factor),4)
gini_decomp((dftmp$Ajustados), dftmp$Contexto, dftmp$Factor)$gini_group$gini_group
#$gini_group$gini_group
# 0         1 
# 0.3729112 0.2986803 

len0 <- c(nrow(dftmp), table(dftmp$Contexto))+1L
len1 <- c()
for(k in seq_along(len0)){
  len1 <- c(len1, rep(c('Ambos', 'Urbano', 'Rural')[k], times=len0[k]))
}
len1 <- factor(len1)

LG_all <- data.frame(
  p=c(lc$p, lc_urb$p, lc_rur$p), L=c(lc$L, lc_urb$L, lc_rur$L), Ámbito=len1)

LG_plot <- ggplot(LG_all, aes(x=p, y=L, group=Ámbito, linetype=Ámbito, color=Ámbito, label=Ámbito)) + 
  geom_line(linewidth=1.2) + 
  # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  coord_fixed(ratio = 1) +          # ensures a square panel
  geom_abline(intercept = 0, slope=1, linewidth=1.2, color='gray75', lty=1) + 
  theme_bw() + 
  theme(text=element_text(family='serif', size=11), legend.position = 'top')

x11(); plot(LG_plot)
ggsave(filename='Lorenz-CDMX.png', plot=LG_plot, width=15, heigh=15, units='cm')

## (2.4) Porcentaje de personas con las LPI y LPEI ====
out <- which(is.na(ym))
fitted_vb <- exp(apply(muestra_vb$draws('y_all'), 2, mean))
fitted_vb <- findInterval(apply(fit_vb$draws('z_all'), 2, mean), 0)
fitted_vb <- findInterval(apply(fit_vb$draws('z_all'), 2, mean), c(0, 1.41))
dftmp <- tibble(Ajustados=fitted_vb[out],
                Ámbito=Xide$rururb[out],
                Factor=Xide$factor[out],
                Tamaño=Xide$tam_hog[out],
                Alcaldia=muni[region][out])
# tipo log normal
dftmp <- dftmp %>%
  mutate(LPI = case_when(
    `Ámbito` == 0 & Ajustados < LPI_urb ~ 1L,
    `Ámbito` == 1 & Ajustados < LPI_rur ~ 1L,
    TRUE                                ~ 0L)) %>%
  mutate(LPEI = case_when(
    `Ámbito` == 0 & Ajustados < LPEI_urb ~ 1L,
    `Ámbito` == 1 & Ajustados < LPEI_rur ~ 1L,
    TRUE                                ~ 0L))
#

# tipo probit
dftmp2 <- dftmp %>%
  group_by(Alcaldia) %>%
  summarise(Total = sum(Factor*Tamaño),
            Total_LPI = sum(Factor*Ajustados*Tamaño), Porcentaje_LPI = 100*Total_LPI/Total)
dftmp2
sum(dftmp2$Total_LPI) / sum(dftmp2$Total) # 25.7 % 
#

# tipo probit ordenado
dftmp2 <- dftmp %>%
  mutate(LPEI=case_when(Ajustados == 0 ~ 1, TRUE ~ 0),
         LPI =case_when(Ajustados == 1 ~ 1, TRUE ~ 0)) %>%
  group_by(Alcaldia) %>%
  summarise(Total = sum(Factor*Tamaño),
            Total_LPI = sum(Factor*LPI*Tamaño), Total_LPEI = sum(Factor*LPEI*Tamaño),
            Porcentaje_LPI = 100*Total_LPI/Total, Porcentaje_LPEI = 100*Total_LPEI/Total)
dftmp2
sum(dftmp2$Total_LPI) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPEI) / sum(dftmp2$Total) # 25.7 % 
#

dftmp2 <- dftmp %>%
  group_by(Alcaldia) %>%
  summarise(Total = sum(Factor*Tamaño),
            Total_LPI = sum(Factor*LPI*Tamaño), Total_LPEI = sum(Factor*LPEI*Tamaño),
            Porcentaje_LPI = 100*Total_LPI/Total, Porcentaje_LPEI = 100*Total_LPEI/Total)

sum(dftmp2$Total_LPI) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPEI) / sum(dftmp2$Total) # 4.5 % 

dftmp2 %>%
  select(Alcaldia, starts_with('Porcentaje')) %>%
  kable(format = 'latex', digits=3,
        label='tab:porcentajes',
        caption='Elaboración propia basada en la muestra \\underline{a posteriori}.')



# save(muni, fit_vb, dftmp2, file='./fit_vb_CPV2020.RData')

##  Cuantificar sesgo ====

# rm(list=ls())
# ruta <- r'(C:\Maestría_CP\Tesis\Articulo)'
# setwd(ruta)
# load('./Xyregion-CDMX.Rdata')
# load('./fit_vb_CPV2020.RData')
# source('./funciones-final.R')
idx <- which(!is.na(ym))
fitted_vb <- apply(fit_vb$draws('y_all'), 2, mean)

dftmp <- tibble(ICTPC_obs = exp(ym), ICTPC_fit=exp(fitted_vb), Alcaldia=muni[region])
dftmp <- dftmp %>%
  group_by(Alcaldia) %>%
  summarise(mean_ICTPC_obs = mean(ICTPC_obs, na.rm=T), # cambiar por median, esta curioso
            mean_ICTPC_fit = mean(ICTPC_fit, na.rm=T),
            diff_ICTPC = abs(mean_ICTPC_obs-mean_ICTPC_fit))

mean(dftmp$diff_ICTPC)

dftmp %>%
  kable(format = 'latex', digits=3,
        label='sesgo',
        caption='Elaboración propia basada en la muestra \\underline{a posteriori}.')

# (4) Descriptivos ====

## (4.1) Tablas de estadísticos ====

datos <- list()
datos$master <- which(!is.na(ym))
dftmp <- data.frame(Observados=ym[datos$master],
                    Ámbito=Xide$rururb[datos$master],
                    Factor=Xide$factor[datos$master],
                    Alcaldia=muni[region[datos$master]])

Res_ambito <- dftmp %>%
  group_by(Ámbito) %>%
  summarise(Minimo=min(exp(Observados)),
            Mediana = median(exp(Observados)),
            Media = mean(exp(Observados)),
            Maximo=max(exp(Observados)),
            `Desv. est.`=sd(exp(Observados)),
            Conteo=n(), 
            .groups = "drop")

Res_ambito %>%
  kable(digits=3, format='latex')

Ambito_agrupado <- dftmp %>%
  mutate(Ámbito=factor(Ámbito, labels=c('Urbano', 'Rural'))) %>%
  group_by(Ámbito) 

Ambito_resumen %>% kable(format = "latex",
                         booktabs = F,
                         caption = "Resumen por Alcaldía y Ámbito",
                         escape = TRUE,
                         linesep = "") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  collapse_rows(columns = 1, valign = "middle")

Ambito_resumen %>% kable(format = "latex",
                         booktabs = TRUE,
                         caption = "Resumen por Alcaldía y Ámbito",
                         escape = TRUE,
                         linesep = "") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  collapse_rows(columns = 1, valign = "middle")


dftmp %>% group_by(Alcaldia) %>%
  summarise(Minimo=min(exp(Observados)),
            Mediana = median(exp(Observados)),
            Media = mean(exp(Observados)),
            Maximo=max(exp(Observados)),
            `Desv. est.`=sd(exp(Observados)),
            Conteo=n(), 
            .groups = "drop")

table_out <- dftmp %>%
  mutate(Ámbito=factor(Ámbito, labels=c('Urbano', 'Rural'))) %>%
  group_by(Alcaldia, Ámbito) %>%
  summarise(Minimo=round(min(exp(Observados)), 3),
            Mediana = round(median(exp(Observados)), 3),
            Media = round(mean(exp(Observados)), 3),
            Maximo=round(max(exp(Observados)), 3),
            `Desv. est.`=round(sd(exp(Observados)), 3),
            Conteo=n(), 
            .groups = "drop") %>%
  complete(Alcaldia, Ámbito,
           fill = list(
             Minimo = NA_real_,
             Mediana = NA_real_,
             Media = NA_real_,
             Maximo = NA_real_,
             `Desv. est.` = NA_real_,
             Conteo = 0)) 
# arrange(Alcaldia, `Ámbito`) %>% #Order rows using column values

# install.packages('kableExtra')
library(kableExtra)
latex_table <- table_out %>%
  kable(format = "latex",
        booktabs = TRUE,
        caption = "Resumen por Alcaldía y Ámbito",
        escape = TRUE,
        linesep = "") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  collapse_rows(columns = 1, valign = "middle")

latex_table

## (4.2) Log-ICTPC grafico de violin ====

Ambito <- ggplot(Ambito_agrupado, aes(x = Ámbito, y = Observados, fill=Ámbito)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.25, outlier.shape = 2, outlier.size = 4) +
  geom_jitter(width = 0.12, alpha = 0.2, size = 1) +
  labs(x = "Ámbito", y = 'log - ingreso corriente total pér cápita', title = "") +
  theme_minimal() + 
  theme(legend.position = 'none', text = element_text(family = "serif", size=11))
x11(); plot(Ambito)
ggsave(filename='./Box-Ambito.png', plot=Ambito,
       width = 15, height = 15, bg = "transparent", units='cm')


Ing <- ggplot(dftmp, aes(x = Ámbito, y = Observados, fill=`Ámbito`)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.25, outlier.shape = 2) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1) +
  labs(x = "Alcaldia", y = 'log-ICTPC', title = "") +
  theme_minimal() + 
  theme(legend.position = 'none', text = element_text(family = "serif"))
x11(); print(Ing)

## (4.3) Lista de covariables ====

library(readxl)
library(tidyr)
library(stringr)

load('./Xyregion-CDMX.Rdata')
descripcion <- read_excel('descripcion_variables.xlsx', range='A1:D308')
head(descripcion)

descripcion <- descripcion %>%
  # borrar espacios invisibles y convertir "" a NA
  mutate(Variable = na_if(str_squish(Variable), "")) %>%
  fill(Variable, .direction = "down") %>%
  mutate(Valores=if_else(Valores=='(1…8)', '(0…999999999)', Valores)) %>%
  mutate(Valores=factor(Valores))
descripcion

nomcon <- descripcion %>%
  filter(Variable %in% colnames(Xcon)) %>%
  distinct(Variable, .keep_all = T)
dim(Xcon)
dim(nomcon)

nombin <- descripcion %>%
  filter(Variable %in% colnames(Xbin)) %>%
  distinct(Variable, .keep_all = T)
dim(Xbin)
dim(nombin)

nomcat <- descripcion %>%
  filter(Variable %in% colnames(Xcat)) %>%
  distinct(Variable, .keep_all = T)
dim(Xcat)
dim(nomcat)

nomall <- rbind(
  nomcon %>% arrange(Variable),
  nombin %>% arrange(Variable),
  nomcat %>% arrange(Variable)
)

# Tipo <- c(rep('Con.', ncol(Xcon)), rep('Bin.', ncol(Xbin)), rep('Cat.', ncol(Xcat)))
Tipo <- c(rep('Continua', ncol(Xcon)), rep('Binaria', ncol(Xbin)), rep('Categórica', ncol(Xcat)))
nomall <- cbind(nomall, Tipo)

tabla_cov <- nomall %>%
  select(!c('Valores', 'Etiquetas')) %>%
  relocate(Tipo) %>%
  tibble %>%
  filter(Variable %in% colnames(Xtmp))

tabla_cov
print(nomcon, n=Inf)

## (4.4) Covariables y beta ====

idx <- which(beta_out$`Frecuencia \\%` > 75)
Parametro <- paste0("$\\beta_{", 1:length(idx), "}$")
Parametro <- paste0('\\makecell[l]{', '$\\beta_{', 1:length(idx), '}$',
                    '\\\\ (', sprintf('%.2f', beta_out$`Frecuencia \\%`[idx]), '\\%)}')

Media_Error <- c()
Lo_Up <- c()
beta_out
for(k in 1:length(idx)){
  a0 <- sprintf('%.4f', beta_out$Media[k])
  a1 <- sprintf('%.4f', beta_out$`Error est.`[k])
  Media_Error[k] <- paste('\\makecell[l]{', a0, '\\\\', '(', a1, ')', '}', sep='')
  b0 <- sprintf('%.4f', beta_out$`0.025 \\%`[k])
  b1 <- sprintf('%.4f', beta_out$`0.975 \\%`[k])
  Lo_Up[k] <- paste('(', b0, ', ', b1, ')', sep='')
}

param_out2 <- tibble(
  `\\makecell[l]{Parámetro\\\\(Frecuencia)}`=Parametro,
  `\\makecell[l]{Media\\\\(Error est.)}`= Media_Error,
  `\\makecell[l]{Intervalo c. \\\\ (2.5 \\%, 97.5 \\%)}`= Lo_Up,
  Descripción = tabla_cov$Descripción[idx]
)

param_out2_ <- kable(
  param_out2,
  format='latex', escape = F,
  label='covariables-beta',
  caption = 'Estimaciones y lista de covariables incluidas, únicamente
  se muestran aquellas covariables con frecuencias de aparición mayores al 75\\%.
  Fuente: elaboración propia.Elaboración propia basada en la muestra \\underline{a posteriori}.'
)

write.table(file='covariables-beta.txt', param_out2_, row.names = F, col.names = F)


## Mas graficos ====

x11()
mcmc_intervals(fit_vb$draws("pr"))
mcmc_intervals(fit_vb$draws("rho"))
mcmc_hist(fit_vb$draws("rho"))
mcmc_hist(fit_vb$draws("z"))
mcmc_hist(fit_vb$draws("mu"))
mcmc_hist(fit_vb$draws("sigma"))

# ajustados, todo el censo
fitted_vb <- apply(fit_vb$draws("y_mis"), 2, mean)[!is.na(ym)]
metrics(fitted_vb, datos$ytest)
metrics(fitted(r1), datos$ytest)
x11(); plot(fitted_vb, datos$ytest); points(fitted(r1), datos$ytest, col='red', pch=19)
#

fitted_vb <- apply(fit_vb$draws("ymis"), 2, mean)
metrics(fitted_vb, datos$ytest)
metrics(ymishat, datos$ytest)
x11(); plot(fitted_vb, datos$ytest); points(ymishat, datos$ytest, col='red', pch=19)

# write.csv(tproc_vb, file='./tiempo_VB_iter125k_grad1_elbo100_eta01.csv')
# save(muestra_VB, file='./muestraVB_iter125k_grad1_elbo100_eta01.Rdata')


## (4.5) Gráfico de distribución del ingreso (histogramas) ----
umbral <- 210 # filtrar dos o tres ingresos muy pequeños
umbral <- 500
tmp <- which(!is.na(ym) & ym>log(umbral)) # o solo 10
dftmp <- data.frame(y = ym[tmp], mun = muni[region[tmp]],
                    amb = factor(Xide$rururb[tmp], labels=c('Urbano', 'Rural')))
y_hist <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray80', fill=4) +
  geom_density(linewidth=1.2, color='gray40') + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif"))
x11(); plot(y_hist)
ggsave(r"(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX\dist-log-ict.pdf)", 
       plot = y_hist, width = 15, height = 15, units = "cm")

y_hist2 <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray80', fill=4) +
  geom_density(linewidth=1, color='gray40') + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif")) + 
  facet_wrap(~ mun, ncol = 4, nrow = 4)
x11(); plot(y_hist2)
ggsave("dist-log-ict-alcaldias.png", plot = y_hist2,
       # device = grDevices::cairo_pdf,
       width = 15, height = 15, units = "cm")

y_hist3 <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray80', fill=4, bins = 30) +
  geom_density(linewidth=1, color='gray50') + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif")) + 
  facet_wrap(~ amb, ncol = 2, nrow = 1, scales = 'fixed')
x11(); plot(y_hist3)
ggsave(r"(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX\dist-log-ict-ambito.pdf)", 
       plot = y_hist3, width = 15, height = 15, units = "cm")


## Gini-Lorenz antes del analisis posterior----
datos <- list()
datos$master <- which(!is.na(ym))
dftmp <- data.frame(Observados=ym[datos$master],
                    Ámbito=Xide$rururb[datos$master],
                    Factor=Xide$factor[datos$master],
                    Alcaldia=muni[region[datos$master]])

lc <- Lc(exp(dftmp$Observados), dftmp$Factor)
lc_urb <- Lc(exp(dftmp$Observados[dftmp$Ámbito==0]), dftmp$Factor[dftmp$Ámbito==0])
lc_rur <- Lc(exp(dftmp$Observados[dftmp$Ámbito==1]), dftmp$Factor[dftmp$Ámbito==1])
round(Gini(exp(dftmp$Observados)),4)
round(gini.wtd(exp(dftmp$Observados), dftmp$Factor),4)
gini_decomp(exp(dftmp$Observados), dftmp$Ámbito, dftmp$Factor)
# $gini_group
# $gini_group$gini_group
# 0         1 
# 0.4804783 0.3116281 



x11();
plot(lc,main="Indice de Gini y curva de Lorenz para EdoMex",lwd=2)
lines(lc_urb$p, lc_urb$L,col="red",lty=2,lwd=2)
lines(lc_rur$p, lc_rur$L,col="blue",lty=2,lwd=2)
legend(0.05,0.9,legend=c("Rural (0.3856)","Urbano (0.3598)","General (0.3656)"),
       col=c("Blue","Red","Black"),lty=c(2,2,1),lwd=2,bty="n",cex=1.1)
grid()



library(ggplot2)
library(scales)
factores <- ggplot(Xide, aes(x = as.factor(rururb), y = factor, fill=as.factor(rururb))) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.25, outlier.shape = 2, outlier.size = 4) +
  geom_jitter(width = 0.08, alpha = 0.05, size = 1) +
  labs(x = "Ámbito", y = bquote(log[10]- ~factor~de~expansión), title = "") +
  theme_minimal() + 
  theme(legend.position = 'none', text = element_text(family = "serif", size=11)) + 
  scale_y_log10(labels = label_log(digits = 2)) +
  annotation_logticks(sides='l')
x11();plot(factores)
ggsave(filename='factores_por_ambito.pdf', plot=factores, width=20, height = 20, units = 'cm')



