
# Este código procesa exclusivamente las simulaciones

rm(list=ls())
# install.packages("devtools")

library(tidyverse)
# library(kableExtra)
library(knitr)

which_closer <- function(a){
  if(is.character(a[1])){
    # quien se acerca mas al promedio de los interceptos
    b <- a[1]
    b <- mean(as.numeric(str_split_1(b, pattern = '\\|')))
    a <- c(b, as.numeric(a[2]), as.numeric(a[3]))
  }
  target <- a[1]
  closer <- numeric(1L)
  dista1 <- abs(a[2]-a[1])
  dista2 <- abs(a[3]-a[1])
  if(dista1 < dista2){
    closer <- 1L # vb is closer
  } else if(dista2 < dista1){
    closer <- 2L # hmc is closer
  } else{
    closer <- 1L # benefit of doubt to vb
  }
  return(closer)
}

which_closer(c(0, 0.1, 0.6))
which_closer(c('1.0|2.0|3.0', 0.8, 0.9))
which_closer(c('0.5', 0.8, 0.9))

metrics <- function(fitted, actual){
  mae <- mean(abs(actual - fitted))
  mse <- mean((actual - fitted)**2)
  mape<- 100*mean(abs((fitted - actual)/actual))
  corr<- cor(actual, fitted)
  return(c(corr, mae, sqrt(mse), mape))
}

metrics_cat <- function(fitted, actual, binary=TRUE){
  cm <- table(fitted, actual)
  total <- length(actual)
  if(binary){
    TN <- cm[1, 1]
    FP <- cm[1, 2]
    FN <- cm[2, 1]
    TP <- cm[2, 2]
    accuracy <- (TP + TN) / total # Overall proportion of correct predictions
    precision <- TP / (TP + FP) # precision: Of predicted positives, how many were correct?
    TPR <- TP / (TP + FN) # recall, sensibilidad: Of actual positives, how many were correctly predicted?
    TNR <- TP / (TP + FN) # specificity: Of actual negatives, how many were correctly predicted?
    FNR <- 1 - TNR # How often the model incorrectly predicts positive
    FPR <- 1 - TPR # 	How often the model misses positive cases
    f1 <- 2 * precision * TPR / (precision + TPR) # Harmonic mean of Precision and Recall
    
    output <- c(accuracy, TPR, TNR, f1)
  } else{
    accuracy <- sum(diag(cm)) / total
    ktau <- cor(fitted, actual, method = "kendall") # kendall's tau
    output <- c(accuracy, ktau)
  }
  return(output)
}

setwd(r"(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Codigos\Simulacion-Modelos)")
master <- getwd()
carpetas <- dir()[endsWith(dir(), 'sim')]
carpetas

# navegar por las carpetas ====
folder <- carpetas[1]

for(folder in carpetas){
  
  setwd(file.path(master, folder, fsep = '/'))
  
  cat('Estamos en la carpeta:', folder, '\n')
  
  # Seccion I: tabla resumen ====
  archivos <- dir()
  archivos <- archivos[startsWith(archivos, 'muestra')]
  # archivos <- list.files(path = '.', pattern = "\\.csv$", full.names = TRUE)

file_info <- tibble(tibble(archivos)) %>%
  mutate(
    type = str_extract(archivos, "^(muestra|predict|SSVS)"),
    intercept = str_extract(archivos, "(varios|unico|no)"),
    method = str_extract(archivos, "(vb|hmc)"),
    percentage = as.integer(str_extract(archivos, "(95|50|25|5)"))
  )

Muestras <- file_info %>%
  filter(type=='muestra') %>%
  filter(percentage%in%c(5, 25)) %>% 
  arrange(percentage, desc(intercept), desc(method))

i=1
tiempos <- c()
for(i in 1:(nrow(Muestras)/2)){
  primero <- 2*i - 1
  segundo <- 2*i
  int <- as.character(Muestras[primero, 3]) # varios, unico, no: interceptos
  per <- as.character(Muestras[primero, 5]) # porcentaje de muestreo: 5% o 25%
  mod <- str_split_1(folder, pattern='-sim')[1]
  vb <- read.csv(paste0('./', Muestras[primero, 1]))
  hmc <- read.csv(paste0('./', Muestras[segundo, 1]))
  L <- nrow(vb)
  tiempos <- c(tiempos, vb[L, 2], hmc[L, 2])
  # Etiquetas para los parametros
  parametros <- c(paste0('$\\rho_{', 1:4,'}$'), '$\\sigma^{2}$','$\\delta_{1}$', 
                  paste0('$\\mu_{', 1:4,'}$'),
                  # paste0('$\\mu_{', 1:4,'}$', collapse = ' \\, '),
                  paste0('$\\mu$'),
                  paste0('$\\beta_{', 1:5,'}$'))
  names(parametros) <- c(paste0('rho', 1:4), 'sigma2', 'delta1',
                         paste0('mu', 1:4), 'mu', paste0('beta', 1:5))
  if(int=='varios'){
    parametros <- parametros[!names(parametros) %in% c('mu')]
  }
  if(int=='unico'){
    parametros <- parametros[!names(parametros) %in% c('mu1', 'mu2', 'mu3', 'mu4')]
  }
  if(int=='no'){
    parametros <- parametros[!names(parametros) %in% c('mu1', 'mu2', 'mu3', 'mu4', 'mu')]
  }
  if(mod == 'log-normal-sesgado'){
    parametros <- parametros[!names(parametros) %in% c('delta1')]
  }
  if(mod == 'probit-sesgado'){
    parametros <- parametros[!names(parametros) %in% c('sigma2', 'delta1')]
  }
  if(mod == 'probit-ordenado-sesgado'){
    parametros <- parametros[!names(parametros) %in% c('sigma2')]
  }
  # El formato se lo daremos en LaTeX, es decir, el redondeo/número de decimales
  valores_reales <- vb[-L, 2]
  vb_mp <- vb[-L, 3] # media posterior
  vb_025 <- vb[-L, 5]
  vb_975 <- vb[-L, 6]
  hmc_mp <- hmc[-L, 3]
  hmc_025 <- hmc[-L, 5]
  hmc_975 <- hmc[-L, 6]
  tmp <- data.frame(parametros, valores_reales,
                    vb_mp, vb_025, vb_975, hmc_mp, hmc_025, hmc_975)
  # identificar quien se acercó más: vb o hmc
  closer <- apply(tmp[, c(2, 3, 6)], MARGIN=1, which_closer)
  for(j in 1:length(closer)){
    if(closer[j] == 1L){
      tmp[j, 3] <- paste0('\\bfseries ', tmp[j, 3])
    } else{
      tmp[j, 6] <- paste0('\\bfseries ', tmp[j, 6])
    }
  }
  
  if(int == 'unico'){
    # tmp['mu', 2] <- 0
    # mu <- paste0(str_split_1(tmp['mu', 2],pattern = '\\|'), collapse = ', ')
    # Le ponemos la media, para quitarnos de problemas
    mu <- mean(as.numeric(str_split_1(tmp['mu', 2],pattern = '\\|')))
    tmp['mu', 2] <- mu
    # tmp$valores_reales <- as.numeric(tmp$valores_reales)
  }

  #colnames(tmp) <- c('Parámetros', 'Valor real', 'Media posterior', 'I.C. 2.5\\%', 'I.C. 97.5\\%',
  #                      'Media posterior', 'I.C. 2.5\\%', 'I.C. 97.5\\%')
  tmp <- as.matrix(tmp); colnames(tmp) <- NULL; rownames(tmp) <- NULL
  tmp <- kable(tmp, format = "latex", escape = FALSE, align = rep("r", 8)) # , digits = 4, pero lo haremos en LaTeX
  # Establecer nombres ====
  ruta0 <- r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi)'
  nombre <- paste0(mod, '-muestreo-', per, '-interceptos-', int, '.txt')
  nombre <- file.path(ruta0, 'Tablas', nombre, fsep = '\\')
  writeLines(tmp, nombre)
  Last <- length(readLines(nombre))
  writeLines(readLines(nombre)[-c(1, 2, 3, Last)], nombre)  # skip: 1 blank line, 2 \begin{tabular}, 3 \hline and last: \end{tabular}
  cat('tabla', i, 'de', nrow(Muestras)/2, '\n -- \n')
}

# Seccion II: métricas de ajuste ====
archivos <- dir()
archivos <- archivos[startsWith(archivos, "muestra") | startsWith(archivos, "predict")]

file_info <- tibble(tibble(archivos)) %>%
  mutate(
    type = str_extract(archivos, "^(muestra|predict|SSVS)"),
    intercept = str_extract(archivos, "(varios|unico|no)"),
    method = str_extract(archivos, "(vb|hmc)"),
    percentage = as.integer(str_extract(archivos, "(95|50|10|25|5)"))
  )

Muestras <- file_info %>%
  filter(type=='predict') %>%
  filter(percentage%in%c(5, 25)) %>% 
  arrange(percentage, desc(intercept), desc(method))

Tiempos <- file_info %>%
  filter(type=='muestra') %>%
  filter(percentage%in%c(5, 25)) %>% 
  arrange(percentage, desc(intercept), desc(method)) %>%
  select(archivos)

i=1
for(i in 1:(nrow(Muestras)/2)){
  primero <- 2*i - 1
  segundo <- 2*i
  int <- as.character(Muestras[primero, 3]) # varios, unico, no: interceptos
  per <- as.character(Muestras[primero, 5]) # porcentaje de muestreo: 5% o 25%
  mod <- str_split_1(folder, pattern='-sim')[1]
  vb <- read.csv(paste0('./', Muestras[primero, 1]))
  hmc <- read.csv(paste0('./', Muestras[segundo, 1]))
  L <- nrow(vb)
  muestra_vb  <- read.csv(paste0('./', Tiempos[primero, 1]))
  muestra_hmc <- read.csv(paste0('./', Tiempos[segundo, 1]))
  if(mod=='log-normal-sesgado'){
    met1 <- metrics(vb$yhat_SN, vb[, 3])
    met2 <- metrics(hmc$yhat_SN, hmc[, 3]) # deberia llamarse yhat_SN
    nombres <- c('Corr.', 'MAE', 'RMSE', 'MAPE', 'Tiempo (s)')
  }
  if(mod=='probit-sesgado'){
    met1 <- metrics_cat(vb$yhat_SN, vb[, 3])
    met2 <- metrics_cat(hmc$yhat_SN, hmc[, 3]) # deberia llamarse yhat_SN
    nombres <- c('Exactitud', '\\makecell[l]{Verdaderos\\\\ positivos}', '\\makecell[l]{Verdaderos\\\\ negativos}', 'F1-Score', 'Tiempo (s)')
  }
  if(mod=='probit-ordenado-sesgado'){
    met1 <- metrics_cat(vb$yhat_SN, vb[, 3], binary = F)
    met2 <- metrics_cat(hmc$yhat_SN, hmc[, 3], binary = F) # deberia llamarse yhat_SN
    delta1_vb <- muestra_vb[nrow(muestra_hmc)-1L, 5]
    delta1_hmc <- muestra_hmc[nrow(muestra_hmc)-1L, 5]
    nombres <- c('Exactitud', '$\\tau$ de Kendall', 'Tiempo (s)')
  }
  s1 <- as.numeric(muestra_vb[nrow(muestra_vb), 2])
  s2 <- as.numeric(muestra_hmc[nrow(muestra_hmc), 2])
  minutos <- c(s1, s2)/60
  segundos <- c(s1, s2)
  # metricas <- data.frame(nombres, formatC(rbind(cbind(met1, met2), minutos), format='f', digits=4))
  metricas <- data.frame(nombres, formatC(rbind(cbind(met1, met2), segundos), format='f', digits=4))
  colnames(metricas) <- NULL; rownames(metricas) <- NULL
  metricas[, 2] <- as.numeric(metricas[, 2])
  metricas[, 3] <- as.numeric(metricas[, 3])
  tmp <- metricas
  if(mod == 'log-normal-sesgado'){
    tmp[nrow(tmp)-1, 2:3] <- -unlist(tmp[nrow(tmp)-1, 2:3]) # mape, mas pequeño: mejor
  }
  tmp[nrow(tmp), 2:3] <- -unlist(tmp[nrow(tmp), 2:3]) # Tiempo (m) # mas pequeño: mejor
  # identificar quien tiene mejores métricas: vb o hmc
  closer <- apply(tmp[, c(2, 3)], MARGIN=1, which.max)
  for(j in 1:length(closer)){
    if(closer[j] == 1L){ # VB is better
      metricas[j, 2] <- paste0('\\bfseries ', metricas[j, 2])
    } else{ # HMC ib better
      metricas[j, 3] <- paste0('\\bfseries ', metricas[j, 3])
    }
  }
  
  metricas <- kable(metricas, format = "latex", escape = FALSE, align = rep("r", 3))
  # Establecer nombres ====
  ruta0 <- r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi)'
  nombre <- paste0(mod, '-muestreo-', per, '-interceptos-', int, '.txt')
  nombre <- file.path(ruta0, 'Metricas', nombre, fsep = '\\')
  writeLines(metricas, nombre)
  Last <- length(readLines(nombre))
  writeLines(readLines(nombre)[-c(1, 2, 3, Last)], nombre)  # skip: 1 blank line, 2 \begin{tabular}, 3 \hline and last: \end{tabular}
  cat('tabla', i, 'de', nrow(Muestras)/2, '\n -- \n')

# Seccion III: Gráfica ====
tmp <- tibble(
Observados=c(vb[, 3], hmc[, 3]), # son identicos
Ajustados=c(vb$yhat_SN, hmc$yhat_SN),
Tipo=c(rep('VB', times=length(vb$X)), rep('HMC', times=length(vb$X))))

if(mod=='log-normal-sesgado'){
  tmp <- tmp %>% mutate(Tipo=factor(Tipo, levels = c('VB', 'HMC'))) # first vb
  g0 <- ggplot(tmp, aes(x=Observados, y=Ajustados, fill=Tipo)) + 
    geom_point(shape=21, color='white', size=3) + 
    geom_abline(slope=1, intercept = 0, linetype=2, col='gray0', linewidth=0.8) + 
    facet_wrap(~Tipo) +
    theme_minimal() + 
    theme(legend.position = 'top', text=element_text(family='serif', size=18)) + 
    scale_fill_manual(values=c('VB'='gray25', 'HMC'='gray75')) + 
    labs(fill='Método', x = expression(Valores ~~ y[ij] ~~ observados),
         y=expression(Valores ~~ y[ij] ~~ ajustados))
  x11(); print(g0)
  nombre <- paste0(mod, '-muestreo-', per, '-interceptos-', int, '.pdf')
  nombre <- file.path(ruta0, 'Scatter', nombre, fsep = '\\')
  ggsave(plot=g0, filename=nombre, units='cm', width=30, height = 15)
} else{
  if(mod == 'probit-ordenado-sesgado'){
    niveles <- c(2, 1, 0)
  } else{
    niveles <- c(1, 0)
  }
  tmp1 <- subset(tmp, Tipo=='VB')
  cm1 <- caret::confusionMatrix(factor(tmp1$Observados, levels=niveles),
                                factor(tmp1$Ajustados, levels=niveles))
  tmp2 <- subset(tmp, Tipo=='HMC')
  cm2 <- caret::confusionMatrix(factor(tmp2$Observados, levels=niveles),
                                factor(tmp2$Ajustados, levels=niveles))
  cm <- rbind(as_tibble(cm1$table), as_tibble(cm2$table))
  
  if(mod=='probit-ordenado-sesgado'){
    cm <- cm %>%
      mutate(Metodo = factor(rep(c('VB', 'HMC'), each=9), levels=c('VB', 'HMC')),
            Tipo=rep(c('Verdaderos (2)', '', '',
                       '', 'Verdaderos (1)', '',
                       '', '', 'Verdaderos (0)'), 2))
  } else{
    cm <- cm %>%
      mutate(Metodo = factor(rep(c('VB', 'HMC'), each=4), levels=c('VB', 'HMC')),
            Tipo=rep(c('Verdaderos\npositivos', 'Falsos\nnegativos',
                       'Falsos\npostitivos', 'Verdaderos\nnegativos'), 2))
  }

  g0 <- ggplot(cm, aes(Prediction, Reference, fill = n)) +
    geom_tile(color='gray85', lwd=1.2) +
    scale_x_discrete(limits = rev)+
    geom_text(aes(label = paste0(Tipo, '\n', n), family = 'serif')) +
    scale_fill_gradient(low = "white", high = "gray75") +
    facet_wrap(~Metodo)+
    theme_minimal() + 
    theme(text=element_text(family='serif', size=18), legend.position = 'none') + 
    labs(x=expression(Valores ~~ y[ij] ~~ ajustados),
         y=expression(Valores ~~ y[ij] ~~ observados),
         fill=expression(n)) + 
    guides(fill = guide_colourbar(barwidth = 0.5,
                                  labels=seq(0, 3000, by=500),
                                  barheight = 20))
  x11(); print(g0)
  nombre <- paste0(mod, '-muestreo-', per, '-interceptos-', int, '.pdf')
  nombre <- file.path(ruta0, 'Scatter', nombre, fsep = '\\')
  ggsave(plot=g0, filename=nombre, units='cm', width=30, height = 15)
}

}

# Seccion IV: SSVS ==== 
#¿lo hare para todos?

archivos <- dir()
archivos <- archivos[startsWith(archivos, "SSVS")]

file_info <- tibble(tibble(archivos)) %>%
  mutate(
    type = str_extract(archivos, "^(muestra|predict|SSVS)"),
    intercept = str_extract(archivos, "(varios|unico|no)"),
    method = str_extract(archivos, "(vb|hmc)"),
    percentage = as.integer(str_extract(archivos, "(95|50|10|25|5)"))
  )

Muestras <- file_info %>%
  filter(type=='SSVS') %>%
  filter(percentage%in%c(5, 25)) %>% 
  arrange(percentage, desc(intercept), desc(method))

i=1
for(i in 1:(nrow(Muestras)/2)){
  primero <- 2*i - 1
  segundo <- 2*i
  int <- as.character(Muestras[primero, 3]) # varios, unico, no: interceptos
  per <- as.character(Muestras[primero, 5]) # porcentaje de muestreo: 5% o 25%
  mod <- str_split_1(folder, pattern='-sim')[1]
  vb <- read.csv(paste0('./', Muestras[primero, 1]))
  hmc <- read.csv(paste0('./', Muestras[segundo, 1]))
  L <- nrow(vb)
  
vb <- vb %>% mutate(Metodo = 'VB', X = NULL)
hmc <- hmc %>%mutate(Metodo = 'HMC', X = NULL)
tmp <- rbind(vb, hmc)
colnames(tmp) <- c('Prob', 'Config', 'Método')

tmp <- tmp %>%
  mutate(Prob = factor(Prob), Config = str_split(Config, "")) %>%
  unnest_wider(Config, names_sep = "") %>%
  rename_with(~ paste0("beta_", seq_along(.)), starts_with("Config")) %>%
  pivot_longer(starts_with("beta_"), names_to = "beta", values_to = "included") %>%
  mutate(included = as.integer(included)) %>%
  complete(Prob, beta, fill = list(included=NA)) %>%
  mutate(Método=factor(Método, levels=c('VB', 'HMC')))

g0 <- ggplot(tmp, aes(x = beta, y = Prob, fill = factor(included))) +
  geom_tile(color = "gray50", linewidth = 0.4) + 
  facet_wrap(~Método, scales = 'free_y') + 
  scale_x_discrete(labels=parse(text = paste0("beta[", 1:5, "]"))) + 
  theme_minimal() + 
  labs(fill='Inclusión', y='Prob. a posteriori', x='') + 
  scale_fill_discrete(labels=c('0'='No', '1'='Sí')) +
  scale_fill_manual(values=c('1'='gray75', '0'='gray25'))+
  theme(legend.position = 'top', text=element_text(family='serif'),
        plot.margin = margin(0, 0, 0, 0),  axis.ticks.length = unit(0, "pt"))
  x11(); print(g0)
  nombre <- paste0(mod, '-muestreo-', per, '-interceptos-', int, '.pdf')
  nombre <- file.path(ruta0, 'SSVS', nombre, fsep = '\\')
  ggsave(plot=g0, filename=nombre, units='cm', width=20, height = 6)

}

}

# ○○○○

# Parte II: identificabilidad de rho ====
# aqui interesan más las estimaciones de rho1 a rho4
# library(bayesplot)
setwd(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Codigos\Simulacion-Modelos)')
master <- getwd()
carpetas <- dir()[endsWith(dir(), 'sim2')]

# navegar por las carpetas ====
folder <- carpetas[1]

for(folder in carpetas){
  
  setwd(file.path(master, folder, fsep = '/'))
  
  cat('Estamos en la carpeta:', folder, '\n')
  
  ## Sección I: muestras de rho1 a rho 4 ====
  
  archivos <- dir()
  archivos <- archivos[startsWith(archivos, 'muestra') & endsWith(archivos, '.csv')]
  
  file_info <- tibble(tibble(archivos)) %>%
    mutate(
      type = str_extract(archivos, "^(muestra|predict|SSVS)"),
      intercept = str_extract(archivos, "(varios|unico|no)"),
      method = str_extract(archivos, "(vb|hmc)"),
      percentage = as.integer(str_extract(archivos, "(95|50|25|5)"))
    )
  
  Muestras <- file_info %>%
    filter(type=='muestra') %>%
    filter(percentage%in%c(5, 25)) %>% 
    arrange(percentage, desc(intercept), desc(method))
  
  i=1
  tiempos <- c()
  for(i in 1:(nrow(Muestras)/2)){
    primero <- 2*i - 1
    segundo <- 2*i
    int <- as.character(Muestras[primero, 3]) # varios, unico, no: interceptos
    per <- as.character(Muestras[primero, 5]) # porcentaje de muestreo: 5% o 25%
    mod <- str_split_1(folder, pattern='-sim')[1]
    vb <- read.csv(paste0('./', Muestras[primero, 1]))
    hmc <- read.csv(paste0('./', Muestras[segundo, 1]))
    L <- nrow(vb)
    tiempos <- c(tiempos, vb[L, 2], hmc[L, 2])
    # Etiquetas para los parametros
    parametros <- c(paste0('$\\rho_{', 1:4,'}$'))
    names(parametros) <- c(paste0('rho', 1:4))
    # El formato se lo daremos en LaTeX, es decir, el redondeo/número de decimales
    valores_reales <- vb[1:4, 2]
    vb_mp <- vb[1:4, 3] # media posterior
    vb_025 <- vb[1:4, 5]
    vb_975 <- vb[1:4, 6]
    hmc_mp <- hmc[1:4, 3]
    hmc_025 <- hmc[1:4, 5]
    hmc_975 <- hmc[1:4, 6]
    tmp <- data.frame(parametros, valores_reales,
                      vb_mp, vb_025, vb_975, hmc_mp, hmc_025, hmc_975)
    # identificar quien se acercó más: vb o hmc
    closer <- apply(tmp[, c(2, 3, 6)], MARGIN=1, which_closer)
    for(j in 1:length(closer)){
      if(closer[j] == 1L){
        tmp[j, 3] <- paste0('\\bfseries ', tmp[j, 3])
      } else{
        tmp[j, 6] <- paste0('\\bfseries ', tmp[j, 6])
      }
    }
    
    if(int == 'unico'){
      # mu <- mean(as.numeric(str_split_1(tmp['mu', 2],pattern = '\\|')))
      # tmp['mu', 2] <- mu
    }
    
    #colnames(tmp) <- c('Parámetros', 'Valor real', 'Media posterior', 'I.C. 2.5\\%', 'I.C. 97.5\\%',
    #                      'Media posterior', 'I.C. 2.5\\%', 'I.C. 97.5\\%')
    tmp <- as.matrix(tmp); colnames(tmp) <- NULL; rownames(tmp) <- NULL
    tmp <- kable(tmp, format = "latex", escape = FALSE, align = rep("r", 8)) # , digits = 4, pero lo haremos en LaTeX
    # Establecer nombres ====
    ruta0 <- r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi)'
    nombre <- paste0(mod, '-muestreo-', per, '-interceptos-', int, '-rho11.txt')
    nombre <- file.path(ruta0, 'Tablas', nombre, fsep = '\\')
    writeLines(tmp, nombre)
    Last <- length(readLines(nombre))
    writeLines(readLines(nombre)[-c(1, 2, 3, Last)], nombre)  # skip: 1 blank line, 2 \begin{tabular}, 3 \hline and last: \end{tabular}
    cat('tabla', i, 'de', nrow(Muestras)/2, '\n -- \n')
    
    ## Sección II: histogramas de rho1 a rho 4====
    archivos <- dir()
    archivos <- archivos[endsWith(archivos, '.RData')]
    
    file_info <- tibble(tibble(archivos)) %>%
      mutate(
        type = str_extract(archivos, "^(muestra|predict|SSVS)"),
        intercept = str_extract(archivos, "(varios|unico|no)"),
        method = str_extract(archivos, "(VB|HMC)"),
        percentage = as.integer(str_extract(archivos, "(95|50|25|5)"))
      )
    
    Muestras <- file_info %>%
      filter(type=='muestra') %>%
      filter(percentage%in%c(5, 25)) %>% 
      arrange(percentage, desc(intercept), desc(method))
    
    i=1
    for(i in 1:(nrow(Muestras)/2)){
      primero <- 2*i - 1
      segundo <- 2*i
      load(paste0(Muestras[primero, 1]))
      load(paste0(Muestras[segundo, 1]))
      rho1 <- as_tibble(muestra_VB$draws('rho')); colnames(rho1)
      rho2 <- as_tibble(muestra_HMC$draws('rho')); colnames(rho2) <- str_remove(colnames(rho2),pattern = '1.')
      rho <- tibble(rbind(rho1, rho2), Metodo=factor(rep(c('VB', 'HMC'), each=nrow(rho1)), levels=c('VB', 'HMC')))
      rho <- rho %>% pivot_longer(cols=starts_with('rho'), names_to='Grupo', values_to = 'rho')
      g0 <- ggplot(rho, aes(x=rho, y=after_stat(density))) + 
        geom_histogram(bins = 30, fill = 'gray75', color='white',) + 
        facet_wrap(Metodo~Grupo, nrow = 2, ncol=4 ,scales = "free", labeller=label_parsed) +
        # facet_grid(Metodo~Grupo, scales = "free", labeller=label_parsed) +
        theme_bw() + theme(text=element_text(family='serif')) + 
        geom_hline(yintercept=0) + labs(x=expression(rho), y='Densidad')
      # x11(); print(g0)
      nombre <- paste0(mod, '-muestreo-', per, '-interceptos-', int, '-rho11.pdf')
      nombre <- file.path(ruta0, 'Scatter', nombre, fsep = '\\')
      ggsave(plot=g0, filename=nombre, units='cm', width=30, height = 15)
      
    }
  }
}

