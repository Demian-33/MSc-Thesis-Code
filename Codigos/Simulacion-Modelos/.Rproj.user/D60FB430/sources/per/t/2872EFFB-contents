
# =============== Analisis a nivel hogares del ictpc =================

# Este código proceso los datos del ICTPC

rm(list=ls())
gc()

# (0) Cargar librerías ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(cmdstanr)
library(kableExtra)
library(stringr)
library(sf) # st_read
library("ggspatial") # annotation_map_tile
# library(showtext)
# showtext_auto()
# showtext_opts(dpi = 300)
# rebuild_cmdstan()

out_stan <- posterior::default_summary_measures()
myquant <- function(a){
  b <- quantile(a, c(0.025, 0.975))
  return(b)
}
out_stan[5] <- 'myquant'

# (1) Datos para el modelo ====
getwd()
load('./Xyregion-CDMX.Rdata')
# source('./funciones-final.R')
source('./funciones-final_.R')

muni <- c(
  'Azcapotzalco', 'Coyoacán', 'Cuajimalpa de Morelos', 'Gustavo A. Madero',
  'Iztacalco', 'Iztapalapa', 'La Magdalena Contreras', 'Milpa Alta',
  'Álvaro Obregón', 'Tláhuac', 'Tlalpan', 'Xochimilco', 
  'Benito Juárez','Cuauhtémoc', 'Miguel Hidalgo', 'Venustiano Carranza') %>% factor

draw_key_cust <- function(data, params, size) {
  data_text <- data
  data_text[c("family")] <- 'serif'
  data_text[c("label")] <- key_label2[names(pal)[match(data$colour, pal)]]
  # data_text$label <- key_label[names(pal)[match(data$colour, pal)]] # similar a arriba
  data_text[c("fill")] <- NULL # ¿no se que hace?
  data_text$colour <- "black" # color dentro del cuadrito
  data_text$alpha <- 1
  data_text$size <- 11 / .pt
  
  grid::grobTree(
    draw_key_rect(data, list()),
    draw_key_text(data_text, list())
  )
}

carto_lightnl_url <- "https://a.basemaps.cartocdn.com/light_nolabels/${z}/${x}/${y}.png" # no labels
carto_light_url <- "https://a.basemaps.cartocdn.com/light_all/${z}/${x}/${y}.png"
shp <- st_read(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX\poligonos_alcaldias_cdmx\poligonos_alcaldias_cdmx.shp)')
shp$CVE_MUN <- str_pad(as.numeric(shp$CVE_MUN)-1, width=3, pad='0')
pal <- scales::grey_pal(start=0, end=0.11)(16)
names(pal) <- shp$NOMGEO
key_label <- shp$NOMGEO
names(key_label) <- shp$CVE_MUN
key_label2 <- shp$CVE_MUN
names(key_label2) <- shp$NOMGEO

LPI_urb <- 4564.97
LPI_rur <- 3296.92
LPEI_urb <- 2354.65
LPEI_rur <- 1800.55

## (1.1) Definir areas pequeñas ====
group <- as.numeric(region)
K <- length(unique(group))

## (1.2) Covariables categoricas (no ordenadas) a matriz diseño ====
for(i in 1:ncol(Xcat)){
  Xcat[, i] <- as.factor(as.integer(Xcat[, i]))
}
tail(tibble(Xcat)) # ya todas son factor
Xcat <- model.matrix( ~ 0 + . , data=Xcat) # generar a matriz diseño

## (1.3) Remover combinaciones lineales y correlaciones ====
Xtmp <- cbind(Xcon, Xbin)
a <- caret::findLinearCombos(Xtmp)
Xtmp <- Xtmp[, -c(a$remove)]
b <- caret::findCorrelation(cov(Xtmp), verbose=T, exact=T, cutoff = 0.8)
Xtmp <- Xtmp[, -c(b)]

dim(Xtmp)

## (1.4) Componentes principales ====
# Eigen <- eigen(cov(Xcon))
# cumsum(Eigen$values) / ncol(Xcon) # 19: 90%, 26: 95%
# Xpc <- Xcon %*% Eigen$vectors[, 1:26]
# a <- caret::findLinearCombos(Xbin)$remove
# Xtmp <- cbind(Xpc, Xbin[, -a])

# (2) Modelo log-normal ====
# 80-20: se calcula un grafico y metricas
# pronostico: se calcula grafico y metricas
# estimacion (todos los datos): cantidades posteriores
set.seed(1)
# T: estimacion, F: 80-20 (establecer prop=0.2) y pronostico
datos <- get_obs(censo=T)
datos <- get_obs(censo=F, prop=0.2)

## (2.1) Cargar el modelo con VB ====
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-log-normal-sesgado-estimacion.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-log-normal-sesgado-8020.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-log-normal-sesgado-pronostico0.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-log-normal-sesgado-pronostico1.RData)')
muestra_vb <- muestra_VB$summary(variables = c('rho', 'sigma', 'mu', 'b'), out_stan)

## (2.2) Cargar el modelo con HMC ====
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-log-normal-sesgado-estimacion.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-log-normal-sesgado-8020.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-log-normal-sesgado-pronostico0.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-log-normal-sesgado-pronostico1.RData)'); muestra_HMC <- fit_hmc; rm(fit_hmc)
muestra_hmc <- muestra_HMC$summary(variables = c('rho', 'sigma', 'mu', 'b'), c(out_stan, 'rhat'))

## (2.3) Estimaciones de los parametros (solo estimación y pronostico) ====
parametros <- c(
  paste0('$\\rho_{', 1:K, '}$'),
  paste0('$\\sigma^{2}$'),
  # paste0('$\\mu_{', 1:K, '}$')
  paste0('$\\mu$')
)

colnames(muestra_hmc) <- paste0(colnames(muestra_hmc), 2)

pred_out <- tibble(parametros,
                   muestra_vb[1:length(parametros), c(2, 6, 7)],
                   muestra_hmc[1:length(parametros), c(2, 6, 7, 8)])

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-tabla.txt')

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-pronostico0-tabla.txt')
pred_out <- pred_out[c(8, 15, 17, 25, 32), ]

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-pronostico1-tabla.txt')
pred_out <- pred_out[c(8, 15, 17, 18), ]

writeLines(pred_out %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

## (2.4) Ajustados vs Observados (solo 80-20 y pronostico) ====
out <- datos$idx # 80-20 y pronostico
# out <- datos$master # estimacion
tproc_vb <- readLines(con=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-tiempos-pronostico1.txt)') %>%
  str_split_1(pattern='~') %>% as.numeric()
tproc_hmc <- readLines(con=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-tiempos-pronostico1.txt)') %>%
  str_split_1(pattern='~') %>% as.numeric()

# ¿metricas en escala real?, sirve más para comparación
yhat_vb <- apply(exp(muestra_VB$draws('y_all')), 2, mean)[out] # ajustados u 80-20
yhat_vb <- apply(exp(muestra_VB$draws('y_mis')), 2, mean) # pronostico
yhat_hmc <- apply(exp(drop(muestra_HMC$draws('y_all'))), 2, mean)[out] # ajustados u 80-20
yhat_hmc <- apply(exp(drop(muestra_HMC$draws('y_mis'))), 2, mean) # pronostico
tmp <- cbind(metrics(yhat_vb, exp(datos$y_test)), metrics(yhat_hmc, exp(datos$y_test)))
tmp <- rbind(tmp, c(tproc_vb[1], tproc_hmc[1]))
rownames(tmp)[5] <- 'Tiempo (s)'
round(tmp, 5)
met <- apply(tmp, 1, which.min)
met[1] <- ifelse(met[1]==1, 2, 1)
for(k in 1:nrow(tmp)){
  tmp[k, met[k]] <- paste0('\\bfseries ', tmp[k, met[k]])
}

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-metricas-pronostico1.txt')
writeLines(tmp %>% kable(format = 'latex', digits = 4, escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:3, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

# grafico en escala log: se ve mejor
yhat_vb <- apply(muestra_VB$draws('y_all'), 2, mean)[out] # ajustados u 80-20
yhat_vb <- apply(muestra_VB$draws('y_mis'), 2, mean) # pronostico
yhat_hmc <- apply(drop(muestra_HMC$draws('y_all')), 2, mean)[out] # ajustados u 80-20
yhat_hmc <- apply(drop(muestra_HMC$draws('y_mis')), 2, mean) # pronostico

tmp <- tibble(Observados=datos$y_test,
              AjustadosVB=yhat_vb, AjustadosHMC=yhat_hmc)
# agregar alcaldia en pronosticos:
tmp <- tmp %>% mutate(Alcaldia=factor(muni[datos$group[datos$idx]]))

tmp <- tmp %>%
  pivot_longer(cols=starts_with('Ajustados'),
    names_to = 'Metodo',
    names_prefix = 'Ajustados',
    values_to = 'Ajustados') %>%
  mutate(Metodo=factor(Metodo, levels=c('VB', 'HMC')))

# para 8020

g0 <- ggplot(tmp%>%filter(Metodo=='VB'), aes(x=Observados, y=Ajustados)) +
  geom_point(color='gray50') +
  geom_abline(slope=1, intercept = 0, color='gray0', linetype=2, , linewidth=0.8)+
  theme_minimal() +
  theme(text=element_text(family='serif', size=20), legend.position = 'top') + 
  guides(shape = "none") +
  labs(x='log-ICTPC observado', y = 'log-ICTPC ajustado')
x11(); print(g0)

g1 <- ggplot(tmp%>%filter(Metodo=='HMC'), aes(x=Observados, y=Ajustados)) +
  geom_point(color='gray40') +
  geom_abline(slope=1, intercept = 0, color='gray0', linetype=2, linewidth=0.8)+
  theme_minimal() +
  theme(text=element_text(family='serif', size=20), legend.position = 'top') + 
  guides(shape = "none") +
  labs(x='log-ICTPC observado', y = 'log-ICTPC ajustado')
x11(); print(g1)

# para pronostico

g0 <- ggplot(tmp%>%filter(Metodo=='VB'), aes(x=Observados, y=Ajustados, color=Alcaldia, shape=Alcaldia, group=Alcaldia)) + 
  geom_abline(intercept = 0, slope=1, linewidth=1.2, color='gray5', lty=3) +
  geom_point(size=2) +
  scale_color_manual(values=c('Milpa Alta'='gray25', 'Miguel Hidalgo'='gray75'))+
  geom_smooth(method = 'lm', formula=y~x, se=F, color='gray25')+
  theme_minimal() + 
  theme(legend.position = 'top', text=element_text(family='serif', size=20)) + 
  labs(x='log-ICTPC ajustado', y='log-ICTPC observado')
x11(); print(g0)

g1 <- ggplot(tmp%>%filter(Metodo=='HMC'), aes(x=Observados, y=Ajustados, color=Alcaldia, shape=Alcaldia, group=Alcaldia)) + 
  geom_abline(intercept = 0, slope=1, linewidth=1.2, color='gray5', lty=3) +
  geom_point(size=2) +
  scale_color_manual(values=c('Milpa Alta'='gray25', 'Miguel Hidalgo'='gray75'))+
  geom_smooth(method = 'lm', formula=y~x, se=F, color='gray25')+
  theme_minimal() + 
  theme(legend.position = 'top', text=element_text(family='serif', size=20)) + 
  labs(x='log-ICTPC ajustado', y='log-ICTPC observado')
x11(); print(g1)

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-scatter8020')

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-scatterpronostico')

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-scatterpronostico1')

ggsave(plot=g0, filename=paste0(nombre, 'VB.pdf'), width=20, height=20, units='cm')
ggsave(plot=g1, filename=paste0(nombre, 'HMC.pdf'), width=20, height=20, units='cm')


## (2.7) Porcentaje de personas ====
## se promedian las muestras a posteriori exponenciadas
out <- which(is.na(ym))
fitted_vb <- apply(exp(muestra_VB$draws('y_all')), 2, mean)
fitted_hmc <- apply(exp(drop(muestra_HMC$draws('y_all'))), 2, mean)
dftmp <- tibble(AjustadosVB=fitted_vb[out],
                AjustadosHMC=fitted_hmc[out],
                Ámbito=Xide$rururb[out],
                Factor=Xide$factor[out],
                Tamaño=Xide$tam_hog[out],
                Alcaldia=muni[region][out])
# tipo log normal
dftmp <- dftmp %>%
  mutate(LPI_vb = case_when(
    `Ámbito` == 0 & AjustadosVB < LPI_urb ~ 1L,
    `Ámbito` == 1 & AjustadosVB < LPI_rur ~ 1L,
    TRUE                                ~ 0L),
    LPI_hmc = case_when(
      `Ámbito` == 0 & AjustadosHMC < LPI_urb ~ 1L,
      `Ámbito` == 1 & AjustadosHMC < LPI_rur ~ 1L,
      TRUE                                ~ 0L)) %>%
  mutate(LPEI_vb = case_when(
    `Ámbito` == 0 & AjustadosVB < LPEI_urb ~ 1L,
    `Ámbito` == 1 & AjustadosVB < LPEI_rur ~ 1L,
    TRUE                                ~ 0L),
    LPEI_hmc = case_when(
      `Ámbito` == 0 & AjustadosHMC < LPEI_urb ~ 1L,
      `Ámbito` == 1 & AjustadosHMC < LPEI_rur ~ 1L,
      TRUE                                ~ 0L))
#
dftmp2 <- dftmp %>%
  group_by(Alcaldia) %>%
  summarise(Total = sum(Factor*Tamaño),
            Total_LPI_vb = sum(Factor*LPI_vb*Tamaño), Total_LPEI_vb = sum(Factor*LPEI_vb*Tamaño),
            Total_LPI_hmc = sum(Factor*LPI_hmc*Tamaño), Total_LPEI_hmc = sum(Factor*LPEI_hmc*Tamaño),
            Porcentaje_LPI_vb = 100*Total_LPI_vb/Total, Porcentaje_LPEI_vb = 100*Total_LPEI_vb/Total,
            Porcentaje_LPI_hmc = 100*Total_LPI_hmc/Total, Porcentaje_LPEI_hmc = 100*Total_LPEI_hmc/Total)

dftmp2 %>% select(Total, starts_with('Porcentaje'))
sum(dftmp2$Total_LPI_vb) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPI_hmc) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPEI_vb) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPEI_hmc) / sum(dftmp2$Total) # 25.7 % 

# ordenar de acuerdo a muni
dftmp2 <- dftmp2[match(muni, dftmp2$Alcaldia), ]
personas <- dftmp2 %>% select(Alcaldia, starts_with('Porcentaje') & ends_with('vb'))
personas <- dftmp2 %>% select(Alcaldia, starts_with('Porcentaje'))
nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-personas.txt')
writeLines(personas %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

## (2.8) Mapas ====
orden <- match(dftmp2$Alcaldia, shp$NOMGEO) # el elemento i (i=1, 2, ...) de #1 esta en la posición indicada de #2
orden2 <- match(muni, shp$NOMGEO) # el elemento i (i=1, 2, ...) de #1 esta en la posición indicada de #2
shp$LPI_vb_porcentaje[orden] <- dftmp2$Porcentaje_LPI_vb
shp$LPEI_vb_porcentaje[orden] <- dftmp2$Porcentaje_LPEI_vb 
shp$LPI_hmc_porcentaje[orden] <- dftmp2$Porcentaje_LPI_hmc
shp$LPEI_hmc_porcentaje[orden] <- dftmp2$Porcentaje_LPEI_hmc
pal <- scales::grey_pal(start=0, end=0.11)(16)
names(pal) <- shp$NOMGEO
key_label <- shp$NOMGEO
names(key_label) <- shp$CVE_MUN
key_label2 <- shp$CVE_MUN
names(key_label2) <- shp$NOMGEO

g1 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPI_vb_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, -1, "cm"),
    text = element_text(family = "serif", size=26),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g1)

g2 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPEI_vb_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPEI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, -1, "cm"),
    text = element_text(family = "serif", size=26),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g2)

g3 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPI_hmc_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, -1, "cm"),
    text = element_text(family = "serif", size=26),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g3)

g4 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPEI_hmc_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPEI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, -1, "cm"),
    text = element_text(family = "serif", size=26),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g4)

# 17.53
nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-')
ggsave(plot=g1, filename=paste0(nombre, 'mapavb-lpi.pdf'), width=24, height = 20, units='cm')
ggsave(plot=g2, filename=paste0(nombre, 'mapavb-lpei.pdf'), width=24, height = 20, units='cm')
ggsave(plot=g3, filename=paste0(nombre, 'mapahmc-lpi.pdf'), width=24, height = 20, units='cm')
ggsave(plot=g4, filename=paste0(nombre, 'mapahmc-lpei.pdf'), width=24, height = 20, units='cm')

## (2.9) Sesgo promedio ====
idx <- which(!is.na(ym))
fitted_vb <- apply(exp(muestra_VB$draws('y_all')), 2, mean)
fitted_hmc <- apply(exp(drop(muestra_HMC$draws('y_all'))), 2, mean)
dftmp <- tibble(
  ICTPC_obs = exp(ym[idx]),
  ICTPC_fit_vb=fitted_vb[idx],
  ICTPC_fit_hmc=fitted_hmc[idx],
  Alcaldia=muni[region][idx])

dftmp <- dftmp %>%
  group_by(Alcaldia) %>%
  summarise(mean_ICTPC_obs = mean(ICTPC_obs, na.rm=T), # cambiar por median, esta curioso
            mean_ICTPC_fit_vb = mean(ICTPC_fit_vb, na.rm=T),
            mean_ICTPC_fit_hmc = mean(ICTPC_fit_hmc, na.rm=T),
            diff_ICTPC_vb = abs(mean_ICTPC_obs-mean_ICTPC_fit_vb),
            diff_ICTPC_hmc = abs(mean_ICTPC_obs-mean_ICTPC_fit_hmc))

mean(dftmp$diff_ICTPC_vb); mean(dftmp$diff_ICTPC_hmc)

dftmp %>% kable(format = 'latex', digits=3,
         label='sesgo',
         caption='Elaboración propia basada en la muestra \\underline{a posteriori}.')

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'log-normal-sesgado-sesgo.txt')

writeLines(dftmp %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

## (2.10) SSVS ====
# vb <- cbind(posterior_prob(muestra_VB), Metodo='VB')
# hmc <- cbind(posterior_prob(muestra_HMC), Metodo='HMC')
# tmp <- rbind(vb, hmc)
# colnames(tmp) <- c('Prob', 'Config', 'Método')
# library(stringr)
# tmp <- tmp %>%
#   mutate(Prob = factor(Prob), Config = str_split(Config, "")) %>%
#   unnest_wider(Config, names_sep = "") %>%
#   rename_with(~ paste0("beta_", seq_along(.)), starts_with("Config")) %>%
#   pivot_longer(starts_with("beta_"), names_to = "beta", values_to = "included") %>%
#   mutate(included = as.integer(included)) %>%
#   complete(Prob, beta, fill = list(included=NA)) %>%
#   mutate(Método=factor(Método, levels=c('VB', 'HMC')))
# 
# g0 <- ggplot(tmp, aes(x = beta, y = Prob, fill = factor(included))) +
#   geom_tile(color = "white", linewidth = 0.4) + 
#   facet_wrap(~Método, scales = 'free_y') + 
#   scale_x_discrete(labels=parse(text = paste0("beta[", 1:5, "]"))) + 
#   theme_minimal() + 
#   labs(fill='Inclusión', y='Prob. a posteriori', x='') + 
#   scale_fill_discrete(labels=c('0'='No', '1'='Sí'))+
#   theme(legend.position = 'top', text=element_text(family='serif'),
#         plot.margin = margin(0, 0, 0, 0),  axis.ticks.length = unit(0, "pt"))
# x11(); print(g0)
# nombre <- paste0(mod, '-muestreo-', per, '-interceptos-', int, '.pdf')
# nombre <- file.path(ruta0, 'SSVS', nombre, fsep = '\\')
# ggsave(plot=g0, filename=nombre, units='cm', width=20, height = 6)

# (3) Modelo probit ====
## (3.1) Seleccionar datos ====
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
    TRUE ~ 0) # no carencias
)
table(tmp$ym_bin)
table(tmp$ym_ord)

set.seed(1)
datos <- get_obs(censo=FALSE, prop=0.2, y = tmp$ym_bin)
# no invertir la proporcion de 0's - 1's
# datos$y_obs <- ifelse(datos$y_obs==1, 0, 1)

## (3.2) Cargar el modelo con VB ====
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-probit-sesgado-estimacion.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-probit-sesgado-8020.RData)')
# muestra_vb <- muestra_probitVB$summary(variables = c('rho', 'b'))
muestra_vb <- muestra_probitVB$summary(variables = c('rho', 'mu'), c(out_stan))

## (3.3) Cargar el modelo con HMC ====
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-probit-sesgado-estimacion.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-probit-sesgado-8020.RData)')
muestra_probitHMC$code()
muestra_hmc <- muestra_probitHMC$summary(variables = c('rho', 'mu'), c(out_stan, 'rhat'))
colnames(muestra_hmc) <- paste0(colnames(muestra_hmc), 2)

## (3.4) Tabla de estimaciones (solo estimación) ====
parametros <- c(
  paste0('$\\rho_{', 1:K, '}$'),
  paste0('$\\mu_{', 1:K, '}$')
)

pred_out <- tibble(parametros,
                   muestra_vb[1:length(parametros), c(2, 6, 7)],
                   muestra_hmc[1:length(parametros), c(2, 6, 7, 8)])

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-sesgado-tabla.txt')
writeLines(pred_out %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

## (3.6) Métricas ====
out <- datos$idx
tproc_vb <- readLines(con=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-tiempos-8020.txt)') %>%
  str_split_1(pattern='~') %>% as.numeric()
tproc_hmc <- readLines(con=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-tiempos-8020.txt)') %>%
  str_split_1(pattern='~') %>% as.numeric()
yhat_vb <- findInterval(apply(muestra_probitVB$draws('z_all'), 2, mean), 0)[out]
yhat_hmc <- findInterval(apply(drop(muestra_probitHMC$draws('z_all')), 2, mean), 0)[out]
# no invertir
# yhat_vb <- ifelse(yhat_vb==1, 0, 1)
# yhat_hmc <- ifelse(yhat_hmc==1, 0, 1)
# yhat_hmc <- yhat_vb
tmp <- cbind(metrics_cat(yhat_vb, datos$y_test), metrics_cat(yhat_hmc, datos$y_test))
tmp <- rbind(tmp, c(tproc_vb[2], tproc_hmc[2]))
rownames(tmp)[5] <- 'Tiempo (s)'

met <- apply(tmp, 1, which.max)
met[5] <- ifelse(met[5]==1, 2, 1)
for(k in 1:nrow(tmp)){
  tmp[k, met[k]] <- paste0('\\bfseries ', tmp[k, met[k]])
}

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-sesgado-metricas8020.txt')
writeLines(tmp %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:3, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

tmp <- tibble(Observados=datos$y_test,
              AjustadosVB=yhat_vb, AjustadosHMC=yhat_hmc)
tmp <- tmp %>%
  pivot_longer(cols=starts_with('Ajustados'),
               names_to = 'Metodo',
               names_prefix = 'Ajustados',
               values_to = 'Ajustados') %>%
  mutate(Metodo=factor(Metodo, levels=c('VB', 'HMC')))

tmp1 <- subset(tmp, Metodo=='VB')
cm1 <- caret::confusionMatrix(factor(tmp1$Observados, levels=c(1, 0)),
                              factor(tmp1$Ajustados, levels=c(1, 0)))
tmp2 <- subset(tmp, Metodo=='HMC')
cm2 <- caret::confusionMatrix(factor(tmp2$Observados, levels=c(1, 0)),
                              factor(tmp2$Ajustados, levels=c(1, 0)))

cm <- rbind(as_tibble(cm1$table), as_tibble(cm2$table))
cm <- cm %>%
    mutate(Metodo = factor(rep(c('VB', 'HMC'), each=4), levels=c('VB', 'HMC')),
    Tipo=rep(c('Verdaderos\npositivos', 'Falsos\nnegativos',
               'Falsos\npostitivos', 'Verdaderos\nnegativos'), 2))

g0 <- ggplot(cm%>%filter(Metodo=='VB'), aes(Prediction, Reference, fill = n)) +
  geom_tile(color='gray85', lwd=1.2) +
  scale_x_discrete(limits = rev)+
  geom_text(aes(label = paste0(Tipo, '\n', n), family = 'serif')) +
  scale_fill_gradient(low = "white", high = "gray50",breaks=seq(0, 400, by=100)) +
  facet_wrap(~Metodo)+
  theme_minimal() +
  theme(text=element_text(family='serif', size=20), legend.position = 'none') +
  labs(x=expression(Valores ~~ y[ij] ~~ ajustados),
       y=expression(Valores ~~ y[ij] ~~ observados),
       fill=expression(n))
  # guides(fill = guide_colourbar(barheight=unit(12, 'cm'), barwidth=unit(1, 'cm'))) 
x11(); print(g0)

g1 <- ggplot(cm%>%filter(Metodo=='HMC'), aes(Prediction, Reference, fill = n)) +
  geom_tile(color='gray85', lwd=1.2) +
  scale_x_discrete(limits = rev)+
  geom_text(aes(label = paste0(Tipo, '\n', n), family = 'serif')) +
  scale_fill_gradient(low = "white", high = "gray50",breaks=seq(0, 400, by=100)) +
  facet_wrap(~Metodo)+
  theme_minimal() +
  theme(text=element_text(family='serif', size=20), legend.position = 'none') +
  labs(x=expression(Valores ~~ y[ij] ~~ ajustados),
       y=expression(Valores ~~ y[ij] ~~ observados),
       fill=expression(n)) +
  guides(fill = guide_colourbar(barheight=unit(12, 'cm'), barwidth=unit(1, 'cm'))) 
x11(); print(g1)

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-sesgado-scatter8020')

ggsave(plot=g0, filename=paste0(nombre, 'VB.pdf'), width=20, height=20, units='cm')
ggsave(plot=g1, filename=paste0(nombre, 'HMC.pdf'), width=20, height=20, units='cm')

## (3.7) Porcentaje de personas ====
out <- which(is.na(ym))
fitted_vb <- findInterval(apply(muestra_probitVB$draws('z_all'), 2, mean), 0)
fitted_hmc <- findInterval(apply(drop(muestra_probitHMC$draws('z_all')), 2, mean), 0)
# no invertir
# fitted_vb <- ifelse(fitted_vb==1, 0, 1)
# fitted_hmc <- ifelse(fitted_hmc==1, 0, 1)
# fitted_hmc <- fitted_vb
dftmp <- tibble(LPI_vb=fitted_vb[out],
                LPI_hmc=fitted_hmc[out],
                Ámbito=Xide$rururb[out],
                Factor=Xide$factor[out],
                Tamaño=Xide$tam_hog[out],
                Alcaldia=muni[region][out])
#
dftmp2 <- dftmp %>%
  group_by(Alcaldia) %>%
  summarise(Total = sum(Factor*Tamaño),
            Total_LPI_vb = sum(Factor*LPI_vb*Tamaño),
            Total_LPI_hmc = sum(Factor*LPI_hmc*Tamaño),
            Porcentaje_LPI_vb = 100*Total_LPI_vb/Total,
            Porcentaje_LPI_hmc = 100*Total_LPI_hmc/Total)

dftmp2 %>% select(Total, starts_with('Porcentaje'))
sum(dftmp2$Total_LPI_vb) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPI_hmc) / sum(dftmp2$Total) # 25.7 % 

personas <- dftmp2 %>% select(Alcaldia, starts_with('Porcentaje') & ends_with('vb'))
personas <- dftmp2 %>% select(Alcaldia, starts_with('Porcentaje'))
# ordenar de acuerdo a muni
personas <- personas[match(muni, personas$Alcaldia), ]

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-sesgado-personas.txt')
writeLines(personas %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

## (3.8) Mapas ====
orden <- match(dftmp2$Alcaldia, shp$NOMGEO) # el elemento i (i=1, 2, ...) de #1 esta en la posición indicada de #2
orden2 <- match(muni, shp$NOMGEO) # el elemento i (i=1, 2, ...) de #1 esta en la posición indicada de #2
shp$rho_vb_estimado[orden2] <- apply(muestra_probitVB$draws('rho'), 2, mean)
shp$LPI_vb_porcentaje[orden] <- dftmp2$Porcentaje_LPI_vb

g1 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPI_vb_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    text = element_text(family = "serif", size=18),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g1)

g2 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPEI_vb_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPEI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    text = element_text(family = "serif", size=18),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g2)

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-sesgado-')
ggsave(plot=g1, filename=paste0(nombre, 'mapavb-lpi.pdf'), width=17.53, height = 17.53, units='cm')
ggsave(plot=g2, filename=paste0(nombre, 'mapahmc-lpi.pdf'), width=17.53, height = 17.53, units='cm')

# (4) Modelo probit-ordenado para estimación ====    
## (4.1) Seleccionar datos ====
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
    TRUE ~ 0)# no carencias
)
table(tmp$ym_ord)
# no invertir
# tmp$ym_ord <- ifelse(tmp$ym_ord == 0, 2, ifelse(tmp$ym_ord == 2, 0, 1))
set.seed(1)
datos <- get_obs(censo=T, y=tmp$ym_ord, prop=0.2)

## (4.2) Cargar el modelo con VB ====
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-probit-ordenado-sesgado-estimacion.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-probit-ordenado-sesgado-8020.RData)')
muestra_vb <- muestra_probitordVB$summary(variables = c('rho', 'delta1', 'mu', 'b'), out_stan)

## (4.4) Cargar el modelo con HMC ====
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-probit-ordenado-sesgado-estimacion.RData)')
load(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-probit-ordenado-sesgado-8020.RData)')
muestra_hmc <- muestra_probitordHMC$summary(variables = c('rho', 'delta1', 'mu', 'b'), c(out_stan, 'rhat'))
colnames(muestra_hmc) <- paste0(colnames(muestra_hmc), 2)

## (4.5) Tabla de estimaciones ====
parametros <- c(
  paste0('$\\rho_{', 1:K, '}$'),
  '$\\delta_{1}$',
  paste0('$\\mu_{', 1:K, '}$')
)

pred_out <- tibble(parametros,
                   muestra_vb[1:length(parametros), c(2, 6, 7)],
                   muestra_hmc[1:length(parametros), c(2, 6, 7, 8)])

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-ordenado-sesgado-tabla.txt')
writeLines(pred_out %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

## (4.6) Métricas ====
# EntornosR
# EntornosR-8020
out <- datos$idx
tproc_vb <- readLines(con=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraVB-tiempos-8020.txt)') %>%
  str_split_1(pattern='~') %>% as.numeric()
tproc_hmc <- readLines(con=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\EntornosR\muestraHMC-tiempos-8020.txt)') %>%
  str_split_1(pattern='~') %>% as.numeric()
delta_vb <- c(0, unlist(muestra_vb[17, 2]))
delta_hmc <- c(0, unlist(muestra_hmc[17, 2]))
yhat_vb <- findInterval(apply(muestra_probitordVB$draws('z_all'), 2, mean), delta_vb)[out]
yhat_hmc <- findInterval(apply(drop(muestra_probitordHMC$draws('z_all')), 2, mean), delta_hmc)[out]
# no invertir
# yhat_vb <- ifelse(yhat_vb==2, 0, ifelse(yhat_vb==0, 2, 1)) # no lo corri asi
# yhat_hmc <- ifelse(yhat_hmc==2, 0, ifelse(yhat_hmc==0, 2, 1)) # omitir aqui
# yhat_hmc <- yhat_vb
tmp <- cbind(metrics_cat(yhat_vb, datos$y_test, F), metrics_cat(yhat_hmc, datos$y_test,F))
tmp <- rbind(tmp, c(tproc_vb[3], tproc_hmc[3]))
rownames(tmp)[2] <- '$\\tau$ de Kendall'
rownames(tmp)[3] <- 'Tiempo (s)'

met <- apply(tmp, 1, which.max)
met[nrow(tmp)] <- ifelse(met[nrow(tmp)]==1, 2, 1)
for(k in 1:nrow(tmp)){
  tmp[k, met[k]] <- paste0('\\bfseries ', tmp[k, met[k]])
}

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-ordenado-sesgado-metricas8020.txt')
writeLines(tmp %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:3, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

tmp <- tibble(Observados=datos$y_test,
              AjustadosVB=yhat_vb, AjustadosHMC=yhat_hmc)
tmp <- tmp %>%
  pivot_longer(cols=starts_with('Ajustados'),
               names_to = 'Metodo',
               names_prefix = 'Ajustados',
               values_to = 'Ajustados') %>%
  mutate(Metodo=factor(Metodo, levels=c('VB', 'HMC')))

tmp1 <- subset(tmp, Metodo=='VB')
cm1 <- caret::confusionMatrix(factor(tmp1$Observados, levels=c(2, 1, 0)),
                              factor(tmp1$Ajustados, levels=c(2, 1, 0)))
tmp2 <- subset(tmp, Metodo=='HMC')
cm2 <- caret::confusionMatrix(factor(tmp2$Observados, levels=c(2, 1, 0)),
                              factor(tmp2$Ajustados, levels=c(2, 1, 0)))

cm <- rbind(as_tibble(cm1$table), as_tibble(cm2$table))
cm <- cm %>%
  mutate(Metodo = factor(rep(c('VB', 'HMC'), each=9), levels=c('VB', 'HMC')),
         Tipo=rep(c('Verdaderos (2)', '', '',
                    '', 'Verdaderos (1)', '',
                    '', '', 'Verdaderos (0)'), 2))

g0 <- ggplot(cm%>%filter(Metodo=='VB'), aes(Prediction, Reference, fill = n)) +
  geom_tile(color='gray85', lwd=1.2) +
  scale_x_discrete(limits = rev)+
  geom_text(aes(label = paste0(Tipo, '\n', n), family = 'serif')) +
  scale_fill_gradient(low = "white", high = "gray50") +
  facet_wrap(~Metodo)+
  theme_minimal() +
  theme(text=element_text(family='serif', size=20), legend.position = 'none') +
  labs(x=expression(Valores ~~ y[ij] ~~ ajustados),
       y=expression(Valores ~~ y[ij] ~~ observados),
       fill=expression(n)) +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                labels=seq(0, 3000, by=500),
                                barheight = 20))
x11(); print(g0)

g1 <- ggplot(cm%>%filter(Metodo=='HMC'), aes(Prediction, Reference, fill = n)) +
  geom_tile(color='gray85', lwd=1.2) +
  scale_x_discrete(limits = rev)+
  geom_text(aes(label = paste0(Tipo, '\n', n), family = 'serif')) +
  scale_fill_gradient(low = "white", high = "gray50") +
  facet_wrap(~Metodo)+
  theme_minimal() +
  theme(text=element_text(family='serif', size=20), legend.position = 'none') +
  labs(x=expression(Valores ~~ y[ij] ~~ ajustados),
       y=expression(Valores ~~ y[ij] ~~ observados),
       fill=expression(n)) +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                labels=seq(0, 3000, by=500),
                                barheight = 20))
x11(); print(g1)

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-ordenado-sesgado-scatter8020')

ggsave(plot=g0, filename=paste0(nombre, 'VB.pdf'), width=20, height=25, units='cm')
ggsave(plot=g1, filename=paste0(nombre, 'HMC.pdf'), width=20, height=25, units='cm')

## (4.7) Porcentaje de personas ====
delta_vb <- c(0.0, unlist(muestra_vb[17, 2]))
delta_hmc <- c(0.0, unlist(muestra_hmc[17, 2]))
out <- which(is.na(ym))
fitted_vb <- findInterval(apply(muestra_probitordVB$draws('z_all'), 2, mean), delta_vb)
fitted_hmc <- findInterval(apply(drop(muestra_probitordHMC$draws('z_all')), 2, mean), delta_hmc)
# no invertir
# fitted_vb <- ifelse(fitted_vb==0, 2, 
#              ifelse(fitted_vb==2, 0, 1))
# fitted_hmc <- ifelse(fitted_hmc==0, 2, 
#                     ifelse(fitted_hmc==2, 0, 1))
# fitted_hmc <- fitted_vb
# dos: sin carencias
# uno: una carencia o LPI
# cero: dos carencias o LPEI

dftmp <- tibble(LPI_vb=ifelse(fitted_vb[out]==1, 1, 0),
                LPI_hmc=ifelse(fitted_hmc[out]==1, 1, 0),
                LPEI_vb=ifelse(fitted_vb[out]==0, 1, 0),
                LPEI_hmc=ifelse(fitted_hmc[out]==0, 1, 0),
                Ámbito=Xide$rururb[out],
                Factor=Xide$factor[out],
                Tamaño=Xide$tam_hog[out],
                Alcaldia=muni[region][out])
#
dftmp2 <- dftmp %>%
  group_by(Alcaldia) %>%
  summarise(Total = sum(Factor*Tamaño),
            Total_LPI_vb = sum(Factor*LPI_vb*Tamaño),
            Total_LPI_hmc = sum(Factor*LPI_hmc*Tamaño),
            Total_LPEI_vb = sum(Factor*LPEI_vb*Tamaño),
            Total_LPEI_hmc = sum(Factor*LPEI_hmc*Tamaño),
            Porcentaje_LPI_vb = 100*Total_LPI_vb/Total,
            Porcentaje_LPEI_vb = 100*Total_LPEI_vb/Total,
            Porcentaje_LPI_hmc = 100*Total_LPI_hmc/Total,
            Porcentaje_LPEI_hmc = 100*Total_LPEI_hmc/Total)

dftmp2 %>% select(Total, starts_with('Porcentaje'))
sum(dftmp2$Total_LPI_vb) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPI_hmc) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPEI_vb) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPEI_hmc) / sum(dftmp2$Total) # 25.7 % 

personas <- dftmp2 %>% select(Alcaldia, starts_with('Porcentaje') & ends_with('vb'))
personas <- dftmp2 %>% select(Alcaldia, starts_with('Porcentaje'))
personas <- personas[match(muni, personas$Alcaldia), ]

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-ordenado-sesgado-personas.txt')
writeLines(personas %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

## (4.8) Mapas ====
orden <- match(dftmp2$Alcaldia, shp$NOMGEO) # el elemento i (i=1, 2, ...) de #1 esta en la posición indicada de #2
orden2 <- match(muni, shp$NOMGEO) # el elemento i (i=1, 2, ...) de #1 esta en la posición indicada de #2
shp$LPI_vb_porcentaje[orden] <- dftmp2$Porcentaje_LPI_vb
shp$LPEI_vb_porcentaje[orden] <- dftmp2$Porcentaje_LPEI_vb
shp$LPI_hmc_porcentaje[orden] <- dftmp2$Porcentaje_LPI_hmc
shp$LPEI_hmc_porcentaje[orden] <- dftmp2$Porcentaje_LPEI_hmc

g1 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPI_vb_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    text = element_text(family = "serif", size=20),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g1)

g2 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPEI_vb_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPEI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    text = element_text(family = "serif", size=20),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g2)

g3 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPI_hmc_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    text = element_text(family = "serif", size=20),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g3)

g4 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPEI_hmc_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPEI (%)', title = "", x='', y='') + 
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    text = element_text(family = "serif", size=20),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g4)

nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-ordenado-sesgado-')
ggsave(plot=g1, filename=paste0(nombre, 'mapavb-lpi.pdf'), width=17.53, height = 17.53, units='cm')
ggsave(plot=g2, filename=paste0(nombre, 'mapavb-lpei.pdf'), width=17.53, height = 17.53, units='cm')
ggsave(plot=g3, filename=paste0(nombre, 'mapahmc-lpi.pdf'), width=17.53, height = 17.53, units='cm')
ggsave(plot=g4, filename=paste0(nombre, 'mapahmc-lpei.pdf'), width=17.53, height = 17.53, units='cm')

# (3* y 4*) Distribucion de 0's, 1's, 2's ====
# 0: urbano, 1: rural
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
    ym_ord <= 1 ~ 1, # carencias
    is.na(ym_ord) ~ NA,
    TRUE ~ 0) # no carencias
)
#renglones: primer argumento: clasificacion
#columnas: segundo argumento: contexto
# hacer dos matrices de conteos
table(tmp$ym_bin)
table(tmp$ym_bin, tmp$contexto)
cm1 <- tibble(
  Categoria = factor(c(0, 0, 1, 1)),
  Contexto = factor(c(0, 1, 0, 1)),
  n = c(1596, 443, 388, 110))

table(tmp$ym_ord, tmp$contexto)
cm2 <- tibble(
  Categoria = factor(c(0, 0, 1, 1, 2, 2)),
  Contexto = factor(c(0, 1, 0, 1, 0, 1)),
  n = c(69, 13, 319, 97, 1596, 443))

g0 <- ggplot(cm1, aes(Contexto, Categoria, fill=n)) +
  geom_tile(color='gray85', lwd=1.2) +
  scale_x_discrete(limits = rev) + 
  scale_fill_gradient(low = "white", high = "gray75") +
  theme_minimal() +
  geom_text(aes(label = paste0(n), family = 'serif'), size=26, size.unit = 'pt') + 
  theme(text=element_text(family='serif', size=26), legend.position = 'none') +
  labs(x=expression(Contexto),
       y=expression(Clasificación),
       fill=expression(n)) +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                labels=seq(0, 2000, by=500),
                                barheight = 20))
x11(); print(g0)
ggsave(plot=g0, width=20, height=20, units='cm',
       filename=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX\dist-probit.pdf)')

g1 <- ggplot(cm2, aes(Contexto, Categoria, fill=n)) +
  geom_tile(color='gray85', lwd=1.2) +
  scale_x_discrete(limits = rev) + 
  scale_fill_gradient(low = "white", high = "gray75") +
  geom_text(aes(label = paste0(n), family = 'serif'), size=26, size.unit = 'pt') + 
  theme_minimal() +
  theme(text=element_text(family='serif', size=26), legend.position = 'none') +
  labs(x=expression(Contexto),
       y=expression(Clasificación),
       fill=expression(n)) +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                labels=seq(0, 2000, by=500),
                                barheight = 20))
x11(); print(g1)
ggsave(plot=g1, width=20, height=20, units='cm',
       filename=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX\dist-probitord.pdf)')

# (5) Estimaciones de beta ====
betas_vb <- muestra_probitordVB$summary(variables = c('b'), out_stan)[, c(2, 6, 7)]
betas_hmc <- muestra_probitordHMC$summary(variables = c('b'), c(out_stan, 'rhat'))[, c(2, 6, 7, 8)]
sum(marginal_pi(muestra_probitordVB)>75)
sum(marginal_pi(muestra_probitordHMC)>75)
posterior_prob(muestra_probitVB)[1:3, ]
tmp <- cbind(betas_vb, marginal_pi(muestra_probitordVB), betas_hmc, c(marginal_pi(muestra_probitordHMC)))
tmp <- cbind(betas_vb, betas_hmc)
idx <- intersect(which(marginal_pi(muestra_probitordVB)>75), which(marginal_pi(muestra_probitordHMC)>75))
rownames(tmp) <- paste0('$\\beta_{',1:nrow(tmp), '}$')
tmp <- tmp[idx, ]
dim(tmp)
colnames(tmp) <- NULL
nombre <- file.path(fsep='\\', r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX)',
                    'probit-ordenado-sesgado-beta.txt')
writeLines(tmp %>% kable(format = 'latex', escape = F), con = nombre)
tmp <- readLines(nombre)
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, nombre)

## (4.4) Precisión y tiempo (ajustados) ====
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

# (6) Estimaciones de los parametros ====

## (3.1) Estimación de beta ====
beta_ <- muestra_VB$draws("b")
Media <- apply(beta_, 2, mean)
Error <- apply(beta_, 2, sd)
Lo025 <- apply(beta_, 2, quantile, 0.025)
Up975 <- apply(beta_, 2, quantile, 0.975)
Frecuencia <- 100*apply(muestra_VB$draws("m_ind"), 2, mean)
Parametro <- paste0("$\\beta_{", 1:ncol(datos$X_obs), "}$")

beta_out <- tibble(Parametro, Media, `Error est.`=Error, `0.025 \\%`=Lo025, `0.975 \\%`=Up975, `Frecuencia \\%`=Frecuencia)
print(beta_out, n=Inf)

# (6) Porcentaje de personas ====

# (7) Lista de covariables ====
library(readxl)
library(tidyr)
library(stringr)
load('./Xyregion-CDMX.Rdata')
descripcion <- read_excel('descripcion_variables.xlsx', range='A1:D308')

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

Tipo <- c(rep('Continua', ncol(Xcon)), rep('Binaria', ncol(Xbin)), rep('Categórica', ncol(Xcat)))
nomall <- cbind(nomall, Tipo)

tabla_cov <- nomall %>%
  select(!c('Valores', 'Etiquetas', 'Variable')) %>%
  relocate(Tipo) %>%
  select('Descripción') %>%
  tibble
# tabla_cov <- tibble(1:nrow(tabla_cov), tabla_cov)
tabla_cov[, 1] <- paste0('\\noindent ', 1:nrow(tabla_cov))
dim(tabla_cov)
tabla_cov
140/4 # sideways, 4 cols de 35?
# 1. bla bla bla
# 2. ble ble ble
# ...
library(kableExtra)
nombre <- r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-v)'
col1 <- tabla_cov[1:70, ] %>% kable(format = 'latex', escape = F)
col1 <- str_replace_all(col1, "\\\\hline", "")
col1 <- str_replace_all(col1, "\\\\\\\\", "")
col1 <- str_replace_all(col1, " & ", ". ")
col2 <- tabla_cov[71:140, ] %>% kable(format = 'latex', escape = F)
col2 <- str_replace_all(col2, "\\\\hline", "")
col2 <- str_replace_all(col2, "\\\\\\\\", "")
col2 <- str_replace_all(col2, " & ", ". ")
writeLines(col1,con = file.path(nombre, 'lista1.txt',fsep = '\\'))
writeLines(col2,con = file.path(nombre, 'lista2.txt',fsep = '\\'))
# lista 1
tmp <- readLines(file.path(nombre, 'lista1.txt',fsep = '\\'))
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, file.path(nombre, 'lista1.txt',fsep = '\\'))
# lista 2
tmp <- readLines(file.path(nombre, 'lista2.txt',fsep = '\\'))
tmp <- tmp[-c(1:5, length(tmp)-1, length(tmp))]
tmp[length(tmp)] <- str_split_1(tmp[length(tmp)], pattern = '\\\\\\\\')[1]
writeLines(tmp, file.path(nombre, 'lista2.txt',fsep = '\\'))

tabla_cov <- nomall %>%
  select(!c('Valores', 'Etiquetas')) %>%
  relocate(Tipo) %>%
  tibble %>% # hasta aqui para todas
  filter(Variable %in% colnames(Xtmp))

tmp <- tibble(Tipo=rep(NA, 26), Variable=rep(NA, 26), Descripción=rep(NA, 26))
tabla_cov <- rbind(tmp, tabla_cov)
print(tabla_cov, n=Inf)

# (8) Covariables y beta ====
muestra_VB
idx <- which(beta_out$`Frecuencia \\%` > 90); length(idx)
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
  Descripción = c(rep('', times=26)[idx[idx<=26]], tabla_cov$Descripción[idx[idx>26]])
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

# 4.6 Descriptivos ====

datos <- list()
datos$master <- which(!is.na(ym))
dftmp <- data.frame(Observados=ym[datos$master],
                    Ámbito=Xide$rururb[datos$master],
                    Factor=Xide$factor[datos$master],
                    Alcaldia=muni[region[datos$master]])

dftmp %>% group_by(Alcaldia) %>%
  summarise(Minimo=min(exp(Observados)),
            Mediana = median(exp(Observados)),
            Media = mean(exp(Observados)),
            Maximo=max(exp(Observados)),
            `Desv. est.`=sd(exp(Observados)),
            Conteo=n(), 
            .groups = "drop")

table_out <- dftmp %>%
  mutate(Ámbito=factor(Ámbito, labels=c('Urbano', 'Rural')),
         Alcaldia=factor(Alcaldia, levels=muni)) %>%
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

Res_ambito <- dftmp %>%
  group_by(Ámbito) %>%
  summarise(Minimo=min(exp(Observados)),
            Mediana = median(exp(Observados)),
            Media = mean(exp(Observados)),
            Maximo=max(exp(Observados)),
            `Desv. est.`=sd(exp(Observados)),
            Conteo=n(), 
            .groups = "drop")

## Cuadro 4.32 Estadisticos por alcaldia ----
print(table_out, n=Inf)
table_out %>%
  mutate(Alcaldia=abbreviate(Alcaldia, minlength = 3, strict = F, dot = T)) %>%
  kable(format='latex')

## Cuadro 4.33. Estadisticos por ámbito ----
print(Res_ambito)


## Figura 4.34 Gráfico de distribución del ingreso (histogramas) ----
umbral <- 400 # filtrar tres ingresos muy pequeños
tmp <- which(!is.na(ym) & ym>log(umbral)) # o solo 10
dftmp <- data.frame(y = ym[tmp], mun = muni[region[tmp]],
                    amb = factor(Xide$rururb[tmp], labels=c('Urbano', 'Rural')))

y_hist <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray25', fill='gray90') +
  geom_density(linewidth=1.2, color='gray50', linetype=2) + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif", size = 18))
x11(); plot(y_hist)

y_hist2 <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray25', fill='gray90', bins = 30) +
  geom_density(linewidth=1, color='gray50', linetype=2) + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif", size = 18)) + 
  facet_wrap(~ amb, ncol = 2, nrow = 1, scales = 'fixed')
x11(); plot(y_hist2)

ggsave(r"(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX\dist-log-ict.pdf)", 
       plot = y_hist, width = 15, height = 15, units = "cm")

ggsave(r"(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-vi\CDMX\dist-log-ict-ambito.pdf)", 
       plot = y_hist2, width = 15, height = 15, units = "cm")

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



