
library(dplyr)
library(ggplot2)
library(tidyr)
library(rlang) # rep_along

setwd(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras)')

source(r'(../Codigos/Simulacion-Modelos/funciones.R)')

# setwd(r'(C:\Maestria\Articulo-GitHub\Figuras)') # guardar en articulo

# una figura de 20 x 20 cm tendrá size = 26, para un buen escalado

# Cap 2. Revisión de literatura ====

# Figura 2.1 normal sesgada ====

dsn <- function(x0, mu=0, sigma=1, lambda=0){
  y0 <- 2/sigma * dnorm((x0-mu)/sigma, 0, 1) * pnorm(lambda * (x0-mu)/sigma)
  return(y0)
}

x0 <- seq(-3, 3, length=250)
y0 <- dsn(x0, lambda=0)
y1 <- dsn(x0, lambda=1)
y2 <- dsn(x0, lambda=-1)
y3 <- dsn(x0, lambda=3)
y4 <- dsn(x0, lambda=-3)

lambda <- rep(c(0, 1, -1, 3, -3), each=length(x0))
graphs1 <- tibble(x0=rep(x0, times=5), y0=c(y0, y1, y2, y3, y4), Lambda=factor(lambda))
graphs1$Lambda <- factor(graphs1$Lambda, levels = c("0", "1", "3", '-1', '-3'))

g1 <- ggplot(data=graphs1 %>% filter(Lambda%in%c(3, 1, 0)),
             aes(x=x0, y=y0, color=Lambda, linetype = Lambda)) + 
  geom_line(linewidth=1.2) + 
  scale_color_manual(values=gray.colors(n=3, end=0.5)) +
  scale_linetype_manual(values=1:3) +
  # scale_color_manual(values = setNames(gray.colors(n=3, end=0.5), c(0, 1, 3)), name = "Lambda") +
  # scale_linetype_manual(values = setNames(1:3,  c(0, 1, 3)), name = "Lambda") +
  geom_hline(yintercept = 0)+
  theme_minimal() + 
  theme(legend.position = 'top', text = element_text(family = "serif", size=26)) + 
  labs(x=expression(x), y = '' , color = bquote(Valor~de~lambda~':')) + 
  guides(
    color = guide_legend(override.aes = list(linetype = 1:3, size = 1.2)),
    linetype = "none"
  ) 
x11(); plot(g1)

g2 <- ggplot(data=graphs1 %>% filter(Lambda%in%c(-3, -1, 0)),
             aes(x=x0, y=y0, color=Lambda, linetype = Lambda)) + 
  geom_line(linewidth=1.2) + 
  scale_color_manual(values=gray.colors(n=3, end=0.5)) +
  scale_linetype_manual(values=1:3) +
  geom_hline(yintercept = 0)+
  theme_minimal() + 
  theme(legend.position = 'top', text = element_text(family = "serif", size=26)) + 
  labs(x=expression(x), y = '' , color = bquote(Valor~de~lambda~':')) + 
  guides(
    color = guide_legend(override.aes = list(linetype = 1:3, size = 1.2)),
    linetype = "none"
  ) 
x11(); plot(g2)

getwd()
ggsave(plot=g1, filename='./c-ii/SN_density_pos.pdf', width=20, height = 20, units='cm')
ggsave(plot=g2, filename='./c-ii/SN_density_neg.pdf', width=20, height = 20, units='cm')

# Figura 2.2 normal sesgada multivariada (contornos) ====

# Figura 2.3 normal sesgada multivariada (superficie) ====
# Se omitió

# Figura 2.4 normal sesgada centrada ====

dcsn <- function(x, mu=0, sigma=1, lambda=0){
  s <- (2/(4-pi))**(1/3)
  r <- sqrt(2/pi)
  # lambda a gamma1
  gamma1 <- (4-pi)/2 * (lambda2rho(lambda)*r)**3 / (1 - 2*lambda2rho(lambda)**2 / pi)**(3/2)
  gamma1 <- abs(gamma1)
  if(lambda!=0){
    sig <- sign(lambda)
  } else{
    sig <- 1
  }
  #
  sigma_star <- sigma**2 * (1 + s**2 * gamma1**(2/3))
  lambda_star <- s * sig*gamma1**(1/3) / sqrt(r**2 + s**2 * gamma1**(2/3) * (r**2 - 1))
  mu_star <- mu - s * sig*gamma1**(1/3)
  f <- 2/sigma_star * dnorm((x-mu_star)/sigma_star)*pnorm(lambda_star * (x-mu_star)/sigma_star)
  return(f)
}

x0 <- seq(-4, 4, length=250)
y0 <- dcsn(x0, lambda=0)
y1 <- dcsn(x0, lambda=1)
y2 <- dcsn(x0, lambda=-1)
y3 <- dcsn(x0, lambda=3)
y4 <- dcsn(x0, lambda=-3)

lambda <- rep(c(0, 1, -1, 3, -3), each=length(x0))
graphs1 <- tibble(x0=rep(x0, times=5), y0=c(y0, y1, y2, y3, y4), Lambda=factor(lambda))
graphs1$Lambda <- factor(graphs1$Lambda, levels = c("0", "1", "3", '-1', '-3'))

g1 <- ggplot(data=graphs1 %>% filter(Lambda%in%c(3, 1, 0)),
             aes(x=x0, y=y0, color=Lambda, linetype = Lambda)) + 
  geom_line(linewidth=1.2) + 
  scale_color_manual(values=gray.colors(n=3, end=0.5)) +
  scale_linetype_manual(values=1:3) +
  # scale_color_manual(values = setNames(gray.colors(n=3, end=0.5), c(0, 1, 3)), name = "Lambda") +
  # scale_linetype_manual(values = setNames(1:3,  c(0, 1, 3)), name = "Lambda") +
  geom_hline(yintercept = 0)+
  theme_minimal() + 
  theme(legend.position = 'top',
        text = element_text(family = "serif", size = 26)) + 
  labs(x=expression(x), y = '' , color = bquote(Valor~de~lambda~':')) + 
  guides(
    color = guide_legend(override.aes = list(linetype = 1:3, size = 1.2)),
    linetype = "none"
  ) 
x11(); plot(g1)

g2 <- ggplot(data=graphs1 %>% filter(Lambda%in%c(-3, -1, 0)),
             aes(x=x0, y=y0, color=Lambda, linetype = Lambda)) + 
  geom_line(linewidth=1.2) + 
  scale_color_manual(values=gray.colors(n=3, end=0.5)) +
  scale_linetype_manual(values=1:3) +
  geom_hline(yintercept = 0)+
  theme_minimal() + 
  theme(legend.position = 'top',
        text = element_text(family = "serif", size=26)) + 
  labs(x=expression(x), y = '' , color = bquote(Valor~de~lambda~':')) + 
  guides(
    color = guide_legend(override.aes = list(linetype = 1:3, size = 1.2)),
    linetype = "none"
  ) 
x11(); plot(g2)

ggsave(plot=g1, filename='./c-ii/CSN_density_pos.pdf', width=20, height = 20, units='cm')
ggsave(plot=g2, filename='./c-ii/CSN_density_neg.pdf', width=20, height = 20, units='cm')


# Figura 2.7: rho y lambda ====
x0 <- seq(-2.5, 2.5, length=250) 
y0 <- lambda2rho(x0)
graphs0 <- tibble(x0=x0, y0=y0)

g0 <- ggplot(data=graphs0) +
  geom_line(aes(x=x0, y=y0, colour = "rho(lambda)"), linewidth=1.4, alpha=1) + # , arrow = arrow(length = unit(0.14, "inches"), type = "closed", ends = "both")
  geom_line(aes(x=y0, y=x0, colour = "lambda(rho)"), linewidth=1.4, alpha=1) + # , arrow = arrow(length = unit(0.14, "inches"), type = "closed", ends = "both")
  geom_abline(slope = 1, intercept = 0, linetype=2, color='gray5', linewidth=1.2) +
  theme_minimal() + 
  scale_colour_manual(
    name = "",                                   # legend title (empty here)
    breaks = c("rho(lambda)", "lambda(rho)"),    # order of legend entries
    values = c("rho(lambda)"    = "gray25",
               "lambda(rho)"    = "gray75"),
    labels = c(expression(rho(lambda)),          # use expression() if you want math
               expression(lambda(rho)))
  ) + 
  theme(legend.position = 'top', text = element_text(family = "serif")) + 
  labs(x='', y='') + 
  coord_fixed(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), ratio = 1)

x11(); print(g0)
# ggsave(plot=g0, filename='./c-ii/rho-lambda.pdf', width=20, height = 20, units='cm')

# Cap 3. Metodología ====

# Figura log-normal sesgada ====
dlogsn <- function(x0, mu=0, sigma=1, lambda=0){
  y0 <- 2/(x0*sigma) * dnorm((log(x0)-mu)/sigma) * pnorm(lambda * (log(x0)-mu)/sigma)
  return(y0)
}

x0 <- seq(0.001, 6, length=200)
dftmp <- tibble(x0=x0, y0=dlogsn(x0, lambda=0), y1=dlogsn(x0, lambda=1), y2=dlogsn(x0, lambda=3),
                       y3=dlogsn(x0, lambda=-1), y4=dlogsn(x0, lambda=-3))
graphs2 <- dftmp %>%
  pivot_longer(
    cols = starts_with("y"),
    names_to = "Lambda",
    values_to = "y"          # the density values
  ) %>%
  mutate(
    Lambda = recode(Lambda, y0 = 0, y1 = 1, y2 = 3, y3 = -1, y4 = -3),
    Lambda = factor(Lambda, labels = c(0, -1, 1, -3, 3))
  )

g3 <- ggplot(data=graphs2 %>% filter(Lambda%in%c(3, 1, 0)),
             aes(x=x0, y=y, color=Lambda, linetype = Lambda)) + 
  geom_line(linewidth=1.2) + 
  scale_color_manual(values=gray.colors(n=3, end=0.5)) +
  scale_linetype_manual(values=1:3) +
  # scale_color_manual(values = setNames(gray.colors(n=3, end=0.5), c(0, 1, 3)), name = "Lambda") +
  # scale_linetype_manual(values = setNames(1:3,  c(0, 1, 3)), name = "Lambda") +
  geom_hline(yintercept = 0)+
  theme_minimal() + 
  theme(legend.position = 'top', text = element_text(family = "serif", size=26)) + 
  labs(x=expression(x), y = '' , color = bquote(Valor~de~lambda~':')) + 
  guides(
    color = guide_legend(override.aes = list(linetype = 1:3, size = 1.2)),
    linetype = "none"
  ) 
x11(); plot(g3)

g4 <- ggplot(data=graphs2 %>% filter(Lambda%in%c(-3, -1, 0)),
             aes(x=x0, y=y, color=Lambda, linetype = Lambda)) + 
  geom_line(linewidth=1.2) + 
  scale_color_manual(values=gray.colors(n=3, end=0.5)) +
  scale_linetype_manual(values=1:3) +
  geom_hline(yintercept = 0)+
  theme_minimal() + 
  theme(legend.position = 'top', text = element_text(family = "serif", size=26)) + 
  labs(x=expression(x), y = '' , color = bquote(Valor~de~lambda~':')) + 
  guides(
    color = guide_legend(override.aes = list(linetype = 1:3, size = 1.2)),
    linetype = "none"
  ) 
x11(); plot(g4)

ggsave(plot=g3, filename='./logSN_density_pos.png', width=20, height = 20, units='cm')
ggsave(plot=g4, filename='./logSN_density_neg.png', width=20, height = 20, units='cm')


# Figura 3.2 SSVS ====
tau2=0.8
c2=4

delta <- function(tau2, c2){
  d <- sqrt(tau2) * sqrt(2*c2*log(sqrt(c2)) / (c2 - 1))
  return(d)
}

delta(tau2, c2)

t0 <- seq(-5, 5, length=250)
y0 <- dnorm(t0, mean=0, sqrt(tau2))
y1 <- dnorm(t0, mean=0, sqrt(tau2*c2))
densidad <- factor(rep(c('Spike', 'Slab'), each=length(t0)))
dftmp <- tibble(beta=c(t0, t0), y=c(y0, y1), Densidad=densidad)

p1 <- delta(tau2, c2)
p2 <- -delta(tau2, c2)

g0 <- ggplot(dftmp, aes(x=beta, y=y, colour = Densidad, linetype = Densidad)) + 
  geom_line(linewidth=1.2) + 
  theme_minimal() + 
  scale_x_continuous(labels = function(breaks){for(b in seq_along(breaks)){pos <- numeric(); if(breaks[b]==0){pos <- b}else{''}; out <- rep_along(breaks, ''); out[b+2] <- '0.0'; return(out)}}) +
  theme(legend.position = 'top', text = element_text(family = "serif", size=26)) + labs(x=expression(beta[k])) + 
  geom_hline(yintercept = 0, color='gray75')+
  labs(y='') + 
  scale_color_manual(values=c('Spike'='gray75', 'Slab'='gray25')) + 
  geom_text(x=p1, y=dnorm(p1, sd=sqrt(tau2)), label='delta[i]', parse=T, color='gray25', hjust=2, vjust=1, size = 8) +
  geom_text(x=p2, y=dnorm(p2, sd=sqrt(tau2)), label='-delta[i]', parse=T, color='gray25', hjust=2, vjust=1, size = 8) 

x11(); print(g0)

getwd()
ggsave(plot=g0, filename='./c-v/SSVS-prior.pdf', width=20, height = 20, units='cm')

# Figura 3.3 SSVS - Spike Slab ====

spike_slab <- function(x, mu1=0, mu2=0, sd1=1, sd2=1, p=0.5){
  y <- p * dnorm(x, mean=mu1, sd=sd1) + (1-p) * dnorm(x, mean=mu2, sd=sd2)
}

y0 <- spike_slab(t0, sd1=tau2*c2, sd2=tau2, p=0.05)
y1 <- spike_slab(t0, sd1=tau2*c2, sd2=tau2, p=0.25)
y2 <- spike_slab(t0, sd1=tau2*c2, sd2=tau2, p=0.50)
y3 <- spike_slab(t0, sd1=tau2*c2, sd2=tau2, p=0.75)
y4 <- spike_slab(t0, sd1=tau2*c2, sd2=tau2, p=0.95)

dftmp <- tibble(beta=t0, y0, y1, y2, y3, y4)

dftmp <- dftmp %>%
  pivot_longer(
    cols = starts_with("y"),        # y0, y1, ..., y4
    names_to = "prop",            # temporary column holding "y0", "y1", ...
    values_to = "y"          # the density values
  ) %>%
  mutate(
    prop = recode(prop, y0 = "p = 0.05", y1 = "p = 0.25", y2 = "p = 0.50", y3 = "p = 0.75", y4 = "p = 0.95"),
    prop = factor(prop)
  )

dftmp

g1 <- ggplot(dftmp%>%filter(!(prop%in%c('p = 0.05', 'p = 0.95'))),
             aes(x=beta, y=y, colour = prop, linetype = prop)) + 
# g1 <- ggplot(dftmp, aes(x=beta, y=y, colour = prop, linetype = prop)) + 
  geom_line(linewidth=1.2) + 
  theme_minimal() + 
  # scale_color_manual(values = grey.colors(n=5, end = 0.5)) + 
  scale_x_continuous(labels = function(breaks){for(b in seq_along(breaks)){pos <- numeric(); if(breaks[b]==0){pos <- b}else{''}; out <- rep_along(breaks, ''); out[b+2] <- '0.0'; return(out)}}) +
  scale_color_manual(values=c('p = 0.05'='gray5', 'p = 0.25'='gray25', 'p = 0.50'='gray50','p = 0.75'='gray75','p = 0.95'='gray95'))+
  theme(legend.position = 'top', text = element_text(family = "serif", size=26)) +
  geom_hline(yintercept = 0, color='gray75')+
  labs(x=expression(beta[k]), y='', color=expression('Proporción' ~~ p[k] ~ ':'), linetype=expression('Proporción' ~~ p[k] ~ ':'))  

x11(); print(g1)

ggsave(plot=g1, filename='./c-v/SSVS-prior-prob.pdf', width=20, height = 20, units='cm')


# Resultados ====
## Figura 4.1. Experimento con mu y rho en el modelo probit ====
n <- 200
mu0 <- seq(-2.5, 2.5, by=0.1)
rho0  <- seq(0.01, 0.99, by=0.01) # mas menos rho: seq o -seq
rho0_mu0 <- expand.grid(rho0, mu0)

get_bin <- function(row, centred=F){
  row <- as.numeric(row)
  if(centred){
    Std <- sqrt(1 - 2*row[1]**2/pi)
    Exp <- row[1]*sqrt(2/pi) / Std
    z_i <- row[2] + (rskewnorm(n=n, mean=0, sigma=Sigma(rho=row[1], N=1)) - Exp) / Std
    pr <- 1 - sn::psn(Exp-row[2]*Std, alpha=rho2lambda(row[1]))
  } else{
    z_i <- row[2] + rskewnorm(n=n, mean=0, sigma=Sigma(rho=row[1], N=1))
    pr <- 1 - sn::psn(-row[2], alpha=rho2lambda(row[1]))
  }
  y1_i <- ifelse(z_i>0, 1, 0)
  y2_i <- rbinom(n=n,size=1, prob= pr)
  n0_1 <- length(which(y1_i==0))
  n1_1 <- length(which(y1_i==1))
  n0_2 <- length(which(y2_i==0))
  n1_2 <- length(which(y2_i==1))
  return(c(n0_1, n1_1, n0_2, n1_2))
}

set.seed(1)
conteos <- t(apply(rho0_mu0, 1, get_bin, centred=F))
head(conteos)

tabla <- tibble(cbind(rho0_mu0, conteos)) %>%
  rename(rho0=Var1, mu0=Var2, n0_1=3, n1_1=4, n0_2=5, n1_2=6)

tabla <- tabla %>%
  pivot_longer(
    cols = starts_with("n"),               # all n* columns
    names_to = c("outcome", "Metodo"),     # split name by "_"
    names_sep = "_") %>%
  pivot_wider(
    names_from = outcome,                  # n0 / n1 become columns
    values_from = value) %>%
  mutate(`Proporción` = n1/(n0+n1)) %>%
  mutate(
    Metodo = recode(Metodo, `1` = "Signo", `2` = "Bernoulli")) %>%
  select(rho0, mu0, Metodo, n0, n1, `Proporción`)         # order columns nicely

head(tabla)
subset(tabla, rho0>0.9 & mu0==0)
subset(tabla, rho0<-0.9 & mu0==0)

g0 <- ggplot(data=tabla, aes(x=(mu0), y=(rho0), fill=`Proporción`)) + 
  geom_raster() +
  # geom_tile(color='gray90', alpha=0.9) + 
  # scale_fill_distiller(palette = "grays", direction = -1) +
  scale_fill_gradient(low='gray95',high = 'gray15')+
  theme_minimal() + 
  labs(x=expression(mu), y=expression(rho), fill='Prop. de unos:       \n') +
  facet_grid(~Metodo) + 
  theme(text = element_text(family = "serif", size = 12),
        legend.position = 'none') 

x11(); print(g0)

ggsave(plot=g0, filename='./c-vi/probit-estudio-rho10-mu.pdf', width=15, height = 15, units='cm') # negative
ggsave(plot=g0, filename='./c-vi/probit-estudio-rho01-mu.pdf', width=15, height = 15, units='cm') # positive

# ggsave(plot=g0, filename='./c-vi/probit-estudio-rho10-mu-centred.pdf', width=20, height = 20, units='cm')
# ggsave(plot=g0, filename='./c-vi/probit-estudio-rho01-mu-centred.pdf', width=20, height = 20, units='cm')

## Figura 4.2 Funcion acumulada SN ====

sn_link <- function(x=0, lambda=0){
  p <- 1 - sn::psn(x=-x, xi=0, omega=1, alpha=lambda) # fn liga
  return(p)
}

t0 <- seq(-2.5, 2.5, length=200)
y0 <- sn_link(t0)
y1 <- sn_link(t0, lambda=-1) # +- 1
y2 <- sn_link(t0, lambda=-3) # +- 3
# x11(); plot(t0, y0, type='l', lwd=2)
# lines(t0, y1, col=4, lwd=2)
# lines(t0, y2, col=3, lwd=2)
# grid()

probit_link <- tibble(x0=t0, y0=y0,y1=y1, y2=y2) %>%
  pivot_longer(
    cols = starts_with("y"),
    values_to = 'y0',
    names_to = 'Densidad') %>%
  mutate(Densidad=factor(Densidad, levels=c('y0', 'y1', 'y2')))

g1 <- ggplot(probit_link, aes(x=x0, y=y0, color=Densidad, linetype=Densidad)) + 
  geom_line(linewidth=1.2) + 
  labs(colour = "Densidad", linetype = "Densidad") +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0, linetype=2, color='gray50')+
  theme_bw(base_size = 12) + 
  scale_linetype_manual(
    name=expression(lambda ~ ':'),
    labels=sprintf("%.1f", c(0, -1, -3)), # c(0, -1, -3)
    values=c('y0'=1, 'y1'=2, 'y2'=3)) + 
  labs(x=expression(eta[ij]), y=expression(Pr(y[ij] == 1 ~ "|" ~ p[ij])))+
  scale_color_manual(
    name=expression(lambda ~ ':'),
    labels=sprintf("%.1f", c(0, -1, -3)), # c(0, -1, -3)
    values=c('y0'='gray50', 'y1'='gray75', 'y2'='gray50'))+
  theme(legend.position = 'top',
        legend.text = element_text(size = 22),
        text=element_text(family='serif'))
x11(); print(g1)

# ggsave(plot=g1, filename='./c-vi/probit-fn-liga-rho01.pdf', width=20, height = 20, units='cm') # positive
# ggsave(plot=g1, filename='./c-vi/probit-fn-liga-rho10.pdf', width=20, height = 20, units='cm') # negative


# (5.3) divergencia KL(q, p) ====
# here, supp(p) < supp(q)
# so, q is normal and p is skew-normal
set.seed(1)
kl_mc <- function(mu_p=0, sd_p=1, mu_q=0, sd_q=1, delta=0, n = 1000){
  x <- sn::rsn(n, xi = mu_q, omega = sd_q, alpha=delta) # documento
  logq <- dnorm(x, mean = mu_q, sd = sd_q, log = TRUE)
  logp <- sn::dsn(x, xi=mu_p, omega=sd_p, alpha=delta, log = T)
  dif <- logp - logq # documento
  est <- mean(dif)
  return(est)
}

kl_mc <- function(mu_p=0, sd_p=1, mu_q=0, sd_q=1, delta=0, n = 1000){
  x <- rnorm(n, mean = mu_q, sd = sd_q)
  logp <- dnorm(x, mean = mu_q, sd = sd_q, log = TRUE)
  logq <- sn::dsn(x, xi=mu_p, omega=sd_p, alpha=delta, log = T)
  dif <- logp - logq
  est <- mean(dif)
  return(est)
}


KL0 <- c()
KL1 <- c()
lambda <- seq(-10, 10, length=200)
rho <- seq(-0.9999, 0.9999, length=200)
for(L in seq_along(lambda)){
  KL0[L] <- kl_mc(delta=lambda[L])
  KL1[L] <- kl_mc(delta=rho2lambda(rho[L]))
}

x11(); plot(lambda, KL0, type='l')
x11(); plot(rho, KL1, col=4)



# (6) rskewnorm ====
sigma2 <- 1.5
rho <- 0.7
n <- 10
ones <- rep(1, times=n)
Jn <- ones %*% t(ones)
A <- diag(sigma2, nrow=n) - rho**2*Jn
eigen(A)$values
sigma <- Sigma(rho=rho, N=5)
L <- 5
covar0 <- cbind(sigma[1:(L-1), L]) %*% solve(sigma[L, L])
covar0 <- (sigma[1:(L-1), L]) %*% solve(sigma[L, L])
sigma[1:(L-1), 1:(L-1)] - covar0 %*% rbind(sigma[L, 1:(L-1)])

n <- 500
rho <- c(0.95)
set.seed(2)
r1 <- rskewnorm(n=1, mean=0, sigma=Sigma(rho=rho, N=n)) # one multivariate skew-normal draw.
r2 <- rskewnorm(n=n, mean=0, sigma=Sigma(rho=rho, N=1)) # n univariate skew-normal draws
r3 <- sn::rsn(n=n, alpha=rho2lambda(rho))

f1 <- sn::selm(r1 ~ 1)
f2 <- sn::selm(r2 ~ 1)
f3 <- sn::selm(r3 ~ 1)
coef(f1, param.type = 'DP')
coef(f2, param.type = 'DP')
coef(f3, param.type = 'DP')


lambda2rho(coef(f1, param.type = 'DP')[3])
lambda2rho(coef(f2, param.type = 'DP')[3])
lambda2rho(coef(f3, param.type = 'DP')[3])

x11(); hist(r1)
x11(); hist(r2)
x11(); hist(r3)


# (7) en dos dimensiones ====

rskewnorm(n=10, mean=0, sigma=Sigma(rho=0.95, N=1, sigma=1.0))
rskewnorm(n=4, mean=0, sigma=Sigma(rho=0.05, N=2, sigma=1.0))
n <- 500
X <- rskewnorm(n=n, mean=0, sigma=Sigma(rho=lambda2rho(2), N=2, sigma=1.0))
X <- rskewnorm(n=n, mean=0, sigma=Sigma(rho=0.0001, N=2, sigma=1.0))
X <- matrix(X, nrow=n, byrow=T)
dim(X)
datosX <- data.frame(x1=X[, 1], x2=X[, 2])
x1 <- datosX$x1
x2 <- datosX$x2
lambda2rho(coef(sn::selm(x1 ~ 1), param.type='DP')[3])
lambda2rho(coef(sn::selm(x2 ~ 1), param.type='DP')[3])

set.seed(1)
datosX <- sn::rmsn(n=n, xi=c(0, 0), Omega=diag(2), alpha=c(2, 2))
colnames(datosX) <- c('x1', 'x2')
x1 <- datosX[, 1]
x2 <- datosX[, 2]
lambda2rho(coef(sn::selm(x1 ~ 1), param.type='DP')[3])
lambda2rho(coef(sn::selm(x2 ~ 1), param.type='DP')[3])

msn <- function(x, mu, Sigma, lambda){
  phi <- mvtnorm::dmvnorm(x=x, mean=mu, sigma=Sigma)
  sigma <- diag(diag(Sigma))
  z <- solve(sigma) %*% (x - mu)
  Phi <- pnorm(t(lambda) %*% z)
  f <- 2 * phi * Phi
  return(as.numeric(f))
}

msn_cent <- function(x, mu, Sigma, lambda){
  Ez  <- lambda2rho(lambda[1])*sqrt(2/pi)
  Sdz <- sqrt(1 - 2*lambda2rho(lambda[1])**2/pi)
  mu <- mu - Ez
  sigma <- diag(Sdz, ncol=ncol(Sigma), nrow=nrow(Sigma))
  Sigma <- solve(sigma) %*% Sigma %*% solve(sigma)
  phi <- mvtnorm::dmvnorm(x=x, mean=mu, sigma=Sigma)
  sigma <- diag(diag(Sigma))
  z <- solve(sigma) %*% (x - mu)
  Phi <- pnorm(t(lambda) %*% z)
  f <- 2 * phi * Phi
  return(as.numeric(f))
}

x1 <- seq(-3, 3, length=100)
x2 <- seq(-3, 3, length=100)
xg <- expand.grid(x1, x2)

## (6.1) Contornos con Sigma = sigma^2 I ====
# la desventaja es que cada entrada de lambda tiene el mismo valor de forma, es decir
# lambda_n = lambda 1_n
# la funcion rskewnorm no puede manejar Sigma = I
Sigma2 <- function(rho, N, sigma2=1.0){
  ones <- matrix(1L, nrow=N)
  J <- tcrossprod(ones)
  I <- diag(x=1L, nrow=N)
  block11 <- sigma2*I
  block12 <- ones*rho*sigma2
  block21 <- t(block12)
  block22 <- sigma2
  MatSigma <- rbind(cbind(block11, block12), cbind(block21, block22))
  return(MatSigma)
}
n <- 1000
X <- rskewnorm(n=n, mean=0, sigma=Sigma2(rho=-lambda2rho(2), N=2, sigma=1.0))
X <- matrix(X, nrow=n, byrow=T)
datosX <- data.frame(x1=X[, 1], x2=X[, 2])
datosY <- sn::rmsn(n=n, xi=c(0, 0), Omega=diag(2), alpha=-c(2, 2))
colnames(datosY) <- c('y1', 'y2')
fg <- apply(xg, 1, msn, mu=c(0, 0), Sigma=diag(2), lambda=-c(2.0, 2.0)) # son iguales
gg <- apply(xg, 1, sn::dmsn, xi=c(0, 0), Omega=diag(2), alpha=-c(2.0, 2.0)) # son iguales
datos <- tibble(xg, fg, gg)
g1 <- ggplot(datos, aes(x=Var1, y=Var2, z=fg, color='A')) +
  geom_contour() + 
  stat_density_2d(datosX, mapping=aes(x=x1, y=x2, color='B'),
                  inherit.aes = F, linewidth = 0.5,  linetype=2) +
  # stat_density_2d(datosY, mapping=aes(x=y1, y=y2),
  #                 inherit.aes = F, linewidth = 0.8, color='gray50') +
  theme_minimal() + labs(
    x = expression(x[1]),
    y = expression(x[2])) + 
  scale_color_manual(name=expression(SN[2](x ~ '|'~  mu ~ ','  ~ Sigma ~ ',' ~ lambda)),
                     values=c(A='gray35', B='gray50'),
                     labels=c('Contornos', 'Estimación kernel'))+
  theme(text=element_text(family='serif'), legend.position = 'top')
x11(); print(g1)
ggsave(plot=g1, filename=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-ii\SN2D-contour2.pdf)',
       width=20, height=20, units='cm')

## (6.2) Sigma = sigma^2 [(1-\rho^2)I + \rho^2 J] ====
n <- 1000
X <- rskewnorm(n=n, mean=0, sigma=Sigma(rho=lambda2rho(2), N=2, sigma=1.0))
X <- matrix(X, nrow=n, byrow=T)
datosX <- data.frame(x1=X[, 1], x2=X[, 2]) # generadas por mi
S <- Sigma(rho = lambda2rho(2), N = 2, sigma2 = 1)[1:2, 1:2]
datosY <- sn::rmsn(n=n, xi=c(0, 0), Omega=S, alpha=c(2, 2))
colnames(datosY) <- c('y1', 'y2')
# contornos de la densidad
fg <- apply(xg, 1, msn, mu=c(0, 0), Sigma=S, lambda=c(2.0, 2.0))
gg <- apply(xg, 1, sn::dmsn, xi=c(0, 0), Omega=S, alpha=c(2.0, 2.0))
datos <- tibble(xg, fg, gg)
tail(datos) # son iguales
# grafica: contornos estimados y teóricos
g1 <- ggplot(datos, aes(x=Var1, y=Var2, z=fg)) +
  geom_contour() + 
  stat_density_2d(datosX, mapping=aes(x=x1, y=x2),
                  inherit.aes = F, linewidth = 0.8, color='red') +
  stat_density_2d(datosY, mapping=aes(x=y1, y=y2),
                  inherit.aes = F, linewidth = 0.8, color='green') +
  theme_minimal() + labs(
    title = "Estimated density contours",
    x = expression(x[1]),
    y = expression(x[2]))
x11(); print(g1)

## (6.3) Contornos con Sigma = sigma^2 I, par. centrada ====
# parametros centrados
Ez  <- lambda2rho(-2)*sqrt(2/pi)
Sdz <- sqrt(1 - 2*lambda2rho(-2)**2/pi)
sigma <- diag(Sdz, ncol=ncol(Sigma), nrow=nrow(Sigma))
# contornos estimados
n <- 1000
S2 <- Sigma2(rho=-lambda2rho(2), N=2, sigma=1.0)
S2[1:2, 1:2] <- solve(sigma)%*%Sigma2(rho=lambda2rho(2), N=2, sigma=1.0)[1:2, 1:2]%*%solve(sigma)
X <- rskewnorm(n=n, mean=0-Ez, sigma=S2)
X <- matrix(X, nrow=n, byrow=T)
datosX <- data.frame(x1=X[, 1], x2=X[, 2])
datosY <- sn::rmsn(n=n, xi=c(0-Ez, 0-Ez), Omega=solve(sigma)%*%diag(2)%*%solve(sigma), alpha=c(2, 2))
colnames(datosY) <- c('y1', 'y2')
# contornos teoricos

fg <- apply(xg, 1, msn_cent, mu=c(0, 0), Sigma=diag(2), lambda=-c(2.0, 2.0))
gg <- apply(xg, 1, sn::dmsn, xi=c(0-Ez, 0-Ez), Omega=solve(sigma)%*%diag(2)%*%solve(sigma), alpha=c(2.0, 2.0))
datos <- tibble(xg, fg, gg)
tail(datos)  # son iguales
g1 <- ggplot(datos, aes(x=Var1, y=Var2, z=fg, color='A')) +
  geom_contour() + 
  stat_density_2d(datosX, mapping=aes(x=x1, y=x2, color='B'),
                  inherit.aes = F, linewidth = 0.5,  linetype=2) +
  # stat_density_2d(datosY, mapping=aes(x=y1, y=y2),
  #                 inherit.aes = F, linewidth = 0.8, color='red') +
  theme_minimal() + labs(
    x = expression(x[1]),
    y = expression(x[2])) + 
  scale_color_manual(name=expression(SN[2]^c*(x ~ '|'~  mu ~ ','  ~ Sigma ~ ',' ~ lambda)),
                     values=c(A='gray35', B='gray50'),
                     labels=c('Contornos', 'Estimación kernel'))+
  theme(text=element_text(family='serif'), legend.position = 'top')
x11(); print(g1)
ggsave(plot=g1, filename=r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\c-ii\SN2C-contour2.pdf)',
       width=20, height=20, units='cm')


## En 3D ====
# create plotly surface (interactive)

fuente <- list(
  family = "serif",
  size = 12,
  color = 'black')

p <- plot_ly(x = ~x1, y = ~x2, z = ~matrix(datos$fg, nrow=100, byrow=F),
             colors=colorRamp(colors = c('gray80', 'gray50')),
             # colorscale = "Greys",
             type = "mesh3d", alpha=1) %>%
  layout(
    font=fuente,
    scene = list(
      xaxis = list(title = "x1", nticks = 8),
      yaxis = list(title = "x2", nticks = 8),
      camera = list(eye = list(x = -3, y = 3, z = 0.4)),
      aspectratio = list(x = 1, y = 1, z = 1),
      zaxis = list(title = "Densidad")
    )) %>% add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      )
    )
p
length(datos$Var1)
length(datos$fg)

x11()
a <- mesh3d(x=x1, y=x2, z=matrix(datos$fg, nrow=100, byrow=F))
a

surface3d(x=x1, y=x2, z=matrix(datos$fg, nrow=100, byrow=F))

## normal simetrica
n <- 500
X <- mvtnorm::rmvnorm(n=1500, mean=c(0, 0))


g0 <- ggplot(datosX, aes(x1, x2)) +
  stat_density_2d(aes(color = after_stat(level)), linewidth = 0.8) +
  theme_minimal() + labs(
    title = "Estimated density contours",
    x = expression(x[1]),
    x1 = expression(x[2]))
x11(); print(g0)

