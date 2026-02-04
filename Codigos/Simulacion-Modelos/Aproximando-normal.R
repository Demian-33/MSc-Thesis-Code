# Ignorar
# library(dplyr)
# x <- runif(200, -3, 3)
# y <- x + runif(200, -3, 3)
# datos <- tibble(x, y)
# x11(); plot(x, y)
# 
# x11()
# ggplot(datos, aes(x = x, y = y)) +
#   geom_jitter() + # Add jittered points
#   geom_path()    # Add a path connecting the points in the order they appear in the data
# 
# x11()
# ggplot(datos, aes(x = x, y = y)) +
#   geom_point() + # draw points, geom_line() los ordena
#   geom_path(position=position_jitter(w=0.02, h=0))    # Add a path connecting the points in the order they appear in the data
# 
# Aproximando una distribución gaussiana en dos dimensiones ====
rm(list=ls())

library(ggplot2)
library(dplyr)

mu0 <- c(-1, 1)
Sigma0 <- matrix(c(1.0, 0.7, 0.7, 2), ncol=2, byrow=T)

# (1) Usando Gibbs ====

complete_cond <- function(mu=mu0, Sigma=Sigma0, cond0=1, cond1=2, xcond=0, n){
  cond_tmp  <- Sigma[cond0, cond1]%*%solve(Sigma[cond1, cond1])
  cond_mean <- mu[cond0] + cond_tmp * (xcond - mu[cond1])
  cond_cova <- cond_tmp %*% Sigma[cond1, cond0]
  xcond0 <- rnorm(n=n, mean=cond_mean, sd=sqrt(cond_cova))
  return(xcond0)
}

set.seed(1)
iter <- 200
x1 <- c(-2.0)
x2 <- c(2.5)
for(k in 1:iter){
  tmp <- complete_cond(cond0=1, cond1=2, xcond=x2[k], n=1) 
  x1 <- c(x1, tmp)
  tmp <- complete_cond(cond0=2, cond1=1, xcond=x1[k], n=1)
  x2 <- c(x2, tmp)
}


df_gibbs <- tibble(x1, x2) %>%
  mutate(likehood=mvtnorm::dmvnorm(., mean=mu0, sigma=Sigma0))
data.grid <- expand.grid(s.1 = seq(-4, 4, length.out=300), s.2 = seq(-4, 4, length.out=300))
q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu0, sigma = Sigma0))
gibbs <- ggplot(q.samp, aes(x=s.1, y=s.2, z=prob)) + 
  geom_contour(bins=15) + # aes(colour = stat(level)), 
  coord_fixed(xlim = c(-3, 1), ylim = c(-1, 3), ratio = 1) + 
  # geom_hex(data=df_gibbs, aes(x=x1, y=x2), inherit.aes = F)+
  # geom_path(data=df_gibbs, aes(x=x1, y=x2, color=likehood), linewidth=1.2, inherit.aes = FALSE,
  #           arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) + 
  # scale_color_gradientn(colors=c("#E5F5E0", "#A1D99B", "#31A354")) +
  geom_point(data=df_gibbs, aes(x=x1, y=x2, size=likehood, color=likehood), inherit.aes = FALSE)+
  geom_path(data=df_gibbs, aes(x=x1, y=x2), linewidth=0.5, alpha=0.75, color='orange', inherit.aes = FALSE,
            arrow = arrow(ends = "last", type = "closed", length = unit(0.4, "cm"))) +
  scale_color_gradientn(colors=grey.colors(n=8, rev=T)) +
  labs(x=expression(x[1]), y=expression(x[2])) + 
  theme_minimal() + 
  theme(aspect.ratio = 1,    text = element_text(family = "serif", size=11),
        legend.position = 'none',
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # plot.background = element_rect(colour = "black", fill = NA, linewidth = 1.2)
  )
x11(); print(gibbs)

# (2) Usando Metropolis con caminata aleatoria ====

set.seed(1)
iter <- 200
x1 <- c(-2.0)
x2 <- c(2.5)
set.seed(1)
for(k in 1:iter){
  xt <- c(x1[k], x2[k])
  # xstar <- mvtnorm::rmvnorm(n=1, mean=c(0, 0))
  xstar <- mvtnorm::rmvnorm(n=1, mean=xt) # candidato
  pi_xt <- mvtnorm::dmvnorm(x=xt, mean=mu0, sigma=Sigma0) # likelihood actual
  pi_star <- mvtnorm::dmvnorm(x=xstar, mean=mu0, sigma=Sigma0) # likelihood candidato
  alpha <- min(c(1, pi_star/pi_xt))
  u <- runif(n=1, 0, 1)
  if(u < alpha){
    x1[k+1] <- xstar[1]
    x2[k+1] <- xstar[2]
  } else{
    x1[k+1] <- xt[1]
    x2[k+1] <- xt[2]
  }
}

df_metropolis <- tibble(x1, x2) %>%
  mutate(likehood=mvtnorm::dmvnorm(., mean=mu0, sigma=Sigma0))
data.grid <- expand.grid(s.1 = seq(-4, 4, length.out=300), s.2 = seq(-4, 4, length.out=300))
q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu0, sigma = Sigma0))
metropolis <- ggplot(q.samp, aes(x=s.1, y=s.2, z=prob)) + 
  geom_contour(bins=15) + # aes(colour = stat(level)), 
  coord_fixed(xlim = c(-3, 1), ylim = c(-1, 3), ratio = 1) + 
  # geom_hex(data=df_gibbs, aes(x=x1, y=x2), inherit.aes = F)+
  # geom_path(data=df_gibbs, aes(x=x1, y=x2, color=likehood), linewidth=1.2, inherit.aes = FALSE,
  #           arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) + 
  # scale_color_gradientn(colors=c("#E5F5E0", "#A1D99B", "#31A354")) +
  geom_point(data=df_metropolis, aes(x=x1, y=x2, size=likehood, color=likehood), inherit.aes = FALSE)+
  geom_path(data=df_metropolis, aes(x=x1, y=x2), linewidth=0.5, alpha=0.75, color='orange', inherit.aes = FALSE,
            arrow = arrow(ends = "last", type = "closed", length = unit(0.4, "cm"))) +
  scale_color_gradientn(colors=grey.colors(n=8, rev=T)) +
  labs(x=expression(x[1]), y=expression(x[2])) + 
  theme_minimal() + 
  theme(aspect.ratio = 1,    text = element_text(family = "serif", size=11),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = 'none',
        # plot.background = element_rect(colour = "black", fill = NA, linewidth = 1.2)
  )
x11(); print(metropolis)

# (3) Usando VB Mean Field ====

q_opt <- function(x, mu=mu0, Sigma=Sigma0, E_x=mu0){
  mu1 <- mu[1] + Sigma[1, 2]%*%solve(Sigma[2, 2]) * (E_x[2] - mu[2])
  mu2 <- mu[2] + Sigma[2, 1]%*%solve(Sigma[1, 1]) * (E_x[1] - mu[1])
  Sigma1 <- Sigma[1, 1] - Sigma[1, 2]%*%solve(Sigma[2, 2])%*%Sigma[2, 1]
  Sigma2 <- Sigma[2, 2] - Sigma[2, 1]%*%solve(Sigma[1, 1])%*%Sigma[1, 2]
  q1 <- dnorm(x[1], mean=mu1, sd=sqrt(Sigma1))
  q2 <- dnorm(x[2], mean=mu2, sd=sqrt(Sigma2))
  q  <- q1 * q2
  return(q)
}

x0 <- seq(-4, 4, length=200)
data.grid2 <- expand.grid(x0, x0)
colnames(data.grid2) <- c('x1', 'x2')
q.samp2 <- cbind(data.grid2, q=apply(data.grid2, 1, q_opt))
head(q.samp2)

meanfield <- ggplot(q.samp, aes(x=s.1, y=s.2, z=prob), color='red') + 
  scale_color_manual(values='red', label='XD') + 
  geom_contour(bins=15) + # aes(colour = stat(level)) + 
  geom_contour(data=q.samp2, mapping=aes(x=x1, y=x2, z=q),
               bins=7, col='orange', linewidth=1.2) + 
  coord_fixed(xlim = c(-3, 1), ylim = c(-1, 3), ratio=1, clip = 'on') +
  labs(x=expression(x[1]), y=expression(x[2])) + 
  theme_minimal() + 
  theme(aspect.ratio = 1,    text = element_text(family = "serif", size=11),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # plot.background = element_rect(colour = "black", fill = NA, linewidth = 1.2)
        )
x11(); print(meanfield)


# (4) Usando HMC ====

gaussian_U <- function(x, mu=mu0, Sigma=Sigma0){
  log_prop <- 0.5 * t(x-mu)%*%solve(Sigma)%*%(x-mu) # el menos lo pone HMC
  # log_prop <- mvtnorm::dmvnorm(x, mean=mu, sigma=Sigma, log=T)
  return(as.numeric(log_prop))
}

grad_gaussian_U <- function(x, mu=mu0, Sigma=Sigma0){
  grad <- solve(Sigma) %*% (x-mu0) # el menos lo pone HMC
  return(grad)
}


HMC = function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q) # accept
  }
  else
  {
    return (current_q) # reject
  }
}

set.seed(1)
iter <- 20
samples <- matrix(NA, iter-1, 2)
q0 <- c(-2, 2.5)
samples[1, ] <- q0
for(i in 2:(iter-1)){
  q <- HMC(U=gaussian_U, grad_U = grad_gaussian_U,
           epsilon = 0.1, L = 10, current_q = q0)
  samples[i, ] <- q
  q0 <- q
}

colMeans(samples)
cov(samples)
Sigma0

df_hmc <- tibble(x1=samples[, 1], x2=samples[, 2]) %>%
  mutate(likehood=mvtnorm::dmvnorm(., mean=mu0, sigma=Sigma0))
data.grid <- expand.grid(s.1 = seq(-4, 4, length.out=300), s.2 = seq(-4, 4, length.out=300))
q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu0, sigma = Sigma0))
hmc <- ggplot(q.samp, aes(x=s.1, y=s.2, z=prob)) + 
  geom_contour(bins=15) + # aes(colour = stat(level)), 
  coord_fixed(xlim = c(-3, 1), ylim = c(-1, 3), ratio = 1) + 
  geom_point(data=df_hmc, aes(x=x1, y=x2, size=likehood, color=likehood), inherit.aes = FALSE)+
  geom_path(data=df_hmc, aes(x=x1, y=x2), linewidth=0.5, alpha=0.75, color='orange', inherit.aes = FALSE,
            arrow = arrow(ends = "last", type = "closed", length = unit(0.4, "cm"))) +
  # geom_path(data=df_hmc, aes(x=x1, y=x2, color=likehood), linewidth=1.2, inherit.aes = FALSE,
  #           arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) +
  # la desventaja de esto, es que las lineas pueden empezar en baja prob. y terminar en alta
  # y tener el mismo color
  scale_color_gradientn(colors=gray.colors(n=2, rev = T))+
  labs(x=expression(x[1]), y=expression(x[2])) + 
  theme_minimal() + 
  theme(aspect.ratio = 1,    text = element_text(family = "serif", size=11),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = 'none',
        # plot.background = element_rect(colour = "black", fill = NA, linewidth = 1.2)
  )
x11(); print(hmc)

# Guardar ====

path <- r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras)'
setwd(path)

ggsave(plot = gibbs, file='./c-iii/normal-gibbs.pdf',
       height = 15, width=15, dpi=300, units = 'cm')

ggsave(plot = gibbs, file='./c-iii/normal-gibbs-after.pdf',
       height = 15, width=15, dpi=300, units = 'cm')

ggsave(plot = metropolis, file='./c-iii/normal-metropolis.pdf',
       height = 15, width=15, dpi=300, units = 'cm')

ggsave(plot = metropolis, file='./c-iii/normal-metropolis-after.pdf',
       height = 15, width=15, dpi=300, units = 'cm')

ggsave(plot = hmc, file='./c-iii/normal-hmc.pdf',
       height = 15, width=15, dpi=300, units = 'cm')

ggsave(plot = hmc, file='./c-iii/normal-hmc-after.pdf',
       height = 15, width=15, dpi=300, units = 'cm')

ggsave(plot = meanfield, file='./c-iv/normal-meanfield.pdf',
       height = 15, width=15, dpi=300, units = 'cm')



