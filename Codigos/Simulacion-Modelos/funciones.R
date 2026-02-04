

rho2lambda <- function(p){
  L <-  p / sqrt(1-p**2)
  return(L)
}

lambda2rho <- function(L){
  p <- L / sqrt(1 + L**2)
  return(p)
}

get_z <- function(y_obs){
  n <- length(y_obs)
  z_out <- c()
  for(i in 1:n){
    z_out[i] <- abs(sn::rsn(n=1, xi=0, omega=1, alpha=rho2lambda(0.5))) * sign(y_obs[i]-0.001)
  }
  return(z_out)
}


Sigma <- function(rho, N, sigma2=1.0){
  ones <- matrix(1L, nrow=N)
  J <- tcrossprod(ones)
  I <- diag(x=1L, nrow=N)
  block11 <- sigma2*rho**2*((rho**(-2) - 1)*I + J)
  block12 <- ones*rho*sigma2
  block21 <- t(block12)
  block22 <- sigma2
  MatSigma <- rbind(cbind(block11, block12), cbind(block21, block22))
  return(MatSigma)
}

rmvnorm <- function(n, mu, Sigma) {
  L <- t(chol(Sigma))  # t(Upper triangular)
  d <- length(mu)
  Z <- matrix(rnorm(d*n), nrow = d, ncol = n)
  X <- mu + (L %*% Z) # A*B = a_{i}*b^{j}, so columns are not linked to each other
  return(X)
}

rskewnorm <- function(n, mean, sigma){
  # implementación eficiente
  # p(U, W) = p(u|w) p(w),
  # so, we generate w from W, then we use it for compute p(u | w)
  k <- 0
  g <- c()
  L <- ncol(sigma)
  covar0 <- sigma[1:(L-1), L]%*%solve(sigma[L, L]) # para evitar problemas con dimensiones
  covar  <- sigma[1:(L-1), 1:(L-1)] - covar0 %*% sigma[L, 1:(L-1)]
  while(k<n){
    # w <- abs(rnorm(n=1, sd=sqrt(sigma[1, 1])))
    # w <- rnorm(n=1, sd=sqrt(sigma[1, 1]))
    w <- truncnorm::rtruncnorm(n=1, a=0, b=Inf, mean=0, sd=sqrt(sigma[1, 1]))
    # if(w>0){
      media <- mean + covar0*(w-0)
      # g <- c(g, MASS::mvrnorm(n=1, mu=media, Sigma=covar))
      #g <- c(g, rmvnorm(n=1, mu=media, Sigma=covar))
      g <- c(g, mvtnorm::rmvnorm(n=1, mean=media, sigma=covar))
      k <- k+1
    # }
  }
  return(g)
}

rskewnorm(n=1, mean=0, sigma=Sigma(rho=0.8, N=2, sigma2=1))

# a <- c()
# for(k in 1:1000){
#   a <- c(a, rskewnorm(n=1, mean=0, sigma=Sigma(rho=0.8, N=1, sigma2=1)))
# }
# fit0 <- sn::selm(a ~ 1)
# a_ <- coef(fit0, param.type = 'DP')
# lambda2rho(a_[3])
# x11(); hist(a)
# abline(v=c(mean(a), median(a)), col=c(2, 4), lwd=2)
# 
# b <- sn::rsn(n = 1000, xi=0, omega=1, alpha=rho2lambda(0.8))
# fit1 <- sn::selm(b ~ 1)
# b_ <- coef(fit1, param.type = 'DP')
# lambda2rho(b_[3])
# x11(); hist(b)
# abline(v=c(mean(b), median(b)), col=c(2, 4), lwd=2)
# 
# d <- rskewnorm(n=1, mean=0, sigma=Sigma(rho=0.8, N=1000, sigma2=1))
# fit2 <- sn::selm(d ~ 1)
# d_ <- coef(fit2, param.type = 'DP')
# lambda2rho(d_[3])
# x11(); hist(d)
# abline(v=c(mean(d), median(d)), col=c(2, 4), lwd=2)


get_idx <- function(region){
  # ¿en donde empiezan los índices de region?
  reg <- unique(region)
  pos <- c()
  for(k in reg){
    pos <- c(pos, which(region==k)[1])
  }
  # agregar la posición del último elemento
  pos <- c(pos, length(region))
  names(pos) <- reg
  return(pos)
}


get_leng <- function(region_tmp){
  leng <- NULL
  for(reg in seq_along(unique(region_tmp))){
    if(reg == length(unique(region_tmp))){
      leng[[reg]] <- get_idx(region_tmp)[reg]:length(region_tmp)
    } else{
      leng[[reg]] <- get_idx(region_tmp)[reg]:(get_idx(region_tmp)[reg+1]-1)
    }
  }
  return(leng)
}

Step <- function(a, delta=0, type=1){
  if(type == 1){
    # probit-like
    b <- ifelse(a>=delta, 1, 0)
  } else{
    print('!')
  }
  return(b)
}

get_sim <- function
(
  N=c(1500, 2000, 1000, 1500), prop=0.1, K=4, rho=c(0.5, 0.75, 0.85, 0.95), mu=c(6, 7, 6.5, 5), b=c(1.9, -0.5, -1.4, 0.0, 0.0), sigma=1.5,
  out=NULL, type = 'LogSN', center=F, noInt=F, restrict_rho=T
)
  {
  if(type%in%c('Probit', 'Ordered')){
    # mu's pequeños, de otro modo todos son cero o uno
    # a menos que se cambie el umbral...
    mu <- c(0.6, 0.7, 0.6, 0.5)
    # b <- c(1.0, -1.5, 0.8, 0.0, 0.0) # suma a 0.3
    sigma <- 1.0
  }
  if(noInt){
    mu <- mu*0
  }
  if(!restrict_rho){
    rho <- c(-0.5, 0.75, -0.85, 0.95)
  }
  n <- floor(N*prop)
  p <- length(b)
  X <- list()
  X_obs <- list()
  X_mis <- list()
  y <- c()
  y_obs <- c()
  y_mis <- c()
  region <- c()
  group <- c()
  groupmis <- c()
  idx_mis <- c()
  for(k in 1:K){
    X[[k]] <- scale(matrix(rnorm(N[k]*p), ncol=p))
    # one multivariate skew-normal draw
    # y_tmp <- rskewnorm(n=1, mean=0, sigma=Sigma(rho=rho[k], N=N[k], sigma2=sigma**2))
    # n univariate skew-normal draws
    y_tmp <- rskewnorm(n=N[k], mean=0, sigma=Sigma(rho=rho[k], N=1, sigma2=sigma**2))
    if(center){
      s2 <- sigma**2
      Exp <- sqrt(s2) * sqrt(2*rho**2 / pi)
      Var <- s2 * (1 - (2*rho**2 / pi))
      y_tmp <- mu[k] + X[[k]] %*% b + sqrt(s2)*(y_tmp-Exp[k])/sqrt(Var[k])
    } else{
      y_tmp <- mu[k] + X[[k]] %*% b + y_tmp
    }
    y <- c(y, y_tmp)
    region <- c(region, rep(k, each=N[k]))
    # idx0: indice de los observados, idx1: indice de los perdidos
    if(!is.null(out) & k%in%out){
      idx0 <- c()
      idx1 <- setdiff(1:N[k], idx0)
      #idx0 <- 0:0; idx1 <- 1:N[k];
      n[k] <- 0
    } else{
      # seleccionar con muestreo
      idx0 <- sort(sample(1:N[k], size=n[k], replace=F))
      idx1 <- setdiff(1:N[k], idx0)
      # seleccionar las primeras n[k] instancias
      # idx0 <- 1:n[k]; idx1 <- (n[k]+1):N[k]
    }
    X_obs[[k]] <- X[[k]][idx0, ]
    X_mis[[k]] <- X[[k]][idx1, ]
    y_obs <- c(y_obs, y_tmp[idx0])
    y_mis <- c(y_mis, y_tmp[idx1])
    tmp <- region[region==k]
    group <- c(group, tmp[idx0])
    groupmis <- c(groupmis, tmp[idx1])
    # almacenar idx1's
    if(k>1){
      idx1 <- idx1 + cumsum(N)[k-1]
    }
    idx_mis <- c(idx_mis, idx1)
  }
  X <- do.call(rbind, X)
  X_obs <- do.call(rbind, X_obs)
  X_mis <- do.call(rbind, X_mis)
  y <- unlist(y)
  y_obs <- unlist(y_obs)
  y_mis <- unlist(y_mis)
  output <- list(group_all=region, y_all=y, X_all=X, X_obs=X_obs, X_mis=X_mis,
                 y_obs=y_obs, y_mis=y_mis, N=N, n=n, b=b, rho=rho, sigma=sigma,
                 group=group, groupmis=groupmis, idx_mis=idx_mis,
                 p=p, K=K, mu=mu, N_obs=length(y_obs), N_mis=length(y_mis))
  return(output)
}

# 
# get_sim_probit <- function
# (
#   N=c(1500, 2000, 1000, 1500), prop=0.1, K=4, rho=c(0.5, 0.75, 0.85, 0.95), mu=c(6, 7, 6.5, 5), b=c(1.9, -0.5, -1.4, 0.0, 0.0), sigma=1.0,
#   out=NULL, center=F
# ){
#   n <- floor(N*prop)
#   p <- length(b)
#   X <- list()
#   X_obs <- list()
#   X_mis <- list()
#   y <- c()
#   y_obs <- c()
#   y_mis <- c()
#   region <- c()
#   group <- c()
#   groupmis <- c()
#   idx_mis <- c()
#   for(k in 1:K){
#     X[[k]] <- scale(matrix(rnorm(N[k]*p), ncol=p))
#     y_tmp <- rskewnorm(n=1, mean=0, sigma=Sigma(rho=rho[k], N=N[k], sigma2=sigma**2))
#     if(center){
#       s2 <- sigma**2
#       Exp <- sqrt(s2) * sqrt(2*rho**2 / pi)
#       Var <- s2 * (1 - (2*rho**2 / pi))
#       y_tmp <- mu[k] + X[[k]] %*% b + sqrt(s2)*(y_tmp-Exp[k])/sqrt(Var[k])
#     } else{
#       y_tmp <- mu[k] + X[[k]] %*% b + y_tmp
#     }
#     y_tmp <- sn::rsn(n = N[k], xi=0, omega=1, alpha=rho2lambda(rho[k]))
#     y_tmp <- c(y_tmp)
#     y_tmp
#     fit0 <- sn::selm(y_tmp ~ 1)
#     coef(fit0, param.type = 'DP')
#     x11(); plot(density(y_tmp))
#     x11(); hist(y_tmp)
#     abline(v=c(mean(y_tmp), median(y_tmp)), col=c(2, 4), lwd=2)
#     y <- c(y, y_tmp)
#     region <- c(region, rep(k, each=N[k]))
#     # idx0: indice de los observados, idx1: indice de los perdidos
#     if(!is.null(out) & k%in%out){
#       idx0 <- c()
#       idx1 <- setdiff(1:N[k], idx0)
#       #idx0 <- 0:0; idx1 <- 1:N[k];
#       n[k] <- 0
#     } else{
#       # seleccionar con muestreo
#       idx0 <- sort(sample(1:N[k], size=n[k], replace=F))
#       idx1 <- setdiff(1:N[k], idx0)
#       # seleccionar las primeras n[k] instancias
#       # idx0 <- 1:n[k]; idx1 <- (n[k]+1):N[k]
#     }
#     X_obs[[k]] <- X[[k]][idx0, ]
#     X_mis[[k]] <- X[[k]][idx1, ]
#     y_obs <- c(y_obs, y_tmp[idx0])
#     y_mis <- c(y_mis, y_tmp[idx1])
#     tmp <- region[region==k]
#     group <- c(group, tmp[idx0])
#     groupmis <- c(groupmis, tmp[idx1])
#     # almacenar idx1's
#     if(k>1){
#       idx1 <- idx1 + cumsum(N)[k-1]
#     }
#     idx_mis <- c(idx_mis, idx1)
#   }
#   X <- do.call(rbind, X)
#   X_obs <- do.call(rbind, X_obs)
#   X_mis <- do.call(rbind, X_mis)
#   y <- unlist(y)
#   y_obs <- unlist(y_obs)
#   y_mis <- unlist(y_mis)
#   output <- list(group_all=region, y_all=Step(y), X_all=X, X_obs=X_obs, X_mis=X_mis,
#                  y_obs=Step(y_obs), y_mis=Step(y_mis), N=N, n=n, b=b, rho=rho, sigma=sigma,
#                  group=group, groupmis=groupmis, idx_mis=idx_mis,
#                  p=p, K=K, mu=mu, N_obs=length(y_obs), N_mis=length(y_mis))
#   return(output)
# }

get_inits <- function(df, group=datos$group, restrict_rho=T){
  # obtiene valores iniciales con la funcion selm de sn
  r0   <- sn::selm(df$y_obs ~ ., data=df)
  y_mishat <- predict(r0, newdata = data.frame(datos$X_mis))
  r0  <- coef(r0, param.type = 'DP')
  b0 <- r0[startsWith(names(r0), 'X')]
  sigma0 <- r0['omega']
  mu0 <- numeric(datos$K)
  rho0 <- numeric(datos$K)
  for(i in unique(group)){
    dftmp <- subset(df, group==i)
    ri  <- sn::selm(dftmp$y_obs ~ ., data=dftmp)
    ri  <- coef(ri, param.type = 'DP')
    mu0[i] <- ri[1]
    if(restrict_rho){
      rho0[i] <- abs(lambda2rho(ri['alpha']))
    } else{
      rho0[i] <- lambda2rho(ri['alpha'])
    }
  }
  return(list(b=b0, rho=rho0, mu=mu0, sigma=sigma0))
}

metrics <- function(fitted, actual){
  mae <- mean(abs(actual - fitted))
  mse <- mean((actual - fitted)**2)
  corr<- cor(actual, fitted)
  mape<- 1/length(actual) * sum(abs(fitted - actual)/actual)
  output <- c(corr, mae, sqrt(mse), mape)
  names(output) <- c('Corr.', 'MAE', 'RMSE', 'MAPE')
  return(output)
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
    names(output) <- c('Acc.', 'TPR', 'TNR', 'F1-Score')
  } else{
    accuracy <- sum(diag(cm)) / total
    ktau <- cor(fitted, actual, method = "kendall") # kendall's tau
    output <- c(accuracy, ktau)
    names(output) <- c('Acc.', 'tau de Kendall')
  }
  return(output)
}

generate_design <- function(N, intercept=F){
  K <- length(N)
  cate <- c()
  for(i in 1:K){
    cate <- c(cate, rep(i, times=N[i]))
  }
  cate <- as.factor(cate)
  if(intercept){
    cate <- as.matrix(model.matrix(~cate))
  } else{
    cate <- as.matrix(model.matrix(~cate-1))
  }
  return(cate[1:dim(cate)[1], 1:dim(cate)[2]])
}

# SSVS ----

posterior_prob <- function(fit){
  posterior_inclusion <- fit$draws('m_ind')
  Mout <- apply(posterior_inclusion, 1, paste, collapse = "") # Etiquetas únicas para cada configuración
  # con un table tambien queda
  # Inicializar contenedores para los conteos y las configuraciones de modelos
  Nm <- numeric() # Conteos para cada configuración de modelo
  Mindex <- character() # Configuraciones únicas
  # Contar la ocurrencia de cada configuración única
  while (length(Mout) > 0) {
    Nm <- c(Nm, sum(Mout %in% Mout[1]))
    Mindex <- c(Mindex, Mout[1])
    Mout <- Mout[Mout != Mout[1]] # Eliminar configuraciones procesadas
  }#
  #Calcular las probabilidades posteriores
  Pm <- Nm / sum(Nm)
  # Ordenar por probabilidad posterior en orden descendente
  ord <- order(Pm, decreasing = TRUE)
  Pm <- Pm[ord]
  Mindex <- Mindex[ord]
  # Combinar en un data frame para visualizar mejor
  Mprob <- data.frame(Probabilidad_Posterior = round(Pm, 4), Configuracion_Modelo = Mindex)
  return(Mprob)
}

marginal_pi <- function(fit){
  #marginal posterior inclusion
  posterior_inclusion <- fit$draws('m_ind')
  ## each variable:
  nrow(posterior_inclusion)
  marginal <- round(colSums(posterior_inclusion) / nrow(posterior_inclusion), 4) * 100
  return(marginal)
}

get_levels <- function(thresholds, Z) {
  if ('list' %in% class(Z)){
    lapply(Z, function(Zk) {
      findInterval(Zk, vec = thresholds)
    })
  } else{
    findInterval(Z, vec = thresholds)
  } 
}

# Mapas ====

# Custom key glyph
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


