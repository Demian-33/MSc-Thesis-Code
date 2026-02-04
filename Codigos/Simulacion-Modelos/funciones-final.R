gamma1 <- function(rho){
  g1 <- (4-pi)/2 * (sqrt(2/pi) * rho )**3 / (1 - 2/pi*rho**2 )**(3/2)
  return(g1)
}

gamma1_inv <- function(g){
  c <- (4-pi)/2 * (2/pi)**(3/2)
  rho <- g**(1/3) / sqrt(c**(2/3) + 2/pi*g**(2/3))
  return(rho)
}

rho2lambda <- function(p){
  L <-  p / sqrt(1-p**2)
  return(L)
}

lambda2rho <- function(L){
  p <- L / sqrt(1 + L**2)
  return(p)
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
  covar0 <- sigma[1:L-1, L]%*%solve(sigma[L, L]) # para evitar problemas con dimensiones
  covar  <- sigma[1:L-1, 1:L-1] - covar0%*%sigma[L, 1:L-1]
  while(k<n){
    w <- rnorm(n=1)
    if(w>0){
      media <- mean + covar0*(w-0)
      #g <- c(g, MASS::mvrnorm(n=1, mu=media, Sigma=covar))
      g <- c(g, rmvnorm(n=1, mu=media, Sigma=covar))
      k <- k+1
    }
  }
  return(g)
}

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


generate_design <- function(N, intercept=TRUE){
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

generate_idx <- function(out=NULL, prop=0.001, master, region, n){
  #prop=0.01
  if(is.null(out)){
    cat('Muestreo normal: n=', n, 'prop.=', prop, '\n')
    idx <- sample(1:n, size=prop*n, replace=F)
  } else{
    # master=idx. de las obs. no nulas
    # region=numero de region de cada obs.
    cat('Muestrear área(s) completa(s):', unique(out), '\n')
    idx <- which(region[master]%in%out)
  }
  idx <- sort(idx, decreasing = F)
  # match(idx, missing) # parece que sort elimina la necesidad de hacer esto
  # a <- c(1, 2, 7); b <- sample(1:9, replace=F, size=9); a; b
  # match(a, b) # en que posición están los elementos de a en b?
  return(idx)
}

get_obs <- function(prop=NA, out=NULL, none=FALSE, censo=FALSE, X=cbind(Xcon, Xbin, Xcat), y, region){
  # none (TRUE/FALSE) : no usar ninguna parte prop del CPV 2020
  # censo (TRUE/FALSE): usar una parte prop del censo / solo perder out areas completas
  # prop (0, 1): proporcion del CPV 2020 que se quiere usar
  # out (1, 16): areas completas que se desean predecir.
  #              de otro modo, se pierden prop al azar
  master <- which(!is.na(y))
  n <- length(master)
  if(none){
    group <- as.numeric(region) # Usar CPV
    observed <- which(!is.na(y))
    missing  <- which(is.na(y))
    # Seleccionar solo una parte del CPV 2020
    if(!is.null(out) & !is.na(prop)){
      set.seed(121999)
      cat('Seleccionar', 100*prop,  '% del censo\n')
      tmp <- tibble(missing, region=region[missing])
      missing <- c()
      for(r in out){
        aux <- subset(tmp, region==r)
        missing <- c(missing, sample(aux$missing, size=nrow(aux)*prop, replace=F))
      }
      idx <- sort(missing) # para que vayan en orden
    } else{
      idx <- rep(TRUE, times=length(observed))
    }
    yobs <- y[observed]  
    ymis <- y[missing]
    Xobs <- X[observed, ]
    Xmis <- X[missing, ]
    output <- list(missing=missing, observed=observed, y_obs=yobs, y_mis=ymis,
                   X_obs=Xobs, X_mis=Xmis, group=group, idx=idx, master=master,
                   y_test=yobs, K=length(unique(region)))
    return(output)
  }
  if(censo){
    # Usar CPV 2020
    idx <- generate_idx(out=out, prop=prop, master=master, region=region, n=n) # perder areas completas
    ytest <- y[master][idx]
    y[master][idx] <- NA
    # X <- X
    group <- as.numeric(region) # Usar CPV
  } else{
    # No usar cpv 2020
    ytest <- y[master]
    unique(region[master])
    idx <- generate_idx(out=out, prop=prop, master=master, region=region, n=n)
    y <- y[master] # solo trabajar con la parte observada
    y[idx] <- NA
    ytest <- ytest[idx]
    X <- X[master, ]
    group <- as.numeric(region[master]) # no usar CPV
  }
  # datos para el modelo
  missing  <- which(is.na(y))
  observed <- which(!is.na(y))
  # Seleccionar solo una parte del CPV 2020
  if(!is.null(out) & !is.na(prop)){
    cat('Seleccionar', 100*prop, '% del censo\n')
    tmp <- tibble(missing, region=region[missing])
    missing <- c()
    for(r in out){
      aux <- subset(tmp, region==r)
      missing <- c(missing, sample(aux$missing, size=nrow(aux)*prop, replace=F))
    }
  }
  #
  yobs <- y[observed]   # or yobs <- y[-missing]
  ymis <- y[missing]
  Xobs <- X[observed, ] # o bien Xobs <- X[-missing, ]
  Xmis <- X[missing, ]
  output <- list(missing=missing, observed=observed, y_obs=yobs, y_mis=ymis, 
                 X_obs=Xobs, X_mis=Xmis, group=group, idx=idx, master=master,
                 y_test=ytest, K=length(unique(region)))
  return(output)
}


dummy_slot <- function(list0, list1, bool){
  if(bool){
    return(list0)
  } else{
    return(list1)
  }
}

# list0 <- list(letters)
# list1 <- list(c(0:9))
# dummy_slot(list0, list1, bool=TRUE)
# for(j in dummy_slot(list0, list1, bool=F)){print(j)}


get_fitted <- function(censo=FALSE, posterior='modelo', Xmis=datos$Xmis, intervalo=FALSE, max_post=250, modelo, out){
  # Existen dos formas de obtener muestras a posteriori para y_{ij} :
  # (i)  directamente de ymis en el modelo
  # (ii) con las muestras a posteriori de los parametros
  # Necesito generar
  # (1) valores ajustados, con o sin usar el censo
  # (2) pronosticos para
  # (2.1) areas completas sin usar nada del censo
  # (2.2) areas completas usando una parte del censo
  
  if(posterior=='modelo'){
    cat('Generando muestras a partir del modelo Stan. \n')
  } else{ if(posterior=='muestras'){
    cat('Generando muestras a partir de muestras a posteriori. \n')
  }
  }
  
  y_ij <- numeric()
  y_labels <- c()
  y_ij_int <- c() # Intervalo =  TRUE
  idx_fitted <- c()
  if(any(is.na(out))){
    out <- unique(region[datos$master][datos$idx])
  }
  
  
  if(censo){
    tmp <- get_leng(region) # Usar CPV, tarda un poco
  } else{
    tmp <- get_leng(region[datos$master][datos$idx]) # No usar CPV
  }
  
  if(posterior!='modelo'){
    # rstan::extract(modelo)$mu[1:max_post, ]
    mu    <- fit_vb$draws("mu")[1:max_post, ]
    beta  <- fit_vb$draws("b")[1:max_post, ]
    sigma <- fit_vb$draws("sigma")[1:max_post]
    rho   <- fit_vb$draws("rho")[1:max_post, ]
    w <- fit_vb$draws("z")[1:max_post, ]
    X <- as.matrix(Xmis)
    set.seed(1)
    # Generar muestras de la a posteriori
    for(k in 1:length(out)){
      # j: indexa cada observacion
      # k: indexa regiones
      # l: indexa renglon de la a posteriori
      # p: indexa covariables
      for(j in dummy_slot(tmp[[out[k]]], tmp[[k]], bool=censo)){
        # tmp[[out[k]]]: usar CPV; tmp[[k]]: no usar CPV.
        cat('Generar observacion y_{', out[k], ',' ,j, '}\n', sep='')
        # Reiniciar para cada observacion nueva
        eta <- NULL
        nu  <- NULL
        for(l in 1:max_post){
          # l es el numero de repeticion, renglon o muestra VB
          # out[k] selecciona el numero de covariable necesaria
          fun_rho0 <- sqrt(2/pi*rho[l, out[k]]**2)
          fun_rho1 <- sqrt(1 - 2/pi*rho[l, out[k]]**2)
          eta[l] <- mu[l, out[k]] + tcrossprod(X[j, ], beta[l, ]) + w[l, out[k]]*rho[l, out[k]]
          eta[l] <- eta[l] - sigma[l]*fun_rho0/ fun_rho1
          nu[l]  <- sigma[l]*sqrt(1-rho[l, out[k]]**2) / fun_rho1
        }
        # Promediar sobre cada renglon de la a posteriori
        samples <- rnorm(n=max_post, mean=eta, sd=nu)
        y_ij <- c(y_ij, mean(samples))
        y_ij_int <- c(y_ij_int, quantile(samples, c(0.025, 0.5, 0.975)))
        y_labels <- c(y_labels, j)
      }
    }
    
    if(censo){
      # idx_fitted <- match((1:length(ym))[datos$master][datos$idx], y_labels)
      # fitted_vb <- y_ij[idx_fitted] # Usar CPV
      idx_fitted <- match(datos$master[datos$idx], y_labels)
      fitted_vb <- y_ij[idx_fitted] # Usar CPV
    } else{
      fitted_vb <- y_ij # No usar CPV
    }
    
  } else{
    # Alternativa (1): ajuste del modelo Stan
    for(k in 1:length(out)){
      for(j in dummy_slot(tmp[[out[k]]], tmp[[k]], bool=censo)){
        y_labels <- c(y_labels, j)
      }
    }
    idx_fitted <- match(datos$master[datos$idx], y_labels)
    # apply(rstan::extract(modelo)$ymis, 2, mean)
    fitted_vb <- apply(modelo$draws("ymis")[1:max_post, ], 2, mean)
  }
  
  output <- list(
    fitted_vb_y_ij = y_ij, # todos los ajustados, por ejemplo, Censo
    fitted_vb_int = y_ij_int, # intervalos de todos los ajustados
    fitted_vb=fitted_vb, # ajustados para los datos observados en out
    # metricas0=metrics(fitted_vb, ytest), # metricas
    lista_reg=tmp, # indice de cada observación por region
    indices=idx_fitted
  )
  return(output)
}

get_fitted_samples <- function(
    region=region,
    censo=censo, intervalo=FALSE, which_out=unique(region[datos$master][datos$idx]), max_post=150, X){
  
  X <- as.matrix(X)
  y_ij <- numeric()
  y_labels <- NULL
  y_ij_int <- c() # Intervalo =  TRUE
  idx_fitted <- NA
  
  if(censo){
    tmp <- get_leng(region) # Usar CPV, tarda un poco
  } else{
    tmp <- get_leng(region[datos$master][datos$idx]) # No usar CPV
  }
  
  mu    <- fit_vb$draws("mu")[1:max_post, ]
  beta  <- fit_vb$draws("b")[1:max_post, ]
  sigma <- fit_vb$draws("sigma")[1:max_post]
  rho   <- fit_vb$draws("rho")[1:max_post, ]
  w <- fit_vb$draws("z")[1:max_post, ]
  X <- as.matrix(X)
  set.seed(1)
  
  # Generar muestras de la a posteriori
  for(k in 1:length(which_out)){
    # j: indexa cada observacion
    # k: indexa regiones
    # l: indexa renglon de la a posteriori
    # p: indexa covariables
    # tmp[[out[k]]]: usar CPV; tmp[[k]]: no usar CPV.
    cat('Generando observaciones de la region: ', which_out[k], '\n', sep='')
    for(j in dummy_slot(tmp[[which_out[k]]], tmp[[k]], bool=censo)){
      # cat('Generar observacion y_{', which_out[k], ',' ,j, '}\n', sep='')
      # Reiniciar para cada observacion nueva
      eta <- NULL
      nu  <- NULL
      for(l in 1:max_post){
        # l es el numero de repeticion, renglon o muestra VB
        # out[k] selecciona el numero de covariable necesaria
        fun_rho0 <- sqrt(2/pi*rho[l, which_out[k]]**2)
        fun_rho1 <- sqrt(1 - 2/pi*rho[l, which_out[k]]**2)
        eta[l] <- mu[l, which_out[k]] + tcrossprod(X[j, ], beta[l, ]) + w[l, which_out[k]]*rho[l, which_out[k]]
        eta[l] <- eta[l] - sigma[l]*fun_rho0/ fun_rho1
        nu[l]  <- sigma[l]*sqrt(1-rho[l, which_out[k]]**2) / fun_rho1
      }
      # Promediar sobre cada renglon de la a posteriori
      samples <- rnorm(n=max_post, mean=eta, sd=nu)
      y_ij <- c(y_ij, mean(samples))
      if(intervalo){
        y_ij_int <- c(y_ij_int, quantile(samples, c(0.025, 0.5, 0.975)))
      }
      y_labels <- c(y_labels, j)
    }
  }
  
  # Devuelve los pronosticos asociados a los datos de prueba
  # if(censo){
  #   # idx_fitted <- match((1:length(ym))[datos$master][datos$idx], y_labels)
  #   # fitted_vb <- y_ij[idx_fitted] # Usar CPV
  #   idx_fitted <- match(datos$master[datos$idx], y_labels)
  #   fitted_vb <- y_ij[idx_fitted] # Usar CPV
  # } else{
  #   fitted_vb <- y_ij # No usar CPV
  # }
  
  output <- list(
    fitted_vb_y_ij = y_ij, # todos los ajustados, por ejemplo, Censo
    fitted_vb_int = y_ij_int, # intervalos de todos los ajustados
    lista_reg=tmp, # indice de cada observación por region
    indices=datos$master[datos$idx]
  )
  return(output)
}

get_fit <- function(object, type='', rem=F){
  out <- object$fitted_vb_int
  Lout <- length(out)
  if(type == 'Lower'){
    out <- out[seq(1, Lout, by=3)]
  } else{
    if(type == 'Upper'){
      out <- out[seq(3, Lout, by=3)]
    } else{
      out <- object$fitted_vb_y_ij
    }
  }
  if(rem){
    idx <- which(!is.na(ym))
    out <- ym[idx]
  }
  return(out)
}

get_fitted_cat <- function(censo=FALSE, posterior='modelo', ytest=datos$ytest, Xmis=datos$Xmis, intervalo=FALSE, max_post, modelo){
  # Existen dos formas de obtener muestras a posteriori para z_{ij} :
  # (1) directamente de zmis en el modelo
  # (2) con las muestras a posteriori de los parametros
  # En general, (1) solo se puede usar con pocos datos faltantes
  if(censo==TRUE & posterior=='modelo'){
    cat('Muy lento!. Establece posterior == "muestras"\n')
    return(0)
  }
  if(posterior!='modelo'){
    # Es decir, vamos a usar la alternativa (2)
    mu    <- rstan::extract(modelo)$mu[1:max_post, ]
    beta  <- rstan::extract(modelo)$b[1:max_post, ]
    sigma <- rep(1.0, times=max_post)
    rho   <- rstan::extract(modelo)$rho[1:max_post, ]
    w <- rstan::extract(modelo)$w[1:max_post, ]
    z_ij <- numeric()
    Xmis <- as.matrix(Xmis)
    set.seed(1)
    z_ij <- NULL
    z_labels <- NULL
    if(intervalo){
      z_ij_int <- c()
    }
    if(censo){
      tmp <- get_leng(region) # Usar CPV, tarda un poco
    } else{
      tmp <- get_leng(region[datos$master][datos$idx]) # No usar CPV
    }
    
    # Generar muestras de la a posteriori
    for(k in 1:length(which_out)){
      # j: indexa cada observacion
      # k: indexa regiones
      # l: indexa renglon de la a posteriori
      # p: indexa covariables
      # tmp[[out[k]]]: usar CPV; tmp[[k]]: no usar CPV.
      for(j in dummy_slot(tmp[[which_out[k]]], tmp[[k]], bool=censo)){
        cat('Generar observacion y_{', which_out[k], ',' ,j, '}\n', sep='')
        # Reiniciar para cada observacion nueva
        eta <- NULL
        nu  <- NULL
        for(l in 1:max_post){
          # l es el numero de repeticion, renglon o muestra VB
          # out[k] selecciona el numero de covariable necesaria
          fun_rho0 <- sqrt(2/pi*rho[l, out[k]]**2)
          fun_rho1 <- sqrt(1 - 2/pi*rho[l, out[k]]**2)
          eta[l] <- mu[l, out[k]] + Xmis[j, ] %*% beta[l, ] + w[l, out[k]]*rho[l, out[k]]
          eta[l] <- eta[l] - sigma[l]*fun_rho0/ fun_rho1
          nu[l]  <- sigma[l]*sqrt(1-rho[l, out[k]]**2) / fun_rho1
        }
        # Promediar sobre cada renglon de la a posteriori
        z_ij <- c(z_ij, mean(rnorm(n=max_post, mean=eta, sd=nu)))
        z_labels <- c(z_labels, j)
        z_ij_int <- c(z_ij_int, quantile(rnorm(n=max_post, mean=eta, sd=nu), c(0.025, 0.5, 0.975)))
      }
    }
    
    y_ij <- ifelse(z_ij>0, 1L, 0L)
    
    if(censo){
      idx_fitted <- match((1:length(ym))[datos$master][datos$idx], z_labels)
      # fitted_vb <- y_ij[match(master[idx], z_labels)] # Usar CPV
      fitted_vb <- y_ij[idx_fitted] # Usar CPV
    } else{
      fitted_vb <- y_ij # No usar CPV
    }
    
  } else{
    # Alternativa (1): ajuste del modelo Stan
    fitted_vb <- ifelse(apply(rstan::extract(modelo)$zmis, 2, mean)>0, 1, 0)
  }
  
  output <- list(
    fitted_vb_z_ij = z_ij,
    fitted_vb_int = z_ij_int,
    fitted_vb=fitted_vb,
    metricas0=metrics_cat(fitted_vb, ytest),
    metricas1=caret::confusionMatrix(table(fitted_vb, ytest)),
    indices=tmp)
  return(output)
}

# df <- tmp
# colnames(tmp)
# group <- datos$group[datos$master]

get_inits <- function(df, group, single_rho=F, single_mu=F, which_in =c(1:13), only_Int=F){
  # obtiene valores iniciales con la funcion selm de sn
  r0   <- sn::selm(df$y_obs ~ ., data=df)
  # y_mishat <- predict(r0, newdata = data.frame(X_mis))
  r0  <- coef(r0, param.type = 'DP')
  b0 <- r0[startsWith(names(r0), 'X')]
  sigma0 <- r0['omega']
  mu0 <- numeric(datos$K)
  rho0 <- numeric(datos$K)
  for(i in unique(group)){
    dftmp <- cbind(df, group=group)
    dftmp <- subset(dftmp, group==i)
    dftmp$group <- NULL
    dftmp <- dftmp[, c(1, 1+which_in)]
    if(only_Int){
      ri  <- sn::selm(dftmp$y_obs ~ 1, data=dftmp)
    } else{
      ri  <- sn::selm(dftmp$y_obs ~ ., data=dftmp)
    }
    ri  <- coef(ri, param.type = 'DP')
    mu0[i] <- ri[1]
    rho0[i] <- abs(lambda2rho(ri['alpha']))
  }
  if(single_rho){
    rho0 <- abs(lambda2rho(r0['alpha']))
  }
  if(single_mu){
    mu0 <- r0[1]
  }
  return(list(b=b0, rho=rho0, mu=mu0, sigma=sigma0))
}

get_coefs <- function(a, intercept=FALSE){
  tmp <- coef(a)
  tmp[which(is.na(tmp))] <- 0.001
  if(!intercept){
    tmp <- tmp[-1]
  }
  return(tmp)
}


# SSVS ----

delta <- function(tau2, c2){
  d <- sqrt(tau2) * sqrt(2*c2*log(sqrt(c2)) / (c2 - 1))
  return(d)
}


posterior_prob <- function(VB_fit){
  # posterior_inclusion <- rstan::extract(VB_fit, "m_ind")$m_ind
  posterior_inclusion <- VB_fit$draws("m_ind")
  Mout <- apply(posterior_inclusion, 1, paste, collapse = "") # Etiquetas únicas para cada configuración
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
  # posterior_inclusion <- rstan::extract(fit, "m_ind")$m_ind
  posterior_inclusion <- fit$draws("m_ind")
  ## each variable:
  nrow(posterior_inclusion)
  marginal <- round(colSums(posterior_inclusion) / nrow(posterior_inclusion), 3) * 100
  return(marginal)
}


