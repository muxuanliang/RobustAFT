library(matrixStats)
library(survival)
library(glmnet)
library(rqPen)  ###rq.lasso.fit QICD

huber_cda <- function(X,Y,lambda,weights,method=c("huber", "tukey"),maxit=10^3,eps=1e-5,k_h=1.345, k_t=4.685,standardize=TRUE){

  method <- match.arg(method)

  n <- dim(X)[1]
  p <- dim(X)[2]

  x <- X
  y <- Y
  w <- weights/(sum(weights))

  if (standardize == TRUE){
    meanx = apply(x, 2, function(t){weighted.mean(t,w)})
    x = scale(x, meanx, FALSE)
    normx = apply(x, 2, function(t){sqrt(sum(w*t^2))})
    x = scale(x, FALSE, normx)
    meany = weighted.mean(y,w)
  }

  if (standardize == FALSE){
    meanx = 0
    normx = 1
  }

  iter <- 0
  dif <- 1

  fit <- rq.lasso.fit(x,y,lambda=0,tau=0.5,intercept=T,weights=w,method="br")
  beta <- fit$coefficients[-1]
  a0 <- fit$coefficients[1]

  w <- weights

  while (iter <= maxit & dif >= eps){

    beta_old <- beta
    a0_old <- a0

    for(j in 1:p){
      r <- (y - a0 - x%*%as.matrix(beta))
      sigma <- weightedMad(r,w)
      r <- r/sigma

      if (method == "huber"){
        beta[j] <- St(2*sum(w*x[,j]^2)*beta_old[j]/sigma^2+sum(w*d.hb(r,k_h)*x[,j])/sigma,lambda/sigma)/(2*sum(w*x[,j]^2)/sigma^2)
      } else if (method == "tukey") {
        beta[j] <- St(2*sum(w*x[,j]^2)*beta_old[j]/sigma^2+sum(w*d.tk(r,k_t)*x[,j])/sigma,lambda/sigma)/(2*sum(w*x[,j]^2)/sigma^2)
      }

    }

    r <- y - a0 - x%*%as.matrix(beta)
    sigma <- weightedMad(r,w)
    r <- r/sigma

    if (method == "huber"){
      a0 <- a0_old + sum(w*d.hb(r,k_h))/(2*sum(w))*sigma
    } else if (method == "tukey") {
      a0 <- a0_old + sum(w*d.tk(r,k_t))/(2*sum(w))*sigma
    }

    iter <- iter+1
    dif <- max(abs(beta-beta_old),abs(a0-a0_old))
  }

  beta <- beta/normx
  a0 <- a0 - sum(beta*meanx)

  return(list(a0=a0, beta=beta, iter=iter, dif=dif))

}


###
bje_ly <- function(x, y, method=c("square", "absolute", "huber", "tukey", "quantile_huber"), k_h=1.345, k_t=4.685, lambda, standardize=TRUE, maxit_em=100, maxit_irls=1000, eps=1e-5){

  method <- match.arg(method)

  nobs <- NROW(x)

  x <- cbind(matrix(1,nobs,1), x)

  p <- NCOL(x)

  beta_pre <- matrix(0, p, 1)
  beta_store <- beta_pre

  beta_local <- NULL

  iter <- 1
  dif <- 1

  while (iter < maxit_em & dif >= eps){

    iter <- iter + 1

    eta <- y[,1] - x %*% beta_pre
    order.eta <- order(eta)
    eta.order <- eta[order.eta]
    y.order <- y[order.eta,1]
    state.order <- y[order.eta,2]
    x.order <- x[order.eta,]

    state.order[nobs] <- 1

    newx <- x.order[state.order==1,]
    newy <- y.order[state.order==1]
    newweight <- newy * 0 + 1

    eta.surv <- Surv(1:length(eta.order), state.order)
    km.eta <- survfit(eta.surv~1)
    weight.low <- km.eta$surv
    weight.up <- diff(-c(1, km.eta$surv))

    ################# data augmentation ###############
    repeats <- rev(cumsum(rev(state.order)))[state.order==0]
    newx.add <- apply(x.order[state.order==0,,drop = FALSE],2,function(t){rep(t,repeats)})
    y.order.list <- as.list(as.data.frame(t(cbind((1:nobs)[state.order==0], repeats))))
    addy <- function(t){
      res <- rep(x.order[t[1],,drop = FALSE] %*% beta_pre, t[2]) + eta.order[((t[1]+1):nobs)[state.order[(t[1]+1):nobs]==1]]
      res
    }
    newy.list.add <- lapply(y.order.list,addy)
    newy.add <- unlist(newy.list.add, use.names = FALSE)
    addw <- function(t){
      res <- weight.up[(t[1]+1):nobs]/weight.low[t[1]]
      res[res>0]
    }
    newweight.list.add <- lapply(y.order.list, addw)
    newweight.add <- unlist(newweight.list.add, use.names = FALSE)
    newx <- rbind(newx, newx.add)
    newy <- c(newy, newy.add)
    newweight <- c(newweight, newweight.add)

    if (method == "square"){

      fit <- glmnet(newx[,-1], newy, family="gaussian", weights=newweight, alpha=1, lambda=lambda/2/sum(newweight), standardize=standardize, intercept=T)
      beta <- c(fit$a0, as.numeric(fit$beta))

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta

    } else if (method == "absolute"){

      fit <- rq.lasso.fit(newx[,-1], newy, lambda=lambda/2/length(newweight), tau=0.5,intercept=T, weights=newweight,method="br")
      beta <- fit$coefficients

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta

    } else if (method == "huber"){

      fit <- huber_cda(X=newx[,-1],Y=newy,lambda=lambda,weights=newweight,method="huber",standardize=standardize)
      beta <- c(fit$a0,fit$beta)

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta

    } else if (method == "tukey"){

      fit <- huber_cda(X=newx[,-1],Y=newy,lambda=lambda,weights=newweight,method="tukey",standardize=standardize)
      beta <- c(fit$a0,fit$beta)

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta

    }

    print(iter)
    print(dif)
  }

  if (dif < eps){
    loc <- which.min(apply(beta_store[,1:(iter-1)], 2, function(t){max(abs(beta - t))}))
    beta_local <- beta_store[,loc:iter]

  } else { beta_local <- beta }


  if (NCOL(beta_local) > 1){
    beta_final <- beta_local[,NCOL(beta_local)]
    beta_final_mean <- rowMeans(beta_local)
  } else { beta_final <- beta_local
  beta_final_mean <- beta_local }

  select <- (abs(beta_final)>eps)
  #print(select)
  #if (lambda!=0 & sum(select[-1])>1){
  #  x.nointercept <- x[,-1]
  #  refit <- bje_ly(x.nointercept[,select[-1]], y, method, k_h, k_t, lambda=0, standardize = FALSE)
  #  beta_refit <- 0 * beta_final
  #  beta_refit[select] <- refit$beta
  #} else {
  #  beta_refit <- beta_final
  #}

  if (sum(select[-1])>1){
    x.nointercept <- x[,-1]
    refit <- bje_refit(x.nointercept[,select[-1]], y, method, k_h=k_h, k_t=k_t, maxit_em=maxit_em, eps=eps)
    beta_refit <- 0 * beta_final
    beta_refit[select] <- refit$beta
    beta_refit_mean <- 0 * beta_final
    beta_refit_mean[select] <- refit$beta_mean
  } else {
    beta_refit <- beta_final
    beta_refit_mean <- beta_final
  }

  select_beta <- (abs(beta_final[-1])>eps)

  return(list(beta=beta_final,beta_mean=beta_final_mean,beta_local=beta_local,
              beta_refit=beta_refit,beta_refit_mean=beta_refit_mean,
              select=select_beta,iter=iter))

}
