library(matrixStats)
library(survival)
library(glmnet)
library(rqPen)  ###rq.lasso.fit QICD
library(ECOSolveR)
library(Matrix)

###
bje_sgl <- function(x, y, index, lam1, lam2, standardize=TRUE, method=c("square", "absolute", "huber", "tukey"), k_h=1.345, k_t=4.685, maxit_em=100, eps=1e-5){

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

      fit <- linear_sgl( newx[,-1], newy, weights=newweight, index, lam1, lam2, standardize=standardize)

      beta <- c(fit$a0,fit$beta)

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta

    } else if (method == "absolute"){

      fit <- lad_socp_sgl( newx[,-1], newy, weights=newweight, index, lam1, lam2)

      beta <- c(fit$a0,fit$beta)

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta

    } else if (method == "huber"){

      fit <- hbtk_cda_sgl(X=newx[,-1],Y=newy, weights=newweight, index, lam1, lam2, method="huber",standardize=standardize)
      beta <- c(fit$a0,fit$beta)

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta


    } else if (method == "tukey"){

      fit <- hbtk_cda_sgl(X=newx[,-1],Y=newy, weights=newweight, index, lam1, lam2, method="tukey",standardize=standardize)
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
