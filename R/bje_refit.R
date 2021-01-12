bje_refit <- function(x, y, method=c("square", "absolute", "huber", "tukey"), k_h=1.345, k_t=4.685, maxit_em=100, maxit_irls=10^3, eps=1e-5){

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

    eta.surv <- survival::Surv(1:length(eta.order), state.order)
    km.eta <- survival::survfit(eta.surv~1)
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

      fit <- lm(newy~newx[,-1],weights=newweight)
      beta <- fit$coefficients

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta

    } else if (method == "absolute"){

      fit <- rqPen::rq.lasso.fit(newx[,-1], newy, lambda=0, tau=0.5,intercept=T, weights=newweight, method="br")
      beta <- fit$coefficients

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta

    } else if (method == "huber"){

      beta_irls_pre <- beta_pre
      dif.irls <- 1
      iter.irls <- 0
      while (iter.irls <= maxit_irls & dif.irls >= eps){

        residuals <- newy-newx %*% beta_irls_pre
        residuals <- residuals/weightedMad(residuals,w=newweight)

        w <- w.hb_all(residuals,k=k_h) * newweight

        fit <- lm(newy~newx[,-1],weights=w)
        beta_irls <- fit$coefficients

        dif.irls <- max(abs(beta_irls - beta_irls_pre))
        beta_irls_pre <- beta_irls
        iter.irls <- iter.irls + 1
      }
      beta <- beta_irls

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta
      #print(iter.irls)

    } else if (method == "tukey"){

      beta_irls_pre <- beta_pre
      dif.irls <- 1
      iter.irls <- 0
      while (iter.irls <= maxit_irls & dif.irls >= eps){

        residuals <- newy-newx %*% beta_irls_pre
        residuals <- residuals/weightedMad(residuals,w=newweight)

        w <- w.tk_all(residuals,k=k_t) * newweight

        fit <- lm(newy~newx[,-1],weights=w)
        beta_irls <- fit$coefficients

        dif.irls <- max(abs(beta_irls - beta_irls_pre))
        beta_irls_pre <- beta_irls
        iter.irls <- iter.irls + 1
      }
      beta <- beta_irls

      dif <- min(apply(beta_store, 2, function(t){max(abs(beta - t))}))

      beta_store <- cbind(beta_store,beta)

      beta_pre <- beta
      #print(iter.irls)

    }

    #print(iter)
    #print(dif)
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

  select_beta <- (abs(beta_final[-1])>eps)

  return(list(beta=beta_final,beta_mean=beta_final_mean,beta_local=beta_local,select=select_beta,iter=iter))

}
