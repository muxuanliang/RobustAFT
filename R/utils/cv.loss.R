cv.loss <- function(surv, predict.x, fit,method=c("square", "absolute", "huber", "tukey"), k_h=1.345, k_t=4.685, sigma=1){

  nobs <- dim(predict.x)[1]
  eta <- surv[,1] - predict.x %*% fit$beta[-1]-fit$beta[1]

  order.eta <- order(eta)
  eta.order <- eta[order.eta]
  y.order <- surv[order.eta,1]
  state.order <- surv[order.eta,2]
  x.order <- predict.x[order.eta,,drop = FALSE]

  state.order[nobs] <- 1

  newx <- x.order[state.order==1,,drop=FALSE]
  newy <- y.order[state.order==1]
  newweight <- newy * 0 + 1

  eta.surv <- Surv(1:length(eta.order), state.order)
  km.eta <- survfit(eta.surv~1)
  weight.low <- km.eta$surv
  weight.up <- diff(-c(1, km.eta$surv))

  ################# data augmentation ###############
  eta.order.obs <- eta.order[state.order==1]
  repeats <- rev(cumsum(rev(state.order)))[state.order==0]
  newx.add <- apply(x.order[state.order==0,,drop = FALSE],2,function(t){rep(t,repeats)})
  y.order.list <- as.list(as.data.frame(t(cbind((1:nobs)[state.order==0], repeats))))
  addy <- function(t){
    res <- rep(x.order[t[1],,drop = FALSE] %*% fit$beta[-1]+fit$beta[1], t[2]) + eta.order[((t[1]+1):nobs)[state.order[(t[1]+1):nobs]==1]]
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

  ############# loss #############################

  if (method == "square"){
    loss <- sum((newy-newx %*% fit$beta[-1]-fit$beta[1])^2*newweight)
  } else if (method == "absolute"){
    loss <- sum(abs(newy-newx %*% fit$beta[-1]-fit$beta[1])*newweight)
  } else if (method == "huber"){
    residuals <- newy-newx %*% fit$beta[-1]-fit$beta[1]
    residuals <- residuals/sigma
    loss <- sum(loss.hb_all(residuals,k=k_h)*newweight)
  } else if (method == "tukey"){
    residuals <- newy-newx %*% fit$beta[-1]-fit$beta[1]
    residuals <- residuals/sigma
    loss <- sum(loss.tk_all(residuals,k=k_t)*newweight)
  }

  loss
}
