######the squared loss + SGL
linear_sgl = function( X, Y, weights, index, lam1, lam2, standardize=FALSE, maxiter=10^3, eps=1e-5) {

  n = dim(X)[1]
  p = dim(X)[2]

  K = length(unique(index))

  pk = NULL

  if (K == 1) {
    pk[1] = sum(index == unique(index)[1])
    K_loc = list(seq(1:length(index))[index == unique(index)[1]])
  }

  if (K > 1){
    pk[1] = sum(index == unique(index)[1])
    K_loc = list(seq(1:length(index))[index == unique(index)[1]])

    for (i in 2:K) {
      pk[i] = sum(index == unique(index)[i])
      K_loc[[i]] = seq(1:length(index))[index == unique(index)[i]]
    }
  }

  x = X
  y = Y

  w = weights/(sum(weights))

  if (standardize == TRUE){
    meanx = apply(x, 2, function(t){weighted.mean(t,w)})
    x = scale(x, meanx, FALSE)
    normx = apply(x, 2, function(t){sqrt(sum(w*t^2))})
    x = scale(x, FALSE, normx)
    meany = weighted.mean(y,w)
    y = y - mean(y)
  }

  if (standardize == FALSE){
    meanx = apply(x, 2, function(t){weighted.mean(t,w)})
    x = scale(x, meanx, FALSE)
    normx = 1
    meany = weighted.mean(y,w)
    y = y - mean(y)
  }

  w = weights

  #beta= rep(1,p)
  fit <- lm(y~x,weights=w)
  beta <- fit$coefficients[-1]
  iter = 0
  dif = 1

  while (iter <= maxiter & dif >= eps) {

    beta_old = beta

    for (i in 1:K) {

      rev = K_loc[[i]]

      r_group = y - x%*%beta + x[,rev]%*%as.matrix(beta[rev])
      sum_k = t(x[,rev])%*%as.matrix(r_group*w)
      con = St(sum_k,lam2)

      if (sqrt(sum(con^2)) <= sqrt(pk[i])*lam1) {beta[K_loc[[i]]] = 0}
      else {
        for (j in 1:pk[i]) {

          r_single = y - x%*%as.matrix(beta) + as.matrix(x[,K_loc[[i]][j]]*beta[K_loc[[i]][j]])

          save = sum(x[,K_loc[[i]][j]]*r_single*w)

          if (abs(save) <= lam2)  {beta[K_loc[[i]][j]]=0}
          else {
            beta[K_loc[[i]][j]]= St(save,lam2)/(sum(x[,K_loc[[i]][j]]^2*w)+sqrt(pk[i])*lam1/max(1e-6,sqrt(sum(beta[K_loc[[i]]]^2))))
          }
        }
      }

    }

    iter=iter+1
    #dif = sqrt(sum((beta_old-beta)^2))
    dif = max(abs(beta-beta_old))
  }

  coefficients = NULL
  beta = beta/normx
  #a0 = meany - sum(beta*meanx)
  a0 = sum((Y - X%*%beta)*weights)/sum(weights)
  varsel = abs(beta)>1e-5

  return(list(a0=a0,beta=beta,varsel=varsel,iter=iter))


}



######the absolute loss + SGL
abs_socp  <-  function(X, Y, w, index, lam1) {

  #X <-  cbind(x)
  nn   <-  length(Y)
  n  <-  nrow(X)
  stopifnot(n == nn)

  p  <-  ncol(X)

  K = length(unique(index))

  pk = NULL

  if (K == 1) {
    pk[1] = sum(index == unique(index)[1])
    K_loc = list(seq(1:length(index))[index == unique(index)[1]])
  }

  if (K > 1){
    pk[1] = sum(index == unique(index)[1])
    K_loc = list(seq(1:length(index))[index == unique(index)[1]])

    for (i in 2:K) {
      pk[i] = sum(index == unique(index)[i])
      K_loc[[i]] = seq(1:length(index))[index == unique(index)[i]]
    }
  }


  f  <-  c(w,rep(0,p),lam1*sqrt(pk))

  N <- c(rep(1,3*n),pk)

  dims <- list(l = 0, q = N+1 , e = 0)

  A <- matrix(0,nrow=3*n,ncol=n+p+K)

  C <- rbind(cbind(diag(n),X),cbind(diag(n),-X),cbind(diag(n),matrix(0,nrow=NROW(X),ncol=NCOL(X))))
  C <- cbind(C,matrix(0,nrow=NROW(C),ncol=K))

  G <- matrix(NA,nrow=6*n,ncol=n+p+K)
  G[seq(1,6*n-1,2),] <- (-C)
  G[seq(2,6*n,2),] <- (-A)

  newG <- NULL
  for (m in 1:K){
    z <- NULL

    for (s in 1:pk[m]){
      com <- rep(0,p-1)
      com[K_loc[[m]][s]] <- 1
      z <- rbind(z,com) }

    newA <- cbind(matrix(0,nrow=pk[m],ncol=n+1),z,matrix(0,nrow=pk[m],ncol=K))

    com1 <- rep(0,K)
    com1[m] <- 1
    newC <- c(rep(0,n+p),com1)

    newG <- rbind(newG,-newC,-newA)
  }

  G <- rbind(G,newG)
  G <- as(G, "dgCMatrix")

  b <- rep(0,3*n+p-1)
  d <- c(-Y,Y,rep(0,n+K))
  h <- c(d[1],b[1:sum(N[1])])
  for (l in 2:length(N)){
    h <- c(h,d[l],b[(sum(N[1:(l-1)])+1):sum(N[1:l])])
  }

  #fit <- ECOS_csolve(c=f,G=G,h=h,dims=dims)
  fit <- ECOSolveR::ECOS_csolve(c=f,G=G,h=h,dims=dims,control=ecos.control(maxit = 200L, feastol = 1e-08, reltol = 1e-08,abstol = 1e-08))
  sol <- fit$x
  coef <- sol[(n+1):(n+p)]
  list(coef=coef)
}



lad_socp_sgl <- function(x,y,weights=weights,index,lam1,lam2,eps=1e-6){

  nobs <- NROW(x)

  x <- cbind(matrix(1,nobs,1), x)

  p <- NCOL(x)

  xnew <- x
  for (i in 2:p){
    z <- rep(0,p)
    z[i] <- lam2
    xnew <- rbind(xnew,z)
  }

  ynew <- c(y,rep(0,p-1))

  weights <- c(weights,rep(1,p-1))


  fit <- abs_socp(xnew,ynew,weights,index,lam1)
  beta <- fit$coef

  beta <- beta*(abs(beta) > eps)

  a0 <- beta[1]
  beta <- beta[-1]


  return(list(a0=a0,beta=beta))

}



######the Huber loss and the Tukey loss + SGL
hbtk_cda_sgl <- function(X, Y, weights, index, lam1, lam2, standardize=TRUE, method=c("huber", "tukey"), maxit=10^3, eps=1e-5, k_h=1.345, k_t=4.685){

  n <- dim(X)[1]
  p <- dim(X)[2]

  K <- length(unique(index))

  pk <- NULL

  if (K == 1) {
    pk[1] <- sum(index == unique(index)[1])
    K_loc <- list(seq(1:length(index))[index == unique(index)[1]])
  }

  if (K > 1){
    pk[1] <- sum(index == unique(index)[1])
    K_loc <- list(seq(1:length(index))[index == unique(index)[1]])

    for (i in 2:K) {
      pk[i] <- sum(index == unique(index)[i])
      K_loc[[i]] <- seq(1:length(index))[index == unique(index)[i]]
    }
  }

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

  #beta <- rep(1,p)
  #a0 <- 0

  fit <- rqPen::rq.lasso.fit(x,y,lambda=0,tau=0.5,intercept=T,weights=w,method="br")
  beta <- fit$coefficients[-1]
  a0 <- fit$coefficients[1]

  w <- weights

  while (iter <= maxit & dif >= eps){

    beta_old <- beta
    a0_old <- a0

    for (i in 1:K) {

      r <- (y - a0 - x%*%as.matrix(beta))
      sigma <- weightedMad(r,w)
      r <- r/sigma

      rev <- K_loc[[i]]

      if (method == "huber"){
        sum_k <- t(x[,rev])%*%as.matrix(w*d.hb(r,k_h)/sigma + 2*w*x[,rev]%*%as.matrix(beta[rev])/sigma^2)
      } else if (method == "tukey") {
        sum_k <- t(x[,rev])%*%as.matrix(w*d.tk(r,k_t)/sigma + 2*w*x[,rev]%*%as.matrix(beta[rev])/sigma^2)
      }

      con <- St(sum_k,lam2/sigma)

      if (sqrt(sum(con^2)) <= sqrt(pk[i])*lam1/sigma) {beta[K_loc[[i]]] = 0}
      else {
        for (j in 1:pk[i]) {

          r <- (y - a0 - x%*%as.matrix(beta))
          sigma <- weightedMad(r,w)
          r <- r/sigma

          if (method == "huber"){
            save <- sum(w*d.hb(r,k_h)*x[,K_loc[[i]][j]])/sigma + sum(2*w*x[,K_loc[[i]][j]]^2)*beta[K_loc[[i]][j]]/sigma^2
          } else if (method == "tukey") {
            save <- sum(w*d.tk(r,k_t)*x[,K_loc[[i]][j]])/sigma + sum(2*w*x[,K_loc[[i]][j]]^2)*beta[K_loc[[i]][j]]/sigma^2
          }


          if (abs(save) <= lam2/sigma)  {beta[K_loc[[i]][j]]=0}
          else {
            beta[K_loc[[i]][j]]= St(save,lam2/sigma)/(sum(2*w*x[,K_loc[[i]][j]]^2/sigma^2)+sqrt(pk[i])*lam1/sigma/max(1e-6,sqrt(sum(beta[K_loc[[i]]]^2))))
          }
        }
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
