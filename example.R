library(smoothmest)
library(MASS)
library(survival)
library(matrixStats)

###sample size
nobs <- 200

###regression parameters
p <- 30
alpha <- 1
beta <- c(1,1.5,2,2.5,3,array(0,10),1,-1,1,-1,1,array(0,10))

###group structure
index <- ceiling(1:p/5)

###covariance structure
#pau <- 0.5
#cov_x <- matrix(NA,p,p)
#for (i in 1:p){
#  for (j in 1:p){
#    cov_x[i,j] <- pau^(abs(i-j))
#  }
#}
#diag(cov_x) <- 1
pau <- 0.5
cov_x <- pau*matrix(1,p,p)
diag(cov_x) <- 1

###error distribution and signal to-noise ratio
error <- 4
snr <- 5
sigma <- as.numeric(sqrt(t(as.matrix(beta))%*%cov_x%*%as.matrix(beta)/snr/23.4))

m <- 100
set.seed(10*m)


###generate covariates
x <- mvrnorm(nobs, rep(0,p), cov_x)

if (error == 1){
  epsilon <- rnorm(nobs, 0, 1)
} else if (error == 2){
  epsilon <- rt(nobs, 3)
} else if (error == 3){
  epsilon <- rdoublex(nobs)
} else if (error == 4){
  epsilon <- NULL
  for (j in 1:nobs){
    z <- rbinom(1, 1, 0.9)
    if (z==1) {epsilon[j] <- rnorm(1, 0, 1)}
    else {epsilon[j] <- rnorm(1, 0, 15) }
  }
}

epsilon <- sigma * epsilon

tlogt.expect <- x %*% beta + alpha
tlogt <- tlogt.expect + epsilon

tau <- 950
tc <- runif(nobs, min = 0, max = tau)
tlogc <- log(tc)
censor <- sum(tlogc < tlogt)/nobs

logt <- pmin(tlogt, tlogc)
indi <- (tlogt <= tlogc)
y <- Surv(logt, indi)


###square loss
#oracle
olsfit1 <- bje_refit(x[,abs(beta) > 1e-5], y, method="square")
olsfit1$beta
#lasso
fitlasso1 <- bje_ly(x, y, method="square", lambda=60, standardize=T)
fitlasso1$beta
#sgl
fitsgl1 <- bje_sgl(x, y, index, 30 , 20, method="square")
fitsgl1$beta



###huber loss
#oracle
olsfit3 <- bje_refit(x[,abs(beta) > 1e-5], y, method="huber")
olsfit3$beta
#lasso
fitlasso3 <- bje_ly(x, y, method="huber", lambda=25, standardize=T)
fitlasso3$beta
#sgl
fitsgl3 <- bje_sgl(x, y, index, 20 , 5, method="huber")
fitsgl3$beta



###absolute loss
#oracle
olsfit2 <- bje_refit(x[,abs(beta) > 1e-5], y, method="absolute")
olsfit2$beta
#lasso
fitlasso2 <- bje_ly(x, y, method="absolute", lambda=20, standardize=T)
fitlasso2$beta
#sgl
memory.limit(size=56000)
fitsgl2 <- bje_sgl(x, y, index, 20 , 5, method="absolute")
fitsgl2$beta



###tukey loss
#oracle
olsfit4 <- bje_refit(x[,abs(beta) > 1e-5], y, method="tukey")
olsfit4$beta
#lasso
fitlasso4 <- bje_ly(x, y, method="tukey", lambda=20, standardize=T)
fitlasso4$beta
#sgl
fitsgl4 <- bje_sgl(x, y, index, 20 , 5, method="tukey")
fitsgl4$beta
