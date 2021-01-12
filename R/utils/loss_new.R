###huber
loss.hb_all <- function(r,k){
  l <- array(0,c(length(r),1))
  l[abs(r) <= k] = r[abs(r) <= k]^2/2
  l[abs(r) > k] = (2*abs(r[abs(r) > k])-k)*k/2
  l
}

w.hb_all <- function(r,k){
  l <- array(0,c(length(r),1))
  l[abs(r) <= k] = r[abs(r) <= k]*0+1
  l[abs(r) > k] = k/abs(r[abs(r) > k])
  l
}


###tukey
loss.tk_all <- function(r,k){
  l <- array(0,c(length(r),1))
  l[abs(r) <= k] = (1-(1-(r[abs(r) <= k]/k)^2)^3)*k^2/6
  l[abs(r) > k] = (0*abs(r[abs(r) > k])+k)^2/6
  l
}

w.tk_all <- function(r,k){
  l <- array(0,c(length(r),1))
  l[abs(r) <= k] = (1-(r[abs(r) <= k]/k)^2)^2
  l[abs(r) > k] = 0*abs(r[abs(r) > k])
  l
}

###soft-thresholding
St <- function(z,lam){
  p <- pmax(abs(z)-lam,0)*sign(z)
  return(p)
}

###the first derivative of the huber loss
d.hb <- function(r,k){
  l <- array(0,c(length(r),1))
  l[abs(r) <= k] = r[abs(r) <= k]
  l[abs(r) > k] = k*sign(r[abs(r) > k])
  l
}

###the first derivative of the tukey loss
d.tk <- function(r,k){
  l <- array(0,c(length(r),1))
  l[abs(r) <= k] = (1-(r[abs(r) <= k]/k)^2)^2*r[abs(r) <= k]
  l[abs(r) > k] = 0*abs(r[abs(r) > k])
  l
}

