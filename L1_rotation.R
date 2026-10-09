proj<- function(R){
  
  res= svd(R)
  res$u%*% t(res$v)
}

IRLS_update <- function(U, R, step=1){
  
  n=nrow(U)
  temp = U%*% R
  W <- 1/(abs(temp) + mean(abs(temp)) * 1e-3)
  
  obj0 = sum(W* temp^2)
  
  temp = W * (U%*% R)
  grad =  2*t(U) %*% temp/n
  
  R.new = proj(R - step*grad)
  temp = U%*% R.new
  obj1 = sum(W* temp^2)
  
  dif = obj0 - obj1
  z = 0
  while(d < 0 | z <10){
    step = step/2
    R.new = proj(R - step*grad)
    temp = U%*% R.new
    obj1 = sum(W* temp^2)
    
    dif = obj0 - obj1
    z = z+1
  }
  R.new
}

IRLS <- function(U, R, tol = 1e-4){
  n = nrow(U)
  obj0 = sum(abs(U%*% R))
  
  R = IRLS_update(U, R)
  obj1 = sum(abs(U%*% R))
  while(abs(obj0-obj1) > tol){
    #print(obj1)
    obj0 = obj1
    R = IRLS_update(U, R)
    obj1 = sum(abs(U%*% R))
  }
  list(R=R, Z= sqrt(n)*U%*% R)
}

mse.fun <- function(Z1, Z2){
  k= ncol(Z2)
  cov.matr = cov(Z1, Z)
  order <- rep(0,k)
  sign <- rep(1, k)
  for(i in 1:k){
    order[i] <- which.max(abs(cov.matr[i,]))
    if(cov.matr[i,order[i]]<0) sign[i] =-1
  }
  Z2.new = Z2[,order] %*% diag(sign)
  mean((Z1-Z2.new)^2)
}
library(GPArotation)


n = 100
d = n
k=5

mse.vsp.100 <- rep(0, 100)
mse.l1.100 <- rep(0, 100)

for(i in 1:100){
  print(i)
  set.seed(i)
  Z = matrix(rnorm(n*k, 0, sd = sqrt(2))*rbinom(n*k, 1, 0.5),n,k)
  Y = matrix(rnorm(d*k),d,k)

  A = Z%*% t(Y) + rnorm(n*d)

  res = svd(A)
  U = res$u[,1:k]
  R = diag(1, k)

  res = IRLS(U, R, tol=1e-5)
  res.vsp = varimax(U, normalize = F, eps = 1e-5)
  res.vsp.Z= sqrt(n) * res.vsp$loadings[1:n,]
    
  mse.l1.100[i] = mse.fun(res$Z, Z)
  mse.vsp.100[i] = mse.fun(res.vsp.Z, Z)
}


n = 200
d = n
k=5

mse.vsp.200 <- rep(0, 100)
mse.l1.200 <- rep(0, 100)

for(i in 1:100){
  print(i)
  set.seed(i)
  Z = matrix(rnorm(n*k, 0, sd = sqrt(2))*rbinom(n*k, 1, 0.5),n,k)
  Y = matrix(rnorm(d*k),d,k)
  
  A = Z%*% t(Y) + rnorm(n*d)
  
  res = svd(A)
  U = res$u[,1:k]
  R = diag(1, k)
  
  res = IRLS(U, R, tol=1e-5)
  res.vsp = varimax(U, normalize = F, eps = 1e-5)
  res.vsp.Z= sqrt(n) * res.vsp$loadings[1:n,]
  
  mse.l1.200[i] = mse.fun(res$Z, Z)
  mse.vsp.200[i] = mse.fun(res.vsp.Z, Z)
}


n = 100
d = n
k=5

mse.vsp.100.t <- rep(0, 100)
mse.l1.100.t <- rep(0, 100)

for(i in 1:100){
  print(i)
  set.seed(i)
  Z = matrix(rt(n*k, df = 5),n,k)/sqrt(5/3)
  Y = matrix(rnorm(d*k),d,k)
  
  A = Z%*% t(Y) + rnorm(n*d)
  
  res = svd(A)
  U = res$u[,1:k]
  R = diag(1, k)
  
  res = IRLS(U, R, tol=1e-5)
  res.vsp = varimax(U, normalize = F, eps = 1e-5)
  res.vsp.Z= sqrt(n) * res.vsp$loadings[1:n,]
  
  mse.l1.100.t[i] = mse.fun(res$Z, Z)
  mse.vsp.100.t[i] = mse.fun(res.vsp.Z, Z)
}


n = 200
d = n
k=5

mse.vsp.200.t <- rep(0, 100)
mse.l1.200.t <- rep(0, 100)

for(i in 1:100){
  print(i)
  set.seed(i)
  Z = matrix(rt(n*k, df = 5),n,k)/sqrt(5/3)
  Y = matrix(rnorm(d*k),d,k)
  
  A = Z%*% t(Y) + rnorm(n*d)
  
  res = svd(A)
  U = res$u[,1:k]
  R = diag(1, k)
  
  res = IRLS(U, R, tol=1e-5)
  res.vsp = varimax(U, normalize = F, eps = 1e-5)
  res.vsp.Z= sqrt(n) * res.vsp$loadings[1:n,]
  
  mse.l1.200.t[i] = mse.fun(res$Z, Z)
  mse.vsp.200.t[i] = mse.fun(res.vsp.Z, Z)
}


boxplot(c(mse.vsp.100, mse.l1.100, mse.vsp.200, mse.l1.200,mse.vsp.100.t, mse.l1.100.t, mse.vsp.200.t, mse.l1.200.t)~ c(rep(1,100), rep(2,100), rep(3,100), rep(4,100),rep(5,100), rep(6,100), rep(7,100), rep(8,100)))
save(mse.vsp.100, mse.l1.100, mse.vsp.200, mse.l1.200,mse.vsp.100.t, mse.l1.100.t, mse.vsp.200.t, mse.l1.200.t, file = "rotation.Rdata")
