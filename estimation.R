rm(list=ls())
set.seed(1)
source(file = "~/functions.R")
library(maxLik)
library(stats4)
library(numDeriv)


T2GWE <- function(u, alpha, beta, gamma) {
  tmp = -alpha/(log(u))
  return(log(tmp^(1/beta)+1)/gamma)
}


alpha1 =2.5; beta1 =0.8; gamma1 =1.3
init_cond = c(2.1, 1.1, 1.5)

ns <- c(50, 100, 250, 500, 1000)
alphas <- c(); betas  <- c(); gammas <- c()
parms1 <- c(); parms2 <- c(); parms3 <- c()
mse1 <- c(); mse2 <- c(); mse3 <- c()

NN = 1000
for (j in 1:length(ns)){
  n <- ns[j]
  columns= c('param1','param2','param3')
  error = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(error) = columns
  
  for (k in 1:NN){
    Fx <- runif(n)
    x <- c(); ls_ins <- c(); wls_ins <- c(); ins <-c(); i21 <- c(); i2n1 <- c()
    for (i in 1:n) {
      ### simulation
      r <- T2GWE(Fx[i], alpha1, beta1, gamma1)
      x <- c(x, r)
      ins <- c(ins, (2*i-1)/(2*n))
      ls_ins <- c(ls_ins, i/(n+1))
      wls_ins <- c(wls_ins, (n+1)**2*(n+1)/(i*(n-i+1)))
      i21 <- c(i21, 2*i-1)
      i2n1 <- c(i2n1, 2*n+1-2*i)
    }
    x = sort(x)
    x_rev = sort(x, decreasing = TRUE)
    est_nlminb <- nlminb(init_cond, ADE_fun_param, lower=0.001)
    error[k,] = est_nlminb$par
  }
  
  parms1 = c(parms1, mean(error$param1)-alpha1)
  parms2 = c(parms2, mean(error$param2)-beta1)
  parms3 = c(parms3, mean(error$param3)-gamma1)

  mse1 = c(mse1, sum((error$param1-alpha1)**2)/NN)
  mse2 = c(mse2, sum((error$param2-beta1)**2)/NN)
  mse3 = c(mse3, sum((error$param3-gamma1)**2)/NN)
}


parms1
parms2
parms3

mse1
mse2
mse3



