#### exponential
# pdf
dexponential <- function(params, x){
  lambda = params[1]
  lambda*exp(-lambda*x)
}

# cdf
pexponential <- function(params,x){
  lambda = params[1]
  1 - exp(-lambda*x)
}



#### weibull
# pdf
dweibull <- function(params, x){
  k = params[1]
  lambda = params[2]
  (k/lambda)*(x/lambda)^(k-1)*exp(-(x/lambda)^k)
}

# cdf
pweibull <- function(params,x){
  k = params[1]
  lambda = params[2]
  1 - exp(-(x/lambda)^k)
}


#### gamma
# pdf
dgamma1 <- function(params, x){
  alpha = params[1]
  beta = params[2]
  dgamma(x, shape=alpha, rate = beta)
}

# cdf
pgamma1 <- function(params,x){
  alpha = params[1]
  beta = params[2]
  pgamma(x, shape=alpha, rate = beta)
}



#### log-logistic
# pdf
dloglogis <- function(params, x){
  alpha = params[1]
  beta = params[2]
  ((beta/alpha)*(x/alpha)^(beta-1))/((1+(x/alpha)^beta)^2)
}

# cdf
ploglogis <- function(params,x){
  alpha = params[1]
  beta = params[2]
  (x^beta)/(x^beta+alpha^beta)
}


####  type2 - gumble weibull exponential
# pdf
dt2gwe <- function(params, x){
  alpha = params[1]
  beta = params[2]
  gamma = params[3]
  alpha*beta*gamma*(exp(gamma*x)-1)^(-beta-1)*exp(gamma*x-alpha*(exp(gamma*x)-1)^(-beta))
}

# cdf
pt2gwe <- function(params,x){
  alpha = params[1]
  beta = params[2]
  gamma = params[3]
  exp(-alpha*(exp(gamma*x)-1)^(-beta))
}


#### weibull generalized exponential
wge_pdf<-function(params, x){
  alpha = params[1]
  theta = params[2]
  gamma = params[3]
  alpha*theta*gamma*exp(gamma*x)*(exp(gamma*x)-1)^(theta-1)*exp(-alpha*(exp(gamma*x)-1)^theta)
}

wge_cdf<-function(params, x){
  alpha = params[1]
  theta = params[2]
  gamma = params[3]
  1 - exp(-alpha*(exp(gamma*x)-1)^theta)
}


#### Exponentiated Gumbel Type-2
eg2_pdf <- function(params,x){
  alpha = params[1]
  phi = params[2]
  theta = params[3]
  alpha*phi*theta*x^(-phi-1)*exp(-theta*x^(-phi))*(1-exp(-theta*x^(-phi)))^(alpha-1)
}

eg2_cdf <- function(params,x){
  alpha = params[1]
  phi = params[2]
  theta = params[3]
  1 - (1-exp(-theta*x^(-phi)))^alpha
}


#### Lomax Gumbel Type -Two
lg2_pdf <- function(params, x){
  alpha = params[1]
  beta = params[2]
  theta = params[3]
  k = params[4]
  alpha*beta^alpha*theta*k*x^(-k-1)*exp(-theta*x^(-k))*(beta-log(1-exp(-theta*x^(-k))))^(-alpha-1)/(1-exp(-theta*x^(-k)))
}

lg2_cdf <- function(params, x){
  alpha = params[1]
  beta = params[2]
  theta = params[3]
  k = params[4]
  1 - beta^alpha*(beta - log(1-exp(-theta*x^(-k))))^(-alpha)
}


## Gumbel type-two distribution
t2g_pdf <- function(params,x){
  alpha = params[1]
  nu = params[2]
  alpha*nu*x^(-nu-1)*exp(-alpha*x^(-nu))
}

t2g_cdf <- function(params,x){
  alpha = params[1]
  nu = params[2]
  exp(-alpha*x^(-nu))
}



#### EXPONENTIATED WEIBULL-LOGISTIC DISTRIBUTION
ewl_pdf <- function(params,x){
  alpha = params[1]
  beta = params[2]
  lambda = params[3]
  theta = params[4]
  theta*(1-exp(-alpha*exp(lambda*beta*x)))^(theta-1)*(lambda*alpha*beta*exp(lambda*beta*x-alpha*exp(lambda*beta*x)))
}

ewl_cdf <-function(params,x){
  alpha = params[1]
  beta = params[2]
  lambda = params[3]
  theta = params[4]
  (1-exp(-alpha*exp(lambda*beta*x)))^theta
}



#################
T2GWE_LL <- function(param){
  alpha = param[1]
  beta = param[2]
  gamma = param[3]
  LnL <- length(x)*log(alpha) + length(x)*log(beta) + sum(log(gamma*exp(-gamma*x))) - (beta+1)*sum(log(1-exp(-gamma*x))) + (beta-1)*sum(-gamma*x) - alpha*sum((exp(gamma*x)-1)**(-beta))
  LnL <- -LnL
  return(LnL)
}


WGE_LL <- function(param){
  alpha = param[1]
  theta = param[2]
  gamma = param[3]
  LL <- length(x)*log(alpha*theta*gamma) + gamma*sum(x) + (theta-1)*sum(log(exp(gamma*x)-1)) - alpha*sum((exp(gamma*x)-1)^theta)
  LL <- -LL
  LL
}


EG2_LL <- function(param){
  alpha = param[1]
  phi = param[2]
  theta = param[3]
  n = length(x)
  LnL <- n*log(alpha*phi*theta) - theta*sum(x^(-phi)) - (phi+1)*sum(log(x))  + sum(log((1-exp(-theta*x^(-phi)))^(alpha-1)))
  LnL <- -LnL
  return(LnL)    
}



LG2_LL <- function(param){
  alpha = param[1]
  beta = param[2]
  theta = param[3]
  k = param[4]
  n = length(x)
  LnL <- n*log(alpha*k*theta)+n*alpha*log(beta) - (k+1)*sum(log(x)) - theta*sum(x^(-k)) - sum(log(1-exp(-theta*x^(-k)))) - (alpha+1)*sum(log(beta-log(1-exp(-theta*x^(-k)))))
  LnL <- -LnL
  return(LnL)    
}



T2G_LL <- function(params){
  alpha = params[1]
  nu = params[2]
  LnL <- length(x)*log(nu) + length(x)*log(alpha) - (nu+1)*sum(log(x)) - alpha*sum(x^(-nu))
  LnL <- -LnL
  return(LnL)
}


############## Estimation function ############## 
############### MLE ############
MLE_fun_param <- function(param){
  alpha = param[1]
  beta = param[2]
  gamma = param[3]
  LnL <- n*log(alpha) + n*log(beta) + sum(log(gamma*exp(-gamma*x))) - (beta+1)*sum(log(1-exp(-gamma*x))) + (beta-1)*sum(-gamma*x) - alpha*sum((exp(gamma*x)-1)**(-beta))
  LnL <- -LnL
  return(LnL)
}


LS_fun_param <- function(param) {
  alpha = param[1]
  beta = param[2]
  gamma = param[3]
  LS <- sum((exp(-alpha*(exp(gamma*x)-1)**(-beta))-ls_ins)**2)
  LS
}


WLS_fun_param <- function(param) { 
  alpha = param[1]
  beta = param[2]
  gamma = param[3]
  WLS <- sum(wls_ins*(exp(-alpha*(exp(gamma*x)-1)**(-beta))-ls_ins)**2)
  WLS
} 


MPS_fun_param <- function(param) { 
  alpha = param[1]
  beta = param[2]
  gamma = param[3]
  x1 = x[1:length(x)-1]
  x2 = x[2:length(x)]
  MPS <- -alpha*(exp(gamma*x[1])-1)^{-beta} + log(1-exp(-alpha*(exp(gamma*x[n])-1)^{beta})) + sum(log(exp(-alpha*(exp(gamma*x2)-1)^{-beta}) - exp(-alpha*(exp(gamma*x1)-1)^{-beta})))
  MPS <- -MPS
  MPS
} 

CVM_fun_param <- function(param) { 
  alpha = param[1]
  beta = param[2]
  gamma = param[3]
  CVM <- 1/n*sum((exp(-alpha*(exp(gamma*x)-1)**(-beta))-ins)**2)
  CVM
} 



ADE_fun_param <- function(param) { 
  alpha = param[1]
  beta = param[2]
  gamma = param[3]
  ADE <- sum(i21*(-alpha*(exp(gamma*x)-1)**(-beta) + log(1-exp(-alpha*(exp(gamma*x_rev)-1)**(-beta))) ))
  ADE <- -ADE
  ADE
} 





