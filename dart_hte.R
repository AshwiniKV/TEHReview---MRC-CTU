# Simulation to introduce Treatment effect heterogeneity

N<-500
numx <- 5
alpha<-0.8
theta<-0.8
beta<- c(1,.8,.6,.4,.2,0,0,0)
gamma <- c(0, 1, 0, 0) #might have to change this depending on scenario
rho<-0 #might have to change this depending on scenario
eta<-0
version="prognostic"
Y <- list()
X <- list()
trt.eff<-list()
X.i <- list()
W <- list()
M <- list()
mu <- list()

set.seed(205)

for(i in 501:1000){
  Z <- rep(c(0, 1), each = N/2)
  
  ### Define the correlation matrix
  sigma <- diag(numx)
  sigma[upper.tri(sigma) | lower.tri(sigma)] <- rho ## Compound symmetry
  
  X.i[[i]] <- mvrnorm(N,mu=rep(0,numx),Sigma=sigma)
  #X.i[[i]] <- matrix(rbinom(N*numx,size=1,prob=0.5), nrow=N)
  
  #might have to change this depending on simulation scenario
  X[[i]] <- cbind(X.i[[i]],
                  X.i[[i]][,1] * X.i[[i]][,2], 
                  sqrt(abs(X.i[[i]][,3])),
                  as.integer(X.i[[i]][,4] > 0))
  
  #might have to change this depending on simulation scenario
  W[[i]] <- cbind(Z * ifelse(X[[i]][,1] > 0, 1, 0),
                  Z * X[[i]][,1],
                  Z * eta*(sin(X[[i]][,1])),
                  Z * (X[[i]][,1] > -0.5 & X[[i]][,1] < 0.5))
  
  #might have to change this depending on simulation scenario
  M[[i]] <- ifelse(X.i[[i]] > 0, 1, 0)
  
  mu[[i]] <- alpha + theta * Z + X[[i]] %*% beta + W[[i]] %*% gamma + M[[i]] %*% phi
  
  Y[[i]] <- rnorm(N, mean=mu[[i]])
  
  #might have to change this depending on simulation scenario
  trt.eff[[i]] <- theta + cbind(ifelse(X[[i]][,1] > 0, 1, 0),
                                X.i[[i]][,1],
                                ifelse(X[[i]][,1]>0 & X[[i]][,2]>0, 1, 0),
                                (X[[i]][,1] > -0.75 & X[[i]][,1] < 0.75)) %*% gamma 
} 

# Causal BART
library(BART)
xtr<-data.frame(x)
xte<-data.frame(Z, x)
xte<-xte[which(xte$Z == 1),]
xte$Z<-0
bart.tot<-wbart(x.train = xtr, y.train = y, x.test = xte)
catt_mat<-bart.tot$yhat.train[,which(Z ==1)] - bart.tot$yhat.test
catite<-apply(catt_mat, 2, mean)
mndiff<-apply(bart.tot$yhat.train[,which(Z ==1)] - bart.tot$yhat.test,1, mean)
credible_int <- apply(catt_mat, 1, function(row) HPDinterval(as.mcmc(row), prob=0.95))

# Causal DART
library(distr)

# DIrichlet prior 
dprior <- function(alpha, F0, n=1e3) { # n should be large since it's an approx for +oo
  
  s <- F0@r(n)                     # step 1: draw from F0
  V <- rbeta(n,1,alpha)            # step 2: draw from beta(1,alpha)
  w <- c(1, rep(NA,n-1))           # step 3: compute 'stick breaking process' into w[i]
  w[2:n] <- sapply(2:n, function(i) V[i] * prod(1 - V[1:(i-1)]))
  
  # return the sampled function F which can be itself sampled 
  # this F is a probability mass function where each s[i] has mass w[i]
  function (size=1e4) {
    sample(s, size, prob=w, replace=TRUE)
  }
}

# Dirichlet posterior
dp_posterior <- function(alpha, F0, data) {
  n <- length(data)
  F_n <- DiscreteDistribution(data) # compute empirical cdf
  
  F_bar <- n/(n+alpha) * F_n + alpha/(n+alpha) * F0
  
  dprior(alpha+n, F_bar)
}

f0 <- function(n) rnorm(n, 0, 3)    # the prior guess (in pdf format)
F0 <- DiscreteDistribution(f0(1e3)) # the prior guess (in cdf format)

y_train<-bart.tot$yhat.train[, which(Z ==1)] - bart.tot$yhat.test
y_bart<-y_train[1,]
distr<-apply(y_train, 2, mean)

data<-y_bart

runs  <- 500
xs    <- seq(-1.5,2.05,len=runs)

xs    <- seq(-1.25,1.75,len=runs)
y_hat <- matrix(nrow=length(xs), ncol=runs)
for(i in 1:runs) {
  Fpost <- dp_posterior(10, F0, data)  
  y_hat[,i] <- ecdf(Fpost())(xs)    # just keeping values of cdf(xs), not the object
  print(i)
}

