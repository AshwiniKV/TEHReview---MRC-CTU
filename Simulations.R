#############################
# Simulation study 
#############################

N<-1000
numx <- 5
alpha<-0.8
theta<-0.8
beta<- c(1,.8,.6,.4,.2,0,0,0)
gamma <- c(0, 1, 0, 0) #might have to change this depending on scenario
phi <- c(0, 0,0, 0, 0) #might have to change this depending on scenario
rho<-0 #might have to change this depending on scenario
version="prognostic"
Y <- list()
X <- list()
trt.eff<-list()
X.i <- list()
W <- list()
M <- list()
mu <- list()

set.seed(205)

for(i in 1:1000){
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
                  Z * ifelse(abs(X[[i]][,1]) < 0.8, 1, 0),
                  Z * ifelse(X[[i]][,1]>0 & X[[i]][,2]>0, 1, 0),
                  Z * (X[[i]][,1] > -0.75 & X[[i]][,1] < 0.75))
  
  #might have to change this depending on simulation scenario
  M[[i]] <- ifelse(X.i[[i]] > 0, 1, 0)
  
  mu[[i]] <- alpha + theta * Z + X[[i]] %*% beta + W[[i]] %*% gamma + M[[i]] %*% phi
  
  Y[[i]] <- rnorm(N, mean=mu[[i]])
  
  #might have to change this depending on simulation scenario
  trt.eff[[i]] <- theta + cbind(X[[i]][,1] * X[[i]][,2], 
                                ifelse(X[[i]][,2] > 0, 1, 0),
                                (X[[i]][,1] > -0.5 & X[[i]][,1] < 0.5), 
                                ifelse(X[[i]][,2] > 0, 1, 0)) %*% gamma 
} 

##########################################################
# Simulations 
##########################################################

#######################################################################################
# Causal Forest 
# Library - grf
# Use the omnibus test for heterogeneity 

library(grf)
cfobject<-lapply(1:100, function(i) causal_forest(X.i[[i]], Y[[i]], Z))
pvalue<-lapply(1:100, function(i) test_calibration(cfobject[[i]])[2,4])
type1error<-ifelse(unlist(pvalue) <= 0.05, 1, 0)
mean(type1error)

cfobject<-lapply(1:100, function(i) causal_forest(X.i[[i]], Y[[i]], Z))
pvalue<-lapply(1:100, function(i) test_calibration(cfobject[[i]])[2,4])
power<-ifelse(unlist(pvalue) <= 0.05, 1, 0)
mean(power)

varimp<-lapply(1:100, function(i) variable_importance(cfobject[[i]]))
varmax<-lapply(1:100, function(i) which.max(varimp[[i]]) == 1)
prob<-sum(unlist(varmax))/100
prob

#######################################################################################
# Support Vector Machine


library(FindIt)

subgroup<-function(out, cov, trt){
F1  <-FindIt(model.treat= out~trt,
             model.main= ~ cov,
             model.int= ~ cov,
             type="continuous",
             treat.type="single",
             search.lambdas=FALSE,
             lambdas = c(-3.8,-4.0)) 
pred1<-predict(F1)
tab1<-table(ifelse(pred1$data$Treatment.effect> 1.05, 1, 0))/1000*100
sbg<-ifelse(tab1[[1]]>15, 1, 0)
sbg
}
rate<-lapply(1:100, function(i) subgroup(Y[[i]], X.i[[i]], Z))
mean(unlist(rate))
#######################################################################################
# TEHTrees 
# Single Sample 
# Double Sample
# Does better at larger sample sizes due to the loss of observations in matching 

Nused <- 200
subjects <- c(1:(Nused/2), 1001:(1000+Nused/2))
res <- lapply(1:1000, function(i) LTfunction(Y[[i]][subjects], Z[subjects], X[[i]][subjects,], X.i[[i]][subjects,]))
subgrp<-lapply(1:1000, function(i) !is.na(res[[i]]$splitleft))
mean(unlist(subgrp))

ttree<-function(out, cov, trt){
matched <- MatchForTree(Y=out,
                        Z = trt,
                        X=cov,
                        version="prognostic")
ymatch <- matched$Y.match
xmatch <- matched$X.match
itrt <- matched$itrt
ictl <- matched$ictl

LT<-LMEtree(ymatch, 
        as.matrix(xmatch), 
        ictl, 
        1:nrow(as.matrix(xmatch)), pval.thresh = 0.05)
LT
}

ttree_res<-lapply(1:100, function(i) ttree(Y[[i]], X.i[[i]], Z))
subgrp<-lapply(1:100, function(i) !is.na(ttree_res[[i]]$splitleft))
mean(unlist(subgrp))

######################################################################################

# Causal Tree

library(causalTree)
subgroup<-function(out, var, trt){
newdata<-data.frame(cbind(out, var))
names(newdata)<-c("Y", "X1", "X2", "X3", "X4", "X5")
tree <- causalTree(out~., data = newdata, treatment = trt,
                   split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5, 
                   cp = 0, minsize = 20, propensity = 0.5)
opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
opfit <- prune(tree, opcp)
#sub<-ifelse(opfit$frame$var[1] != "<leaf>", 1, 0)
#sub<-ifelse(opfit$frame$var[1] == "X1"| opfit$frame$var[1] == "X2", 1, 0)
sub<-ifelse(opfit$frame$var[1] == "X1", 1, 0)
sub
}
subgrp<-lapply(1:100, function(i) subgroup(Y[[i]], X.i[[i]], Z))
mean(unlist(subgrp))
#######################################################################################