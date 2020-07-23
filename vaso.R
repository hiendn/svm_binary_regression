library(catdata)
library(sandwich)
library(colorspace)
library(autodiffr)
library(JuliaCall)
library(rHanso)

data("Wells")
Wells[,1] <- as.numeric(Wells[,1])-1
Wells[,5] <- as.numeric(Wells[,5])-1

glm_ <- glm(switch~.,data=as.data.frame(Wells),family = binomial)
summary(glm_)
coeftest(glm_,vcov=sandwich)

ad_setup(JULIA_HOME = '/Applications/Julia-1.3.app/Contents/Resources/julia/bin/')

n_ <- dim(Wells)[1]
d_ <- 4
X_ <- cbind(1,as.matrix(Wells[,2:5]))
beta_ <- rep(1,d_+1)
prob_fun <- function(x_,beta_) {
  exp(-max(1-sum(x_*beta_),0))/(
    exp(-max(1-sum(x_*beta_),0)) + 
      exp(-max(1+sum(x_*beta_),0))
  )
}
y_ <- (Wells[,1]-0.5)*2

like_fun <- function(beta_) {
  reg_ <- X_%*%beta_
  like_ <- sum(-sapply(1-y_*reg_,function(x){max(x,0)})) -
    sum(log(exp(-sapply(1-reg_,function(x){max(x,0)})) + 
              exp(-sapply(1+reg_,function(x){max(x,0)}))))
  -like_/length(y_)
}
ad_like_fun <- ad_variant(like_fun)
grad_fun <- makeGradFunc(ad_like_fun,mode='reverse',
                         use_tape = TRUE,compiled = TRUE,
                         debug = FALSE)
hess_fun <- makeHessianFunc(ad_like_fun)

opt_ <- optim(rep(0,d_+1),fn=like_fun,gr=grad_fun,
              method='BFGS',control = list(maxit = 1000,reltol=1e-10))
opt_
# hanso_ <- hanso(fn=ad_like_fun,gr=grad_fun,x0=c(0,0),
#                 maxitgs = 10,strongwolfe = 1,normtol = 1e-16,
#                 scale=1,wolfe1 = 1e-16,wolfe2 = 1e-6,
#                 maxit = 10000)
# hanso_
# 
# class_fun <- function(beta_) {
#   reg_ <- X_%*%beta_
#   return(sapply(1-reg_,function(x){max(x,0)}))
# }

library(numDeriv)
one_like_fun <- function(beta_,x1_,y1_) {
  reg_ <- x1_%*%beta_
  like_ <- -max(1-y1_*reg_,0) -
    log(exp(-max(1-reg_,0))+exp(-max(1+reg_,0)))
  like_
}
B_mat_ <- matrix(0,length(beta_),length(beta_))
for (ii in 1:length(y_)) {
  grad_ <- grad(one_like_fun,opt_$par,x1_=X_[ii,],y1_=y_[ii])
  B_mat_ <- grad_%*%t(grad_) + B_mat_
}

se_ <- sqrt(diag(solve(hess_fun(opt_$par))%*%B_mat_%*%solve(hess_fun(opt_$par))/(n_^2)))

Wells2_ <- Wells
names(Wells2_) <- c('y','arsen','dist','edu','assoc')
Wells2_[,5] <- jitter(Wells2_[,5],1)
plot(Wells2_[,-1],
     col=rainbow_hcl(2)[c(Wells2_[,1]+1)],
     pch=c(1,4)[c(Wells2_[,1]+1)],lwd=2,cex=2)
