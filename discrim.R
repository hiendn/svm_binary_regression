library(DAAG)
library(autodiffr)
library(JuliaCall)
library(rHanso)

set.seed(1000)
partition_ <- sample(c(rep(1,920),
                       rep(2,920),
                       rep(3,920),
                       rep(4,920),
                       rep(5,921)))

part_ <- 5

data("spam7")
data_ <- spam7
data_[,7] <- as.numeric(data_[,7])-1

train_ <- as.matrix(data_[-which(partition_==part_),])
test_ <- as.matrix(data_[which(partition_==part_),])

svm_ <- svm(spam7[-which(partition_==part_),-7],spam7[-which(partition_==part_),7],kernel = 'linear')

ad_setup(JULIA_HOME = '/Applications/Julia-1.3.app/Contents/Resources/julia/bin/')

n_ <- dim(train_)[1]
d_ <- 6
X_ <-cbind(1,train_[,-7])
beta_ <- rep(1,d_+1)
prob_fun <- function(x_,beta_) {
  exp(-max(1-sum(x_*beta_),0))/(
    exp(-max(1-sum(x_*beta_),0)) + 
      exp(-max(1+sum(x_*beta_),0))
  )
}
y_ <- (train_[,7]-0.5)*2

X2_ <-  cbind(1,test_[,-7])
y2_ <- (test_[,7]-0.5)*2

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
prob_fun <- function(x_,beta_) {
  exp(-max(1-sum(x_*beta_),0))/(
    exp(-max(1-sum(x_*beta_),0)) + 
      exp(-max(1+sum(x_*beta_),0))
  )
}

prob_ <- c()
for (ii in 1:dim(test_)[1]) {
  prob_[ii] <- prob_fun(X2_[ii,],opt_$par)
}
preds_ <- 2*(round(prob_)-0.5)
mean(preds_==y2_)

mean(predict(svm_,spam7[which(partition_==part_),-7])==spam7[which(partition_==part_),7])

