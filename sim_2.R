library(autodiffr)
library(JuliaCall)
library(rHanso)
library(MixSim)
library(e1071)

ad_setup(JULIA_HOME = '/Applications/Julia-1.3.app/Contents/Resources/julia/bin/')

results_ <- matrix(NA,100,3)

for (rr in 1:100) {
  n_ <- 1000
  d_ <- 2
  N_ <- 1000
  model_ <- MixSim(BarOmega = 0.5,K=2,p=d_)
  data_ <- simdataset(n = n_, 
                      Pi = model_$Pi, 
                      Mu = model_$Mu, 
                      S = model_$S, n.out = 0)
  data2_ <- simdataset(n = N_, 
                       Pi = model_$Pi, 
                       Mu = model_$Mu, 
                       S = model_$S, n.out = 0)
  X2_ <- cbind(1,data2_$X)
  y2_ <- 2*(data2_$id-1.5)
  y_ <- 2*(data_$id-1.5)
  X_ <- cbind(1,data_$X)
  beta_ <- rep(1,d_+1)
  
  like_fun <- function(beta_) {
    reg_ <- X_%*%beta_
    like_ <- sum(-sapply(1-y_*reg_,function(x){max(x,0)})) -
      sum(log(exp(-sapply(1-reg_,function(x){max(x,0)})) + 
                exp(-sapply(1+reg_,function(x){max(x,0)}))))
    -like_/length(y_)
  }
  ad_like_fun <- ad_variant(like_fun)
  grad_fun <- makeGradFunc(ad_like_fun)
  
  opt_ <- optim(rep(0,d_+1),fn=like_fun,gr=grad_fun,
                method='BFGS',control = list(maxit = 20000,reltol=1e-16))
  # opt_
  
  prob_fun <- function(x_,beta_) {
    exp(-max(1-sum(x_*beta_),0))/(
      exp(-max(1-sum(x_*beta_),0)) + 
        exp(-max(1+sum(x_*beta_),0))
    )
  }
  
  prob_ <- c()
  for (ii in 1:N_) {
    prob_[ii] <- prob_fun(X2_[ii,],opt_$par)
  }
  
  preds2_ <- 2*(round(prob_)-0.5)
  holder_ <- X_
  glm_ <- glm((y_+1)/2~holder_-1,family=binomial(link = 'logit'))
  holder_ <- X2_
  glm_preds2_ <- 2*(round(predict(glm_,as.data.frame(holder_),type='response'))-0.5)
  svm_ <- svm(X_[,-1],as.factor(y_),kernel = 'linear')
  svm_preds2_ <- 2*(as.numeric(predict(svm_,X2_[,-1]))-1.5)
  
  results_[rr,1] <- mean(preds2_==y2_)
  results_[rr,2] <- mean(glm_preds2_==y2_)
  results_[rr,3] <- mean(svm_preds2_==y2_)
  
  print(c(rr,results_[rr,]))
} 
# save(results_,file='class_n100d5w5.rdata')
