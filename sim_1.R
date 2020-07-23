library(autodiffr)
library(JuliaCall)
library(rHanso)

ad_setup(JULIA_HOME = '/Applications/Julia-1.3.app/Contents/Resources/julia/bin/')


n_ <- 2000
d_ <- 1

results_ <- matrix(NA,100,d_+1)

for (rr in 1:100) {
  X_ <- cbind(1,matrix(rnorm(n_*5),n_,d_))
  beta_ <- rep(1,d_+1)
  prob_fun <- function(x_,beta_) {
    exp(-max(1-sum(x_*beta_),0))/(
      exp(-max(1-sum(x_*beta_),0)) + 
        exp(-max(1+sum(x_*beta_),0))
    )
  }
  y_ <- c()
  for (ii in 1:n_) {
    rando_ <- runif(1)
    y_[ii] <- (rando_<prob_fun(X_[ii,],beta_)) -
      (rando_>=prob_fun(X_[ii,],beta_))
  }
  
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
                method='BFGS',control = list(maxit = 20000,reltol=1e-6))
  results_[rr,] <- opt_$par
  print(c(rr,results_[rr,]))
}
save(results_,file='sim_n2000d1.rdata')

