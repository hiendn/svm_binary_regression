library(colorspace)

plot(c(-4,2),c(0,1),type='n',xlab='x',ylab='f(1|x)',main='n=1000')
grid()
for (ii in 1:100) {
  dens_fun <- function(x) {
    exp(-pmax(1-results_[ii,1]-results_[ii,2]*x,0))/
      (
        exp(-pmax(1-results_[ii,1]-results_[ii,2]*x,0)) +
          exp(-pmax(1+results_[ii,1]+results_[ii,2]*x,0))
      )
  }
  curve(dens_fun,-4,2,col=rainbow_hcl(100)[ii],add=T)  
}
dens_fun <- function(x) {
  exp(-pmax(1-1-1*x,0))/
    (
      exp(-pmax(1-1-1*x,0)) +
        exp(-pmax(1+1+1*x,0))
    )
}
curve(dens_fun,-4,2,col='black',add=T,lwd=3,lty=2)  
