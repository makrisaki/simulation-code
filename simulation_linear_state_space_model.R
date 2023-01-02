library(foreach)
library(yuima)

# dX_t=-a*X_tdt+b*dW_1(t)
# dY_t=c*X_tdt+sigma*W_2(t)
# Y_t is observable

#True values
a<-1.5
b<-0.3
c<-1
sigma<-0.02


trueparam <- list(a=a, b=b, c=c, sigma=sigma)

mod <- setModel(drift = c("-a*X","c*X"),diffusion = matrix(c("b","0","0","sigma"),2,2),solve.variable = c("X","Y"))


gamma <- function(a,b,c,sigma){
  ((sigma^2)*a/c^2)*(sqrt(1+b^2*c^2/(sigma^2*a^2))-1)
}

alpha <- function(a,b,c,sigma){
  a+c^2*gamma(a,b,c,sigma)/sigma^2
}

#Kalman-Bucy filter
m_vec <-function(a,b,c,sigma,obs.data,DeltaY,n,h,m_start){
  n <- length(obs.data)
  x <- seq(0,by=h,length=n)
  alpha <- alpha(a,b,c,sigma)
  gamma <- gamma(a,b,c,sigma)
  
  
  M <- numeric(n)
  M[1] <- m_start
  
  time_seq <- (0:n)*h
  
  #for(i in 2:n){
  #  M[i] <- exp(-alpha*i*h)*m_start+(c*gamma/sigma^2)*exp(-alpha*(i*h-time_seq[1:(i-1)]))%*%DeltaY[1:(i-1)]
  #}
  
  #obj <- cumsum(exp(alpha * time_seq[1:(n-1)]) * DeltaY)
  #M[-1] <- exp(-alpha * (2:n) * h) * (m_start + (c*gamma/sigma^2) * obj)
  M <- filter(c(m_start, exp(-2*alpha*h) * (c*gamma/sigma^2) * DeltaY), 
              filter = exp(-alpha*h), method = "r")
  
  return(M)
}

get_est <- function(n,h){
  samp <- setSampling(delta=h, n=n)
  sim <- simulate(mod, sampling = samp,true.parameter = trueparam)
  #sim <- simulate(mod, sampling = samp,true.parameter = trueparam,xinit=c(1,1))
  obs.data <- as.vector(sim@data@zoo.data$"Series 2")
  DeltaY <- diff(obs.data)
  
  #estimation of
  sigma_est <- sqrt(sum(DeltaY^2)/(n*h))
  
  #quasi-likelihood function
  qlf <- function(a,b,c,sigma,m_start=0,remove=0){
    m <- m_vec(a,b,c,sigma,obs.data,DeltaY,n,h,m_start)[(1+remove):n]
    DeltaY <- DeltaY[(1+remove):n]
    H <- -sum((DeltaY-h*c*m)^2)
    #H <- m%*%DeltaY-h*c*sum(m^2)/2
    if(is.finite(H)==TRUE){return(H)}else{return(0)}
  }
  
  #result of estimation 1
  opt1<- function(x){return(qlf(x[1],x[2],c,sigma_est))}
  
  #result of estimation 2
  opt2<- function(x){return(qlf(x[1],x[2],c,sigma_est,m_start=1))}
  
  #result of estimation 3
  opt3<- function(x){return(qlf(x[1],x[2],c,sigma_est,m_start=1,remove=100))}
  
  #result of estimation 4
  opt4<- function(x){return(qlf(x[1],x[2],c,sigma_est,m_start=1,remove=1000))}
  
  res1 <- optim(c(0.1,0.1),opt1,method = "L-BFGS-B",lower = c(0.01,0.01),upper=c(10,10),control = list(fnscale = -1))
  res2 <- optim(c(0.1,0.1),opt2,method = "L-BFGS-B",lower = c(0.01,0.01),upper=c(10,10),control = list(fnscale = -1))
  res3 <- optim(c(0.1,0.1),opt3,method = "L-BFGS-B",lower = c(0.01,0.01),upper=c(10,10),control = list(fnscale = -1))
  res4 <- optim(c(0.1,0.1),opt4,method = "L-BFGS-B",lower = c(0.01,0.01),upper=c(10,10),control = list(fnscale = -1))
  
  return(c(sigma_est,res1$par,res2$par,res3$par, res4$par))
  #return(res3$par)
}


#if(0){
#res1 <- replicate(10000,get_est(n=100000,h=0.001,m_start=0))
t1 <- system.time(res_all <- foreach(i = 1:10000, .combine = "rbind") %dopar% 
                    get_est(n=10^6,h=0.0001))
print(t1)
#}

#res_all is the result.

