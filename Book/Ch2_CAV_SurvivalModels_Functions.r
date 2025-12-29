# Functions to estimate univariate hazard models:

# Ardo, UCL 2016.


UniSurvModel <- function(idN,t,delta,X=NULL,model="exp",
               start.values=NULL,digits=3,method="Nelder-Mead",b=NULL){

# Checks:
N <- length(idN)
if(length(t)!=N){
   stop("\nERROR. Length <t> should be length <idN>.\n\n")
}
if(length(delta)!=N){
  stop("\nERROR. Length <delta> should be length <idN>.\n\n")
  }
if(!model%in%c("exp","weib","gomp")){
   stop("\nERROR. Choose model from <exp>, <weib>, <gomp>.\n\n")
} 
if(all(X[,1]==1)){
  stop("\nERROR. Don't include intercept column.\n\n")
} 
 
# Hazard, density, and Survivor functions:
# Exponential:
h.exp  <- function(t,lambda){rep(lambda, length(t))}
d.exp  <- function(t,lambda){lambda*exp(-lambda*t)}
S.exp  <- function(t,lambda){exp(-lambda*t)}
# Weibull:
h.weib <- function(t,lambda,tau){lambda*tau*t^(tau-1)}
d.weib <- function(t,lambda,tau){lambda*tau*t^(tau-1)*exp(-lambda*t^tau)}
S.weib <- function(t,lambda,tau){exp(-lambda*t^tau)}
# Gompertz:
eps    <- 0.000001
h.gomp <- function(t,lambda,xi){
  ifelse(abs(xi)<eps,h.exp(t,lambda),lambda*exp(xi*t))
}
d.gomp <- function(t,lambda,xi){ 
  lambda*exp(xi*t)*exp(-lambda*xi^(-1)*(exp(xi*t)-1))
}
S.gomp <- function(t,lambda,xi){ 
  ifelse(abs(xi)<eps,S.exp(t,lambda),exp(-lambda*(xi^-1)*(exp(xi*t)-1)))
}

# Function to display MLE:
display <- function(max,names){
 cat("\n-2loglik =", 2*max$val,". AIC = ",2*max$val+2*length(max$par),"\n")
 cat("Optimisation method =", method,"\n")
 cat("Convergence code =", max$con,"\n")
 cat("Parameters and SEs:\n")
 p    <- max$par
 p.se <- sqrt(diag(solve(max$hessian)))
 out  <-  round(cbind(p,p.se,"  Wald ChiSq"=(p/p.se)^2,
                      "Pr>ChiSq"=1-pchisq((p/p.se)^2,df=1)),digits)
 if(nrow(out)==(1+ncov)){out[1,3:4]   <- "-"}
 if(nrow(out)==(1+ncov+1)){out[1,3:4] <- "-"; out[ncov+2,3:4]<-"-"}
 print(cbind(names,out),quote=FALSE)
}

# Delta method for exp(p):
delta.exp <- function(max,index,name){
  p <- max$par 
  variance.t <- diag(solve(max$hessian))[index]; 
  variance   <-  variance.t*exp(p[index])^2
  cat(name," = ",round(exp(p[index]),digits)," with SE",
      round(sqrt(variance),digits),"\n")
}

# Number of covariates:
ncov    <- ncol(X)
x.names <- c("x1","x2","x3","x4","x5","x6")[1:ncov] 

# Add intercept and define <data>:
X    <- cbind(1,X)
data <- as.data.frame(cbind(idN,t,delta,X))

##################
# Exponential model:
if(model=="exp"){
loglikelihood<-function(p){
  beta   <- p
  loglik <- rep(NA,N)
  for(i in 1:N){
   lambda.i  <- exp(beta%*%X[i,1:(ncov+1)])
   loglik[i] <- log(S.exp(t[i],lambda.i)*h.exp(t[i],lambda.i)^delta[i])
  }
  -sum(loglik)
}
if(is.null(start.values)){p <- c(-5,rep(0,ncov))}else{p<-start.values}
max   <- optim(par=p, fn=loglikelihood, method = method,
               control=list(maxit=20000),hessian=TRUE)
names <- c("(Intercept)", x.names)
display(max,names)
}

################
# Weibull model:
if(model=="weib"){
loglikelihood <- function(p){
  beta   <- p[1:(ncov+1)]
  tau    <- exp(p[ncov+2])
  loglik <- rep(NA,N)
  for(i in 1:N){
   lambda.i  <- exp(beta%*%X[i,1:(ncov+1)])
   loglik[i] <- log(S.weib(t[i],lambda.i,tau)*h.weib(t[i],lambda.i,tau)^delta[i])
  }
  -sum(loglik)
}
if(is.null(start.values)){p <- c(-5,rep(0,ncov),.1)}else{p<-start.values}
max   <- optim(par=p, fn=loglikelihood, method = method,
               control=list(maxit=20000),hessian=TRUE)
names <- c("(Intercept)", x.names,"log(tau)")
display(max,names)
delta.exp(max,ncov+2,"tau")
}

##################
# Gompertz model:
if(model=="gomp"){
loglikelihood<-function(p){
  beta   <- p[1:(ncov+1)]
  xi     <- p[ncov+2]
  loglik <- rep(NA,N)
  for(i in 1:N){
   lambda.i  <- exp(beta%*%X[i,1:(ncov+1)])
   loglik[i] <- log(S.gomp(t[i],lambda.i,xi)*h.gomp(t[i],lambda.i,xi)^delta[i])
   }
  -sum(loglik)
}
if(is.null(start.values)){p <- c(-5,rep(0,ncov),.1)}else{p <- start.values}
max<-optim(par=p, fn=loglikelihood, method = method,
           control=list(maxit=20000),hessian=TRUE)
names <- c("(Intercept)", x.names,"xi")
display(max,names)
}

# Optional: Return estimated parameters:
# list(par=max$par)
}
