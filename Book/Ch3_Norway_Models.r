# Three-state Gompertz and Weibull models with covariates

# Ardo, UCL 2016

# Using <data> in the long format:
#  id state age sex
#  1     1   0    0
#  1     1   1    1
#  1     1   2    2
#  1     2   3    3
#  2     1   0    0
#  2     1   1    1  etc.

# Load <Pmatrix.Gomp> and <Pmatrix.Weib> first by running script
# Ch3_Norway_Models_Functions.r

# Prelim:
library(msm)
set.seed(12345)
digits <- 2

# Choose model (2 for Gomp, 3 for Weib, 5  for restricted Weib):
M <- 2
qnames <- c("q12","q13","q23")

# Parameters for P matrix computation:
h.int  <- 1/12
eps.xi <- 0.001

# Data info:
subjects <- as.numeric(names(table(data$id)))
N <- length(subjects)
cat("Sample size =",N,"\n")
cat("Frequencies observed state:"); print(table(data$state))
cat("State table:"); print(statetable.msm(state,id,data=data))

# Prepare data:
da <- new.env()
for(i in 1:N){
  data.i <- data[data$id==subjects[i],]
  eval(parse(text=paste(paste("da$d",as.character(i),sep=""),"<-data.i",sep="")))
}
da <- as.list(da)
cat("Data prepared.\n")

# Prelim model:
Q <- rbind(c(0,1,1), c(0,0,1),c(0,0,0))
dead      <- 3
censored  <- -2

####################################
# Gompertz model:                  #
####################################
if(M==2){

# Loglikelihood:
loglikelihood <-function(p){
 # Parameters:
 beta0  <- p[1:3]
 xi     <- p[4:6]
 beta   <- p[7:12]
 # Vector for regression part:
 lambda <- rep(NA,3)
 # Contribution per subject:
 loglik <- 0
 for(i in 1:N){
  # Data and prelim for subject i:
  eval(parse(text=paste("data.i<-",paste("da$d",as.character(i),sep=""),sep="")))
  dur <- data.i$duration[1]
  sex <- data.i$sex[1]
  O   <- data.i$state
  age <- data.i$age
  # Loop over individual follow-up:
  for(j in 2:length(O)){
    # Regression:
    lambda[1] <- exp(beta0[1]+beta[1]*sex+beta[4]*dur)
    lambda[2] <- exp(beta0[2]+beta[2]*sex+beta[5]*dur)
    lambda[3] <- exp(beta0[3]+beta[3]*sex+beta[6]*dur)
    # Time interval:
    t1 <- age[j-1]
    t2 <- age[j]
    # Even number of nodes:
    nnodes <- max(4,round((t2-t1)/h.int))
    nnodes <- ifelse(nnodes/2==round(nnodes/2),nnodes,nnodes+1)
    h <- (t2-t1)/nnodes
    P <- Pmatrix.Gomp(t1=t1,t2=t2,lambda=lambda,xi,nnodes,h)
    # Likelihood contribution:
    if(O[j]==dead){
      q13 <- lambda[2]*exp(xi[2]*(t2-h))
      q23 <- lambda[3]*exp(xi[3]*(t2-h))
      contribution <- P[O[j-1],1]*q13+P[O[j-1],2]*q23
    }
    if(O[j]==censored){ contribution <- P[O[j-1],1]+P[O[j-1],2]}
    if(O[j]==1|O[j]==2){contribution <- P[O[j-1],O[j]]}
    # Update likelihood:
    loglik <- loglik+log(contribution)
  }
 }
 # Monitoring:
 cat("-2*Loglik = ", -2*loglik,"\n")
 -loglik
}

# Maximise:
p <- c(-5,-5,-5,0,0,0,rep(0,6))
max <- optim(par=p, fn=loglikelihood, method = c("BFGS"),
             control=list(maxit=3000, fnscale=2000),hessian=TRUE)
cat("\nModel",M,"\n")
cat("-2loglik =", 2*max$value,"\n")
conv <- max$convergence
cat("Convergence code =", conv,"\n")
p    <- max$par
p.se <- sqrt(diag(solve(max$hessian)))
print(cbind(q=qnames,p=round(p,digits),
        se=round(p.se,digits),"Wald ChiSq"=round((p/p.se)^2,digits),
      "Pr>ChiSq"=round(1-pchisq((p/p.se)^2,df=1),digits)),quote=FALSE)

}

####################################
# Weibull model:                  #
####################################
if(M%in%c(3,5)){
# Loglikelihood:
loglikelihood <- function(p){
 # Parameters:
 beta0 <- p[1:3]  
 if(M==3){
   tau <- exp(p[4:6]); beta <- p[7:12]
 }
 if(M==5){ 
   tau  <- c(exp(p[4]),1,exp(p[5])); 
   beta <- c(p[6],0,p[7],p[8],0,p[9])
  }
 # Vector for regression part:
 lambda <- rep(NA,3)
 # Contribution per subject:
 loglik <- 0
 for(i in 1:N){
  # Data and prelim for subject i:
  eval(parse(text=paste("data.i<-",paste("da$d",as.character(i),sep=""),sep="")))
  dur <- data.i$duration[1]
  sex <- data.i$sex[1]
  O   <- data.i$state
  age <- data.i$age
  # Loop over individual follow-up:
  for(j in 2:length(O)){
    # Regression:
    lambda[1] <- exp(beta0[1]+beta[1]*sex+beta[4]*dur)
    lambda[2] <- exp(beta0[2]+beta[2]*sex+beta[5]*dur)
    lambda[3] <- exp(beta0[3]+beta[3]*sex+beta[6]*dur)
    # Time interval:
    t1 <- age[j-1]
    t2 <- age[j]
    # Even number of nodes:
    nnodes <- max(4,round((t2-t1)/h.int))
    nnodes <- ifelse(nnodes/2==round(nnodes/2),nnodes,nnodes+1)
    h <- (t2-t1)/nnodes
    P <- Pmatrix.Weib(t1=t1,t2=t2,lambda=lambda,tau,nnodes,h)
     # Likelihood contribution:
    if(O[j]==dead){
       q13 <- lambda[2]*tau[2]*(t2-h)^(tau[2]-1)
       q23 <- lambda[3]*tau[3]*(t2-h)^(tau[3]-1)
       contribution <- P[O[j-1],1]*q13+P[O[j-1],2]*q23
    }
    if(O[j]==censored){contribution <- P[O[j-1],1]+P[O[j-1],2]}
    if(O[j]==1|O[j]==2){contribution <- P[O[j-1],O[j]]}
    # Update likelihood:
    loglik <- loglik+log(contribution)
  }
 }
 # Monitoring:
 cat("-2*Loglik = ", -2*loglik,"\n")
 -loglik
}

# Maximise:
if(M==3){p <- c(-15,-10,-10,log(4),log(2),log(3),rep(0,6))}
if(M==5){p <- c(-16.4, -10.1, -10.1, 1.5, 1.1,-0.5, -0.4, 0.1, 0.0)}
max <- optim(par=p, fn=loglikelihood, method = c("BFGS"),
             control=list(maxit=3000, fnscale=2000),hessian=TRUE)
cat("Model",M,"\n")
cat("-2loglik =", 2*max$value,"\n")
conv <- max$convergence
cat("Convergence code =", conv,"\n")
p <- max$par
variance <- diag(solve(max$hessian))
# Delta method:
if(M==5){
 for(i in 4:5){variance[i] <- variance[i]*(exp(p[i])^2)}
 p[4:5] <- exp(p[4:5])
}else{
 for(i in 4:6){variance[i] <- variance[i]*(exp(p[i])^2)}
 p[4:6] <- exp(p[4:6])
}
p.se <- sqrt(variance)
if(M==5){qnames <- NA}
print(cbind(q=qnames,p=round(p,digits),
       se=round(p.se,digits),"Wald ChiSq"=round((p/p.se)^2,digits),
      "Pr>ChiSq"=round(1-pchisq((p/p.se)^2,df=1),digits)),quote=FALSE)
}

