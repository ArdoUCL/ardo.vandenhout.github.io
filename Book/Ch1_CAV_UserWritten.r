# Chapter 1: CAV-history data analysis. 
# Userwritten code for parameter estimation

# Ardo, 2016

# Prelim:
library(msm)
digits <- 3

# Rename data, and give info:
dta <- cav
subjects <- as.numeric(names(table(dta$PTNUM)))
N <- length(subjects)
cat("Sample size =",N,"\n")
cat("Frequencies observed state:"); print(table(dta$state))
cat("State table:"); print(statetable.msm(state,PTNUM,data=dta))

# Add baseline age:
bage <- rep(NA,nrow(dta))
for(i in 1:N){
  select.i <- dta$PTNUM==subjects[i]
  bage[select.i] <- dta$age[select.i][1]
}
dta <- cbind(cav,bage=bage)

# Create history-state data:
Ostate <- dta$state
for(i in 1:N){
  dta.i <- dta[dta$PTNUM==subjects[i],]
  Hstate.i <- rep(NA, nrow(dta.i))
  for(j in 1:nrow(dta.i)){
    Hstate.i[j] <- max(dta.i$state[1:j])
  }
  dta$state[dta$PTNUM==subjects[i]] <- Hstate.i
}
dta <- cbind(dta,Ostate=Ostate)
# History data info:
cat("\nFrequencies observed history of state:"); print(table(dta$state))
cat("State table:"); print(statetable.msm(state,PTNUM,data=dta))

# Prepare data for quick access:
dta.split<-split(dta,dta$PTNUM)
cat("Data pre-formatted.\n")

# Prelim:
dead <- 4

# Likelihood:
loglikelihood<-function(p){
 # Intensities parameters:
 beta0  <- p[1:5]
 beta1  <- p[6:7]
 beta2  <- p[8:9]
 beta3  <- p[10:14]
 # Contribution per unit:
 loglik <- 0
 for(i in 1:N){
   # Data and prelim for subject i:
   data.i <- dta.split[[i]]
   O    <- data.i$state
   t    <- data.i$years
   bage <- data.i$bage[1]
   dage <- data.i$dage[1]
   # Loop over observations for subject i: 
   for(j in 2:length(O)){
     # Q and P matrix:
     Q <- matrix(0,4,4)
     Q[1,2] <- exp(beta0[1]+beta1[1]*t[j-1]+beta2[1]*bage+beta3[1]*dage)
     Q[1,4] <- exp(beta0[2]+beta1[2]*t[j-1]+beta2[2]*bage+beta3[2]*dage)
     Q[1,1]<- -sum(Q[1,])
     Q[2,3] <- exp(beta0[3]+beta3[3]*dage)
     Q[2,4] <- exp(beta0[4]+beta3[4]*dage)
     Q[2,2]<- -sum(Q[2,])
     Q[3,4] <- exp(beta0[5]+beta3[5]*dage)
     Q[3,3]<- -sum(Q[3,])
     P <- MatrixExp(mat=Q,t=t[j]-t[j-1])
     # Likelihood contribution:
     if(O[j]!=dead){contribution<-P[O[j-1],O[j]]}
     if(O[j]==dead){
        contribution <- P[O[j-1],1]*Q[1,dead]+
              P[O[j-1],2]*Q[2,dead]+P[O[j-1],3]*Q[3,dead]
      }
    # Update likelihood:
    loglik <- loglik+log(contribution)
  }
 }
 # Monitoring:
 cat("-2*Loglik = ", -2*loglik,"\n")
 -loglik
}

# Starting values for the maximisation:
p0 <- c(rep(-5,5),rep(0,9))
# If <model> is the object produced by msm(), then
# you can choose  
p0 <- model$opt$par 
# to compare results

# Maximise:
max <- optim(par=p0, fn=loglikelihood, method = c("Nelder-Mead"),
             control=list(maxit=5000),hessian=TRUE)

# Generate output:
cat("\nModel with covariates age and dage:")
cat("\n-2loglik =", 2*max$value,"\n")
conv<-max$convergence
cat("Convergence code =", conv,"\n")
p<-max$par; p.se<-sqrt(diag(solve(max$hessian)))
qnames<-"q"
print(cbind(q=qnames,p=round(p,digits),
        se=round(p.se,digits),"Wald ChiSq"=round((p/p.se)^2,digits),
       "Pr>ChiSq"=round(1-pchisq((p/p.se)^2,df=1),digits)),quote=FALSE)

# Transition probs for 1 year:
beta0 <- p[1:5]
beta1 <- p[6:7]
beta2 <- p[8:9]
beta3 <- p[10:14]
cat("\nTransition probs for 1 year\n")
bage0  <- median(data$bage[data$firstobs==1])
dage0  <- median(data$dage[data$firstobs==1])
years0 <- 0
Q <- matrix(0,4,4)
Q[1,2] <- exp(beta0[1]+beta1[1]*years0+beta2[1]*bage0+beta3[1]*dage0)
Q[1,4] <- exp(beta0[2]+beta1[2]*years0+beta2[2]*bage0+beta3[2]*dage0)
Q[1,1] <- -sum(Q[1,])
Q[2,3] <- exp(beta0[3]+beta3[3]*dage0)
Q[2,4] <- exp(beta0[4]+beta3[4]*dage0)
Q[2,2] <- -sum(Q[2,])
Q[3,4] <- exp(beta0[5]+beta3[5]*dage0)
Q[3,3] <- -sum(Q[3,])
P <- MatrixExp(mat=Q,t=1) 
print(round(P,digits))
