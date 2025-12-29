# MLE for Exponential, Weibull and Gompertz univariate hazard models:

# Ardo, UCL 2016.


# Prelim:
library(survival)
library(msm)
digits <- 3

# Load the functions:
source('Ch2_CAV_SurvivalModels_Functions.r')


# Data:
idN   <- as.numeric(names(table(cav$PTNUM)))
N     <- length(idN)
t     <- rep(NA,N)
delta <- rep(NA,N)
X <- matrix(NA,N,3)
for(i in 1:N){
  data.i   <- cav[cav$PTNUM==idN[i],]
  t[i]     <- data.i$years[nrow(data.i)]
  delta[i] <- as.numeric(data.i$state[nrow(data.i)]==4)
  X[i,] <- c(data.i$sex[1],data.i$dage[1],data.i$age[1])
}

# Exponential model:
cat("\nExponential model:")
UniSurvModel(idN,t,delta,X=X,model="exp",digits=digits)

# Weibull model:
cat("\nWeibull model:")
UniSurvModel(idN,t,delta,X=X,model="weib",digits=digits)

# Gompertz model:               
cat("\nGompertz model:")
UniSurvModel(idN,t,delta,X=X,model="gomp",digits=digits)

############################
# using survival()
# Function to print models:
print.m <- function(m){
   cat("-2Loglik = ",-2*m$loglik[2],"\n")
   cat("Hazard parameters:\n")
   print(round(c(-m$coeff/m$scale,shape=1/m$scale),digits))

}
# Prepare data:
X <- cbind(1,X)
data <- as.data.frame(cbind(idN,t,delta,X))
names(data) <- c("idN","t","delta","intercept","x1","x2","x3")
# Using survreg():
cat("\n\nUsing R package survival:\n")
cat("Exponential model:\n")
m.exp <- survreg(Surv(t,delta)~x1+x2+x3,data=data,dist="exponential")
print.m(m.exp)
cat("\nWeibull model:\n")
m.weib <- survreg(Surv(t,delta)~x1+x2+x3,data=data,dist="weibull")
print.m(m.weib)

# Fit Cox model:
cox.model <- coxph(Surv(t,delta)~x1+x2+x3,data=data)
cat("\n\nCox model. ")
print(cox.model)

