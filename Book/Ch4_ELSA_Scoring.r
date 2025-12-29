# Scoring algorithm for ELSA data in Chapter 4

# Ardo , UCL 2016


##################################################
# Data are not available. Access to ELSA data can 
# be requested via the Economic and Social 
# Data Service (www.esds.ac.uk).
#
# Data format needed for this code:
# 
#   id  state  age  sex  educ 
#   1    2      1    1   0
#   1    3      3    1   0
#   1    4      5    1   0
#   1    4      7    1   0
#   2    2      1    0   0
#   2    2      3    0   0    etc. 
# 
#
# Typically age is a shifted version of age in years
#
# The models below are defined by
# 1. Two-column matrices <Qind..> to determine which 
#    transitions are modelled by the intercept, age,
#    and covariates
# 2. Vectors <..free> to identify the rows of the 
#    free parameters in <Qind..>
# 3. Matrices <..fixed> with no. rows equal to the 
#    length of vector <..free> to specify the 
#    link between free and fixed parameters 
#
#
####################################################

# Prelim:
digits <- 3
set.seed(12345)

# Prepare data for quick access:
da <- new.env()
for(i in 1:N){
  data.i <- dta[dta$id==subjects[i],]
  eval(parse(text=paste(paste("da$d",
        as.character(i),sep=""),"<-data.i",sep="")))
}
da <- as.list(da)
cat("Data prepared.\n")

# Data info:
subjects <- as.numeric(names(table(dta$id)))
N <- length(subjects)
cat("Sample size =",N,"\n")
cat("Frequencies observed state:"); print(table(dta$state))
cat("State table:"); print(statetable.msm(state,id,data=dta))

# Definitions for states and transitions:
# Dead state:
D <- 5
# Intercepts:
Qind   <- cbind(c(1,1,2,2,2,3,3,3,4,4),c(2,D,1,3,D,2,4,D,3,D))
npar   <- nrow(Qind)
qnames <- rep("q0",npar)
# Intercept only?:
NULLMODEL <- FALSE
if(NULLMODEL){
 # Total no. of parameters:
 nparTotal <- npar
 MODEL <- 0
}else{
 # Age:
 QindA <- cbind(c(1,1,2,2,3,3,4),c(2,D,3,D,4,D,D))
 # Choose baseline distributions (only for age-dependent transitions!):
 # Select by using/not-using #-sign: 
 MODEL <- 1; baseline <- rep("Gomp",7)
 #MODEL <- 2; baseline <- c("Gomp","Weib", "Gomp","Weib","Gomp","Weib","Weib")
 #MODEL <- 3; baseline <- rep("Weib",7)
 #MODEL <- 4; baseline <- c("Weib", "Gomp","Weib","Gomp","Weib","Gomp","Gomp")
 # Check:
 if(length(baseline)!=nrow(QindA)){
   stop("\nERROR. Number of specified baseline distr. not correct.\n")
 }
 A.free <- c(1,2,3,5)
 nparA  <- length(A.free)
 A.fixed      <- matrix(0,length(A.free),nrow(QindA))
 A.fixed[1,1] <- 1
 A.fixed[2,]  <- c(0,1,0,1,0,1,1)
 A.fixed[3,3] <- 1
 A.fixed[4,5] <- 1
 qnames <- c(qnames,rep("qAge",nparA))
 # Total no. of parameters:
 nparTotal <- npar+nparA
}
# Number of covariates:
ncov <- 2
# Enforce ncov=0 for NULLMODEL:
if(NULLMODEL){ncov <- 0}
# Define covariates effects if any:
if(ncov>0){
 # X1:
 QindX1  <- cbind(c(1,1,2,2,3,3,4),c(2,D,3,D,4,D,D))
 X1.free <- c(1,2,3,5)
 nparX1  <- length(X1.free)
 X1.fixed      <- matrix(0,length(X1.free),nrow(QindX1))
 X1.fixed[1,1] <- 1
 X1.fixed[2,]  <- c(0,1,0,1,0,1,1)
 X1.fixed[3,3] <- 1
 X1.fixed[4,5] <- 1
 qnames <- c(qnames,rep("qX1",nparX1))
 # X2:
 QindX2   <- cbind(c(1,2,3),c(2,3,4))
 nparX2   <- nrow(QindX2)
 X2.fixed <- diag(nrow(QindX2)) 
 qnames   <- c(qnames,rep("qX2",nparX2))
 # Total no. of parameters:
 nparTotal <- npar+nparA+nparX1+nparX2
}


#############################################################
# Functions for scoring:
# Building Q matrix:
Qmatrix <- function(beta,x,x1,x2){
 Q <- matrix(0,D,D)
 for(k in 1:npar){ Q[Qind[k,1],Qind[k,2]] <- exp(beta[k]) }
 if(!NULLMODEL){
  for(k in 1:nparA){
   for(l in 1:ncol(A.fixed)){ 
    if(baseline[l]=="Gomp"){
     if(A.fixed[k,l]==1){
      Q[QindA[l,1],QindA[l,2]] <- Q[QindA[l,1],QindA[l,2]]*exp(beta[npar+k]*x)
     }
    }
    if(baseline[l]=="Weib"){
     if(A.fixed[k,l]==1){
        Q[QindA[l,1],QindA[l,2]] <- Q[QindA[l,1],QindA[l,2]]*exp(beta[npar+k])*x^(exp(beta[npar+k])-1)
     }
    }
   }       
  }
 }
 # Covariates:
 if(ncov>0){
  for(k in 1:nparX1){
   for(l in 1:ncol(X1.fixed)){ 
    if(X1.fixed[k,l]==1){
      Q[QindX1[l,1],QindX1[l,2]] <- Q[QindX1[l,1],QindX1[l,2]]*exp(beta[npar+nparA+k]*x1)
    }
   }       
  }
  for(k in 1:nparX2){
   for(l in 1:ncol(X2.fixed)){ 
    if(X2.fixed[k,l]==1){
      Q[QindX2[l,1],QindX2[l,2]] <- Q[QindX2[l,1],QindX2[l,2]]*exp(beta[npar+nparA+nparX1+k]*x2)
    }
   }       
  }
  }
  for(k in 1:(D-1)){Q[k,k] <- -sum(Q[k,])}
  # Return:
  return(Q)
}

# P matrix using eigenvalue decomposition:
matrixexp <- function(Q,t){
  Eigen <- eigen(Q)
  A <- Eigen$vectors
  d <- Eigen$values
  D <- diag(exp(d*t))
  # Return:
  return(A%*%D%*%solve(A))
}
# Building V:
V <- function(G,d,t){
 V0 <- matrix(0,D,D)
 for(i in 1:D){
  for(j in 1:D){
    if(i==j){
      V0[i,j] <- G[i,i]*t*exp(d[i]*t)
    }else{
      V0[i,j] <- G[i,j]*(exp(d[i]*t)-exp(d[j]*t))/(d[i]-d[j])
    }
  }
  }
  # Return:
  return(V0)
}
# Derivative P matrix using eigenvalue decomposition:
DerivP <- function(Q,t,x,x1,x2){
  # Decomposition:
  Eigen <- eigen(Q, symmetric=FALSE)
  A <- Eigen$vectors
  d <- Eigen$values
  # Derivative matrices:
  Deriv <- array(NA,c(nparTotal,D,D))
  for(k in 1:npar){
     derivQ <- matrix(0,D,D)
     derivQ[Qind[k,1],Qind[k,2]] <-        Q[Qind[k,1],Qind[k,2]]
     derivQ[Qind[k,1],Qind[k,1]] <-  -derivQ[Qind[k,1],Qind[k,2]]
     G <- solve(A)%*%derivQ%*%A
     Deriv[k,1:D,1:D] <- A%*%V(G,d,t)%*%solve(A)
  }
  if(!NULLMODEL){
  # Derivative for A:
  for(k in 1:nparA){
   derivQ <- matrix(0,D,D)
   for(l in 1:ncol(A.fixed)){ 
    if(A.fixed[k,l]==1){ 
     if(baseline[l]=="Gomp"){
      derivQ[QindA[l,1],QindA[l,2]] <-   x*Q[QindA[l,1],QindA[l,2]]
      derivQ[QindA[l,1],QindA[l,1]] <-  -derivQ[QindA[l,1],QindA[l,2]]
     }
     if(baseline[l]=="Weib"){
      derivQ[QindA[l,1],QindA[l,2]] <-  (1+exp(beta[npar+k])*log(x))*Q[QindA[l,1],QindA[l,2]]
      derivQ[QindA[l,1],QindA[l,1]] <-  -derivQ[QindA[l,1],QindA[l,2]]
     }
    }
   }
   G <- solve(A)%*%derivQ%*%A
   Deriv[npar+k,1:D,1:D] <- A%*%V(G,d,t)%*%solve(A)
  }
  }
  # Covariates:
  if(ncov>0){
  # Derivative for X1:
  for(k in 1:nparX1){
   derivQ <- matrix(0,D,D)
   for(l in 1:ncol(X1.fixed)){ 
    if(X1.fixed[k,l]==1){
     derivQ[QindX1[l,1],QindX1[l,2]] <-  x1*Q[QindX1[l,1],QindX1[l,2]]
     derivQ[QindX1[l,1],QindX1[l,1]] <- -derivQ[QindX1[l,1],QindX1[l,2]]
    }
   }
   G <- solve(A)%*%derivQ%*%A
   Deriv[npar+nparA+k,1:D,1:D] <- A%*%V(G,d,t)%*%solve(A)
  }
  # Derivative for X2:
  for(k in 1:nparX2){
   derivQ <- matrix(0,D,D)
   for(l in 1:ncol(X2.fixed)){ 
    if(X2.fixed[k,l]==1){
     derivQ[QindX2[l,1],QindX2[l,2]] <-  x2*Q[QindX2[l,1],QindX2[l,2]]
     derivQ[QindX2[l,1],QindX2[l,1]] <- -derivQ[QindX2[l,1],QindX2[l,2]]
    }
   }
   G <- solve(A)%*%derivQ%*%A
   Deriv[npar+nparA+nparX1+k,1:D,1:D] <- A%*%V(G,d,t)%*%solve(A)
  }
  }       
  # Return:
  return(Deriv)
}

##################################
# Scoring algorithm:

# Starting values:
# For models without covariates:
if(NULLMODEL){ p0 <- rep(-3,npar) }
if(!NULLMODEL & ncov==0 & MODEL==1){ 
  p0 <- c(rep(-3,npar),rep(0,nparA))
}
# For models with two covariates:
if(ncov>0){
 # For Gompertz model:
 if(MODEL==1){ 
  p0 <- c(rep(-3,npar),rep(0,nparA),rep(0,nparX1),rep(0,nparX2))
 } 
 # For model with Weibull mortality:
 if(MODEL==2){
   p0 <- c(-3,-10,-3,-3,-10,-3,-3,-10,-3,-10, 0,.5,0,0, rep(0,nparX1),rep(0,nparX2))
 }
 # For Weibull model:
 if(MODEL==3){
   p0 <- c(rep(-10,npar), rep(.5,nparA), rep(0,nparX1),rep(0,nparX2))
 }
 # For model with Gompertz mortality:
 if(MODEL==4){ 
   p0 <- c(-10,-3,-10,-10,-3,-10,-10,-3,-10,-3, 0.5,0,.5,.5, rep(0,nparX1),rep(0,nparX2))
 }
}
beta <- p0
cat("\nStarting values",beta,"\n")

# Parametrs for no. of interations:
iter <- 1
max.iter <- 50
diff <- rep(Inf,length(beta))

# Precision:
abstol <- 1e-06

# Function to build score and info:
Scoring <- function(beta){
  S <- rep(0,nparTotal)
  M <- matrix(0,nparTotal,nparTotal)
  loglik <- 0
  for(i in 1:N){ 
   # Data and prelim for subject i:
   eval(parse(text=paste("data.i<-",paste("da$d",as.character(i),sep=""),sep="")))
   O   <- data.i$state
   age <- data.i$age 
   x1  <- data.i$sex[1]
   x2  <- data.i$educ[1]
   # Loop over individual follow-up:
   for(j in 2:length(O)){
    t <- age[j]-age[j-1]
    x <- age[j-1]
    # Q matrix:
    Q <- Qmatrix(beta,x=x,x1=x1,x2=x2)
    # P matrix and derivative:
    P <- matrixexp(Q,t=t)
    Deriv <- DerivP(Q,t,x=x,x1=x1,x2=x2)
    # Likelihood contribution:
    if(O[j]==D){
      contribution <- P[O[j-1],1:(D-1)]%*%Q[1:(D-1),D]
    }else{
      contribution <- P[O[j-1],O[j]]
    }
    loglik <- loglik+log(contribution)
    # Scoring contribution:
    S.i <- rep(0,nparTotal)
    if(O[j]==D){
       # For all parameters:
       denom <- P[O[j-1],1:(D-1)]%*%Q[1:(D-1),D]
       for(k in 1:nparTotal){
         S.i[k] <- Q[1:(D-1),D]%*%Deriv[k,O[j-1],1:(D-1)]/denom
       }
       # For intercept:
       for(k in 1:npar){
         if(Qind[k,2]==D){ S.i[k] <- S.i[k] + Q[Qind[k,1],D]*P[O[j-1],Qind[k,1]]/denom}
       }
       if(!NULLMODEL){
       # For age:
       for(k in 1:nparA){ 
        for(l in 1:ncol(A.fixed)){ 
           if(A.fixed[k,l]==1){
              if(QindA[l,2]==D){
                   if(baseline[l]=="Gomp"){
                      S.i[npar+k] <- S.i[npar+k] + x*Q[QindA[l,1],D]*P[O[j-1],QindA[l,1]]/denom
                   }
                   if(baseline[l]=="Weib"){
                      S.i[npar+k] <- S.i[npar+k] + (1+exp(beta[npar+k])*log(x))*Q[QindA[l,1],D]*P[O[j-1],QindA[l,1]]/denom
                   }
               }
            }
           }                      
         }
        }
       # Covariates:
       if(ncov>0){
         # For X1:
         for(k in 1:nparX1){ 
            for(l in 1:ncol(X1.fixed)){ 
              if(X1.fixed[k,l]==1){
                if(QindX1[l,2]==D){
                   S.i[npar+nparA+k] <- S.i[npar+nparA+k] + x1*Q[QindX1[l,1],D]*P[O[j-1],QindX1[l,1]]/denom
                }
              }
           }
         }
         # For X2:
         for(k in 1:nparX2){ 
           for(l in 1:ncol(X2.fixed)){ 
              if(X2.fixed[k,l]==1){
                if(QindX2[l,2]==D){
                   S.i[npar+nparA+nparX1+k] <- S.i[npar+nparA+nparX1+k] + x2*Q[QindX2[l,1],D]*P[O[j-1],QindX2[l,1]]/denom
                }
              }
           }
        }       
      }
    }
    if(O[j]!=D){
       # For all parameters:
       denom <- P[O[j-1],O[j]]
       for(k in 1:nparTotal){ S.i[k] <- Deriv[k,O[j-1],O[j]]/denom }
    }
    # Update scoring and information:
    S <- S+S.i
    for(k in 1:nparTotal){
      for(l in 1:nparTotal){
        M[k,l] <- M[k,l] + S.i[k]*S.i[l]
      }
    }
    }
  }
  return(list(S=S,M=M,minus2loglik=-2*loglik))
}

# Scoring algorithm:
while(iter<=max.iter & sum(diff) > abstol){
 # Scoring iteration
 SandM <- Scoring(beta)
 S <- SandM$S
 M <- SandM$M
 minus2loglik <- SandM$minus2loglik
 # Update beta:
 beta.old <- beta
 beta <- beta+solve(M)%*%S
 
 # Monitor convergence
 cat("Iteration =",iter,"\n-2Loglik =",minus2loglik,"\n")
 diff <- abs(beta.old-beta)
 cat("Updated parameters:\n",round(c(beta),digits),"\n")
 cat("Sum of absolute difference with old ones =",sum(diff),"\n\n")
 
 # Next iteration:
 iter <- iter+1
}

# Results:
cat("\nResults using",iter-1,"iterations (transformed parameters)\n")
cat("and starting values: ",p0,"\n")
p    <- beta
p.se <- sqrt(diag(solve(M)))
print(cbind(q=qnames,p=round(p,digits),
            se=round(p.se,digits)),quote=FALSE)

# Inference for Weibull tau-parameter in Model 2:
if(MODEL==2){
   variance <- diag(solve(M))
   variance[12] <- variance[12]*(exp(p[12])^2)
   cat("For tau:\n")
   print(c(exp(p[12]),sqrt(variance[12])))
}

# Information criterions:
cat("With -2*loglik = ",minus2loglik,"\n")
cat("and AIC = ",minus2loglik+2*length(beta),"\n")
cat("and BIC = ",minus2loglik+log(N)*length(beta),"\n")

