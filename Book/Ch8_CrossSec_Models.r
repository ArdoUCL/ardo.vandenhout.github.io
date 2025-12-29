# Fitting MSSMs to cross-sectional data from the Netherlands: models

# Ardo, UCL 2016


# Prelim:
digits <- 2

# Model:
Parameters <- function(p){
    beta <- c(p[1],p[2],p[1],p[2])
    xi   <- c(p[3],0,p[3],0)
    list(beta=beta,xi=xi)
}
pnames <- c("qF","qB","xiF")
p0     <- c(-2,-2,0)

# Q matrix:
Qmatrix <- function(beta,xi,t){
  q12 <- exp(beta[1]+xi[1]*t)
  q21 <- exp(beta[2]+xi[2]*t)
  q23 <- exp(beta[3]+xi[3]*t)
  q32 <- exp(beta[4]+xi[4]*t)
  matrix(c(-q12,q12,0, q21,-(q21+q23),q23, 0,q32,-q32),3,3,byrow=TRUE)
}

# Loglikelihood:      
loglikelihood <- function(p){
   # Extract parameters:
   param <- Parameters(p)
   beta  <- param$beta
   xi    <- param$xi
   # Contribution per year:
   loglik <- 0
   for(i in 2:n){
    # Prevalence:
    pie <- c(freq1[i-1],freq2[i-1],freq3[i-1])/size[i-1]
    x   <- c(freq1[i],freq2[i],freq3[i])
    # Qmatrix:
    Q <- Qmatrix(beta,xi,age[i-1])
    # One-step probs:
    P <- MatrixExp(Q,t=1)
    # Update:
    loglik <- loglik+log(dmultinom(x=x,prob=pie%*%P))
   }
 # Return:
 -loglik
}

# Minimise:
max <- optim(par=p0, fn=loglikelihood, method = "Nelder-Mead",
             control=list(maxit=5000,trace=FALSE),hessian=TRUE)
cat("\nStarting values =", p0,"\n")
cat("-2log(likelihood) =", round(2*max$value,digits),"\n")
conv <- max$convergence
cat("Convergence code =", conv,"\n")
p    <- max$par 
p.se <- sqrt(diag(solve(max$hessian)))
print(cbind(parameter=pnames,p=round(p,digits),se=round(p.se,digits)),
      quote=FALSE)

# Extract parameters:
param <- Parameters(p)
beta  <- param$beta
xi    <- param$xi
# Construct estimated one-year P matrices:
P <- array(NA,c(n,3,3))
for(i in 1:n){
    # Qmatrix:
    Q <- Qmatrix(beta,xi,age[i])
    # One-step probs:
    P[i,,] <- MatrixExp(Q,t=1)
 }

# Prediction:
pie <- c(freq1[1],freq2[1],freq3[1])/size[1]
pie.pred     <-matrix(NA,n,3)
pie.pred[1,] <-pie
P.pred       <-array(NA,c(n,3,3))
P.pred[1,,] <-diag(3)
for(i in 2:n){
    # Qmatrix:
    Q <- Qmatrix(beta,xi,age[i-1])
    # One-step probs:
    P.pred[i,,] <- P.pred[i-1,,]%*%MatrixExp(Q,t=1)
    # Prediction:
    pie.pred[i,] <- pie%*%P.pred[i,,]
}

# Plot prediction?:
cex.lab <- 1.2
plot(c(age[1]-shift,age[n]-shift),c(0,1),type="n",xlab="Age",
       ylab="Prevalence",cex.lab=cex.lab, main="")
pch <- c(15,17,19)
greys <- gray.colors(n=3, start = 0.2, end = 0.7, gamma = 2.2)
for(i in 1:3){
  lines(age-shift,pie.pred[,i],lwd=lwd,col=greys[i],lty=1)
  y <- dta[,i+1]/size
  points(age-shift,y,pch=pch[i],col=greys[i])
}
legend(34, 1, c("Normal weight","Overweight ","Obese"), col = greys,
       text.col = 1, lty = NULL, pch = pch, cex=1, bg = "white",
       text.width = strwidth("0000000000000"))
  

# Example P-matrix:
age.i <- c(25,35) + shift
for(i in 1:length(age.i)){
  cat("\nFor age", age.i[i] -shift,"one-year P-matrix is: \n")
  Q <- Qmatrix(beta,xi,age.i[i])
  print(round(MatrixExp(Q,t=1),digits))
}



