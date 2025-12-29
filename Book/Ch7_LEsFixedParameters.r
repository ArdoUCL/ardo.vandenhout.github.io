# Chapter 7. LEs for a predefined process

# Ardo, UCL 2016

# Prelim:
library(msm)
digits <- 3
logit  <- function(lp){exp(lp)/(1+exp(lp))}

# Number of states:
D <- 3

# Set paremeters and covariance matrix:
# Multi-state model:
beta0 <- c(-4.2,-4,-3)
xi    <- c(0.02,0.005,0.015)
cat("\nIntercept parameters:",beta0,"\n")
cat("Age-effect parameters:",xi,"\n")
# Logistic regression model:
gamma <- c(6,-0.07)
cat("\nLogistic model parameters:",gamma,"\n")


# Set up plot framework:
opar <- par(mfrow=c(3,1), mex=0.8,mar=c(5,5,1,1))

######################
# Plot from State 1:
# Plot graphs twice with different value for h:
for(run in 1:2){

# Grid parameter:
if(run==1){h <- 1}else{h <- 1/12}

# Age grid:
age0     <- 65
max.age  <- 115
age.grid <- seq(from=age0,to=max.age,by=h)
L <- length(age.grid)

# Plot step function survival:
P  <- array(NA,c(L,D,D))
PM <- array(NA,c(L,D,D))

# Loop over age grid:
# P-matrix per time interval:
P[1,1:D,1:D] <- diag(D)
for(j in 2:L){
    t1 <- age.grid[j-1]
    t2 <- age.grid[j]
    # Q and P matrix:
    Q <- matrix(0,D,D)
    Q[1,2] <- exp(beta0[1]+xi[1]*t1)
    Q[1,3] <- exp(beta0[2]+xi[2]*t1); Q[1,1] <- -(Q[1,2]+Q[1,3])
    Q[2,3] <- exp(beta0[3]+xi[3]*t1); Q[2,2] <- -Q[2,3]
    P[j,1:D,1:D] <- MatrixExp(mat=Q,t=t2-t1)
}
# P-matrices for whole time grid. And their rows:
PM[1,1:D,1:D] <- P[1,1:D,1:D]
for(l in 2:L){
    PM[l,1:D,1:D] <- PM[l-1,1:D,1:D]%*%P[l,1:D,1:D]
}

# Plot param for run:
if(run==1){
  lwd <- 5; col <- "grey"; lty <- c(1,1,1)
}else{
  lwd <- 3; col <- 1; lty <- c(1,2,3)
}
cex.lab <- 1.2
# Plot framework only once:
if(run==1){
  plot(c(age0,max.age),c(0,1),type="n", 
     ylab="Transition probability",xlab="",cex.lab=cex.lab)
}
# Plot transition probs:
lines(age.grid,PM[,1,1],lwd=lwd,col=col,lty=lty[1],type="s")
lines(age.grid,PM[,1,2],lwd=lwd,col=col,lty=lty[2],type="s")
lines(age.grid,PM[,1,3],lwd=lwd,col=col,lty=lty[3],type="s")

# Add legend only once:
if(run==2){
  legend(x=115,y=0.8, legend = c("to state 1", "to state 2", "to state 3"),
               text.width = strwidth("XXXXXXXXXXXX"),
               title="From baseline state 1",
               lty = lty, col=c(col,col,col), lwd=2,cex=1.2,
               xjust = 1, yjust = 1)
}
}

#######################
# From State 2:
# Plot graphs twice with different value for h:
for(run in 1:2){

# Grid parameter:
if(run==1){h <- 1}else{h <- 1/12}

# Age grid:
age0     <- 65
max.age  <- 115
age.grid <- seq(from=age0,to=max.age,by=h)
L <- length(age.grid)
  

# Plot step function survival:
P  <- array(NA,c(L,D,D))
PM <- array(NA,c(L,D,D))

# Loop over age grid:
# P-matrix per time interval:
P[1,1:D,1:D] <- diag(D)
for(j in 2:L){
    t1 <- age.grid[j-1]
    t2 <- age.grid[j]
    # Q and P matrix:
    Q<-matrix(0,D,D)
    Q[1,2] <- exp(beta0[1]+xi[1]*t1)
    Q[1,3] <- exp(beta0[2]+xi[2]*t1); Q[1,1] <- -(Q[1,2]+Q[1,3])
    Q[2,3] <- exp(beta0[3]+xi[3]*t1); Q[2,2] <- -Q[2,3]
    P[j,1:D,1:D] <- MatrixExp(mat=Q,t=t2-t1)
}
# P-matrices for whole time grid. And their rows:
PM[1,1:D,1:D] <- P[1,1:D,1:D]
for(l in 2:L){
    PM[l,1:D,1:D] <- PM[l-1,1:D,1:D]%*%P[l,1:D,1:D]
}

# Plot param for run:
if(run==1){
  lwd <- 5; col<-"grey"; lty <- c(1,1)
}else{
  lwd <- 3; col <- 1; lty <- c(1,2)
}
cex.lab <- 1.2
# Plot framework only once:
if(run==1){
  plot(c(age0,max.age),c(0,1),type="n", 
     ylab="Transition probability",xlab="",cex.lab=cex.lab)
}
# Plot transition probs:
lines(age.grid,PM[,2,2],lwd=lwd,col=col,lty=lty[1],type="s")
lines(age.grid,PM[,2,3],lwd=lwd,col=col,lty=lty[2],type="s")

# Add legend only once:
if(run==2){
  legend(x=115,y=0.8, legend = c("to state 2", "to state 3"),
               text.width = strwidth("XXXXXXXXXXXX"),
               title="From baseline state 2",
               lty = lty, col=c(col,col), lwd=2,cex=1.2,
               xjust = 1, yjust = 1)
}
}
               
###################################
# Plot baseline state probability
plot(c(age0,max.age),c(0,1),type="n", xlab="Age",
     ylab="Probability to be in state 1",cex.lab=cex.lab)
p1 <- logit(gamma[1]+gamma[2]*age.grid)
lines(age.grid,p1,lwd=lwd,col=col,lty=1,type="s")     
 
    
################################
# Using logistic regression:
p1 <- logit(gamma[1]+gamma[2]*age0)
cat("Probability to be in state 1 at age ",age0," is",round(p1,digits),"\n")

##################################
# LEs

# Estimation LEs step function:
e11 <- h*sum(PM[,1,1])
e12 <- h*sum(PM[,1,2])
e1  <- p1*e11
e22 <- h*sum(PM[,2,2])
e2  <- p1*e12+(1-p1)*e22
LEs <- c(e11,e12,e22,e1,e2)
labels <- c("e11","e12","e22","e1","e2")
cat("\nLEs using grid approximation with h =",round(h,digits),":\n")
print(cbind(labels,round(LEs,digits)),quote=FALSE)
e <- e1+e2
cat("Total LEs:",round(e,digits),"\n")


# Estimation LEs step Simpson's rule:
# Work with even number of intervals (so with uneven number of nnodes):
nnodes <- ifelse(round(L/2)!=L/2,L,L-1)

# Composite Simpson's rule:
# Defined for indexing by 0,1,2,..,n, with n even.
# Here adapted using <adapt> for indexing by 1,2,3,...,n+1:
simpson <- function(integrand,nnodes,h){
   n <- nnodes-1
   adapt <- 1
   S <- integrand[0+adapt] 
   for(j in seq(1,n/2-1,by=1)){  S <- S + 2 * integrand[2*j+adapt]} 
   for(j in seq(1,n/2,by=1)){    S <- S + 4 * integrand[2*j-1+adapt]} 
   S <- S+integrand[n+adapt]
   h*S/3
}

e11 <- simpson(PM[,1,1], nnodes,h)
e12 <- simpson(PM[,1,2], nnodes,h)
e1  <- p1*e11
e22 <- simpson(PM[,2,2], nnodes,h)
e2  <- p1*e12+(1-p1)*e22
LEs <- c(e11,e12,e22,e1,e2)
labels <- c("e11","e12","e22","e1","e2")
cat("\nLEs using Simpson's rule with h =",round(h,digits),":\n")
print(cbind(labels,round(LEs,digits)),quote=FALSE)
e <- e1+e2
cat("Total LEs:",round(e,digits),"\n")
