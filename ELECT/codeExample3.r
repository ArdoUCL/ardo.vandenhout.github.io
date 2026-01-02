# A multi-state survival model: msm() and ELECT. Example 3

# Ardo, UCL 2012

# Prelim:
library(msm)
digits <- 3

# Adapt the two following commands to your computer settings:
source("C:\\mydirectory\\ELECT.r")

# Data info:
dta <- cav
subjects <- as.numeric(names(table(dta$PTNUM)))
N <- length(subjects)
cat("Sample size =",N,"\n")
cat("Frequencies observed state:"); print(table(dta$state))
cat("State table:"); print(statetable.msm(state,PTNUM,data=dta))

# Choose model:
Model<-1

# Model definitions:
# Generator matrix Q:
q <-0.01
Q <- rbind(c(0,q,0,q), c(0,0,q,q),c(0,0,0,q),c(0,0,0,0))
qnames <- c("q12","q14","q23","q24","q34")
# Covariates:
covariates <- as.formula("~age+dage+sex")
constraint <- NULL 
fixedpars  <- c(6,8,9,10)
# MC matrix E:
ematrix <- rbind( c( 0, 0.1, 0 ,0), c( 0.1, 0, 0.1,0 ), 
                    c( 0, 0.1, 0,0 ), c( 0, 0, 0, 0) ) 
# Control:
method <- "BFGS"

# Fit the model:
model <- msm(state~age, subject=PTNUM, data=dta, center=FALSE, 
            qmatrix=Q, death=TRUE, covariates=covariates, constraint=constraint,
            ematrix=ematrix,fixedpars=fixedpars,  method=method, 
            control=list(trace=0,REPORT=1,maxit=1000,fnscale=100000))

# Generate output:
cat("\n-2loglik =", model$minus2loglik,"\n")
cat("Convergence code =", model$opt$convergence,"\n")
nbeta <- max(which(names(model$estimates)=="qcov"))
p    <- model$estimates[1:nbeta]
p.se <- sqrt(diag(model$covmat)[1:nbeta])
print(cbind(q=qnames,p=round(p,digits),
            se=round(p.se,digits)),quote=FALSE)
cat("\nMisclassification matrix (given mean of covariate values):\n")
print(round(ematrix.msm(model,ci="none"),digits))


# LEs:
age0   <- round(mean(dta[dta$firstobs==1,]$age))
dage0  <- round(mean(dta[dta$firstobs==1,]$dage))
sex0   <- median(dta[dta$firstobs==1,]$sex)
sddata <- dta[dta$firstobs==1,]
LEs <- elect(model=model, b.covariates=list(age=age0,dage=dage0,sex=sex0),
                 statedistdata=sddata, time.scale="years",
                 h=0.5, age.max=100, S=500,setseed=12345)
summary.elect(LEs, digits=2)



