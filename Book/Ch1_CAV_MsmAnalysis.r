# Chapter 1: CAV-history data analysis
# Using msm() for parameters estimation

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

# Choose model:
Model <- 1
# Model formulation:
if(Model==1){
  # Generator matrix Q:
  q <- 0.01
  Q <- rbind(c(0,q,0,q), c(0,0,q,q),c(0,0,0,q),c(0,0,0,0))
  qnames <- c("q12","q14","q23","q24","q34")
  # Covariates:
  covariates <- as.formula("~years+bage+dage")
  constraint <- NULL
  fixedpars  <- c( 8:10,13:15)
  # Control:
  method <- "BFGS"
}

# Fit model using msm:
model <- msm(state~years, subject=PTNUM, data=dta, center=FALSE, 
            qmatrix=Q, death=TRUE, covariates=covariates, constraint=constraint,
            fixedpars=fixedpars, method=method, control=list(trace=0,REPORT=1))

# Generate output:
cat("\nModel",Model," with covariates: "); print(covariates)
cat("and constraints:\n"); print(constraint)
cat("and fixedpars:\n"); print(fixedpars)
cat("\n-2loglik =", model$minus2loglik,"\n")
conv <- model$opt$convergence; cat("Convergence code =", conv,"\n")
p    <- model$estimates; p.se <- sqrt(diag(model$covmat))
print(cbind(q=qnames,p=round(p,digits),
            se=round(p.se,digits),"Wald ChiSq"=round((p/p.se)^2,digits),
            "Pr>ChiSq"=round(1-pchisq((p/p.se)^2,df=1),digits)),quote=FALSE)

# Transition probs for 1 year:
cat("\nTransition probs for 1 year\n")
bage0 <- median(dta$bage[dta$firstobs==1])
dage0 <- median(dta$dage[dta$firstobs==1])
pmat  <- pmatrix.msm(model, t=1, t1=0, covariates=list(years=0,
                    bage=bage0,dage=dage0), ci="normal")
print(pmat)
