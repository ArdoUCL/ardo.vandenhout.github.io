# A multi-state survival model: msm() and ELECT. Example 1

# Ardo, UCL 2012

# Preliminaries:
library(msm)
# Adapt the two following commands to your computer settings:
#load("C:\\mydirectory\\dataExample1.RData")
#source("C:\\mydirectory\\ELECT.r")

# Data info:
cat("Sample size:"); print(length(table(data$id)))
cat("Frequencies observed state:"); print(table(data$state))
cat("State table:"); print(statetable.msm(state,id,data=data))
q<-0.001; Q<-rbind(c(0,q,q), c(q,0,q),c(0,0,0))

# Model fit:
model<-msm(state~age, subject=id, data=data, center=FALSE, 
             qmatrix=Q, death=TRUE, covariates=~age+ybirth,
             censor= -2, censor.states=c(1,2),
             method="BFGS", control=list(trace=1,REPORT=1,maxit=1000,fnscale=100000))
qnames<-c("q12","q13","q21","q23")
p<-model$estimates; p.se<-sqrt(diag(model$covmat))
print(cbind(q=qnames,p=round(p,3),se=round(p.se,3)),quote=FALSE)

# Estimate life expectancies:
sddata<-data[data$state%in%c(1,2),]
LEs.pnt<-elect(model=model, b.covariates=list(age=0,ybirth=20),
                 statedistdata=sddata, time.scale.msm="years",
                 h=0.5, age.max=40, S=0)
summary.elect(LEs.pnt, digits=2)

# Plot LEs for an age range:
age.range<- -10:15
probs<-c(.025,.5,.975)
L<-length(age.range)
LEs<-array(NA,c(L,length(LEs.pnt$pnt),length(probs)))
for(i in 1:L){
  age0<-age.range[i]; ybirth0<-1990-(75+age0)-1900
  cat("Running simulation for age ",age0+75,"and year of birth",ybirth0,"\n")
  results<-elect(model=model, b.covariates=list(age=age0,ybirth=ybirth0),
                 statedistdata=sddata, time.scale.msm="years",h=0.5, age.max=40, S=500,setseed=12345)
  for(j in 1:7){
   for(k in 1:length(probs)){
     LEs[i,j,k]<-quantile(results$sim[,j],probs=probs[k])
   } 
  }
}
x.axis<-c(min(age.range),max(age.range))+75
y.axis<-c(0,20)
plot(x.axis,y.axis,ylab="Life Expectancy",xlab="Age in 1990",type="n",cex.lab=1.5)
# Total LEs:
lines(age.range+75,LEs[,7,1],col="red",lwd=1)
lines(age.range+75,LEs[,7,2],col="red",lwd=3)
lines(age.range+75,LEs[,7,3],col="red",lwd=1)
# Margina LEs in state 1:
lines(age.range+75,LEs[,5,1],col="blue",lwd=1,lty=2)
lines(age.range+75,LEs[,5,2],col="blue",lwd=3,lty=2)
lines(age.range+75,LEs[,5,3],col="blue",lwd=1,lty=2)

