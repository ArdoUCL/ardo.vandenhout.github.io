# A multi-state survival model: msm() and ELECT. Example 2

# Ardo, UCL 2012


# Note: data are not provided for this example. Sorry about that. 

# Preliminaries:
library(msm)
# Adapt the following command to your computer settings:
#source("C:\\mydirectory\\ELECT.r")

# Model fit:
cat("Sample size:"); print(length(table(data$id)))
cat("Frequencies observed state:"); print(table(data$state))
cat("State table:"); print(statetable.msm(state,id,data=data))
q<-0.05
Q<-rbind(c(0,q,q), c(0,0,q),c(0,0,0))
covariates<-as.formula("~age+ybrth+sex+educ")
model<-msm(state~age, subject=id, data=data, center=FALSE, qmatrix=Q,death=TRUE,
              covariates=covariates, censor=c(-1,-2),
              censor.states=list(c(1,2),c(1,2)),method="BFGS",
              control=list(trace=0, REPORT=1,maxit=10000, fnscale=80000))
qnames<-c("q12","q13","q23")
p<-model$estimates; p.se<-sqrt(diag(model$covmat))
print(cbind(q=qnames,p=round(p,3),se=round(p.se,3)),quote=FALSE)

# Life expectancy estimation:
sddata<-data[data$state%in%c(1,2),]
age<- 70 - 78.5
age.max<- 115 - 78.5
ybrth<- 1920 - 1900
educ <-1
# For women:
LEsW<-elect(model=model, b.covariates=list(age=age,ybrth=ybrth,sex=0,educ=educ),
                statedistdata=sddata, h=.5,time.scale="years",
                age.max=age.max, S=1000)
summary.elect(LEsW,probs=c(.025,.975),digits=2)
# For men:
LEsM<-elect(model=model, b.covariates=list(age=age,ybrth=ybrth,sex=1,educ=educ),
                statedistdata=sddata, h=.5, time.scale="years",
                age.max=age.max, S=1000)

# Graphical respresentation of estimated distributions:
lwd<-3; cex.lab<-2
ylab<-c("Density"," "," ","Density"," ", " ")
tekst<-names(LEsW$pnt)
opar<-par(mfrow=c(2,3), mex=0.8,mar=c(5,5,2,1)+.1)
index<-1
for(i in c(1,2,4,5,6,7)){
     plot(c(0,14),c(0,1.5),main="",ylab=ylab[index],xlab="Years",type="n",cex=1,cex.axis=1.5,cex.lab=cex.lab)
     lines(density(LEsW$sim[,i]),col=1,lwd=lwd,cex.lab=cex.lab)
     lines(density(LEsM$sim[,i]),col=2,lwd=lwd,cex.lab=cex.lab)
     text(12,1.4,tekst[i],cex=3,col="blue")
     if(index==1){
       legend(0, 1.5, c("Women", "Men"), col = c(1,2), text.width=4.5,
       text.col = c(1,2), lty=c(1,1),lwd=c(3,3),merge=TRUE,cex=1.5, bg = 'white')
     }
     index<-index+1

 }
opar<-par(mfrow=c(1,1), mex=0.8,mar=c(5,5,2,1)+.1)
