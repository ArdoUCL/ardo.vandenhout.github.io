# Fitting MSSMs to cross-sectional data from the Netherlands: data

# Ardo, UCL 2016

 
# Prelim:
digits<-3

# (Extended) Data taken from van de Kasteele et al (SIM)
# (Only subset so that mortality is no problem):
shift <- -15
age   <- 15:40 + shift
freq1 <- c(114,95,86,88,68,74,70,60,64,53,50,68,51,51,56,59,49,47,57,
           61,68,73,59,62,79,58)
freq2 <- c(7,3,4,10,9,7,10,17,9,9,12,15,20,17,14,16,24,24,22,36,43,48,
           35,72,56,36)
freq3 <- c(2,4,3,1,1,2,3,2,4,3,4,8,4,9,3,6,3,13,5,4,9,14,18,10,17,11)
size  <- freq1+freq2+freq3
dta   <- cbind(age=age,freq1=freq1,freq2=freq2,freq3=freq3,size=size)
n     <- length(age)
   
# Print data:
cat("\nAge transformed by age + shift with shift = ",shift,"\n\n")
cat("Data set:\n")
print(dta)

# Plot data:
# Prelim:
lwd <- 4
pch <- 16
col <- c(1,2,3)
# Plot framwork:
plot(c(age[1]-shift,age[n]-shift),c(0,1),type="n",xlab="Age",
      ylab="Distribution 3 States")
# Plot lines:
for(i in 1:3){
  y <- dta[,i+1]/size
  points(age-shift,y,pch=16,col=col[i])
  lines(age-shift,y,lwd=lwd,col=col[i])
}
# Add legend:
legend(16, .7, c("State 1", "State 2", "State 3"), col = col,
       text.col = "black", lty = 1, lwd=lwd ,
       merge = TRUE, bg = "gray")



