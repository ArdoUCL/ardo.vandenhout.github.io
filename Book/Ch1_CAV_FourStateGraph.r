# Chapter 1: Diagram four-state model for CAV data

# Ardo, 2016

# Prelim plot (line widths, text, and colours):
lwd1 <- 3
lwd2 <- 3
states <- c("1","4","2","3")
states.txt0 <- c("  "," ","  ", " ")
states.txt1 <- c("No","Dead","Mild", "Moderate and")
states.txt2 <- c("CAV"," ","CAV", "severe CAV")
kleur<-"black"
kleur.txt<-"blue"

# Parameters (sizes for boxes and arrows):
nodes <- rbind(c(0,85),c(70,30),c(70,85),c(140,85))
rib   <- 30
space <- 5
arrow1 <- 4
arrow2 <- 3

# Frame Plot:
op <- par(mar=c(0,0,0,0))
plot(c(0,180),c(0,100),type="n",ylab="",xlab="", axes=F,main="",asp=.8)

##################
# Plot boxes for four states:
for(i in 1:4){
 a <- nodes[i,]
 lines(c(a[1],a[1]+rib),c(a[2],a[2]),lwd=lwd1,lty=1,col=kleur)
 lines(c(a[1]+rib,a[1]+rib),c(a[2]-rib,a[2]),lwd=lwd1,lty=1,col=kleur)
 lines(c(a[1]+rib,a[1]+rib),c(a[2]-rib,a[2]),lwd=lwd1,lty=1,col=kleur)
 lines(c(a[1],a[1]+rib),c(a[2]-rib,a[2]-rib),lwd=lwd1,lty=1,col=kleur)
 lines(c(a[1],a[1]),c(a[2],a[2]-rib),lwd=lwd1,lty=1,col=kleur)

 # Text in boxes:
 eps <- 3
 text(a[1]+1/2*rib,a[2]-1/4*rib, states[i], cex=2,col=kleur)
 text(a[1]+1/2*rib,a[2]-3/4*rib+2.5*eps, states.txt0[i], cex=1, col=kleur.txt)
 text(a[1]+1/2*rib,a[2]-3/4*rib+eps, states.txt1[i], cex=1, col=kleur.txt)
 text(a[1]+1/2*rib,a[2]-3/4*rib-eps, states.txt2[i], cex=1, col=kleur.txt)
}

########################
# Plot arrows:
# Shortcuts for plotting lines:
a <- nodes[1,]
b <- nodes[2,]
c <- nodes[3,]
d <- nodes[4,]

# q_12:
lines(c(a[1]+rib+space,c[1]-space),c(a[2]-1/3*rib,c[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)
lines(c(c[1]-space-arrow1,c[1]-space),c(c[2]-1/3*rib-arrow2,c[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)
lines(c(c[1]-space-arrow1,c[1]-space),c(c[2]-1/3*rib+arrow2,c[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)

# q_21:
lines(c(a[1]+rib+space,c[1]-space),c(a[2]-2/3*rib,c[2]-2/3*rib),lwd=lwd2,lty=1,col=kleur)
lines(c(a[1]+rib+space,a[1]+rib+space+arrow1),c(a[2]-2/3*rib,a[2]-2/3*rib+arrow2),lwd=lwd2,lty=1,col=kleur)
lines(c(a[1]+rib+space,a[1]+rib+space+arrow1),c(a[2]-2/3*rib,a[2]-2/3*rib-arrow2),lwd=lwd2,lty=1,col=kleur)

# q_14
lines(c(a[1]+1/2*rib,b[1]-space),c(a[2]-rib-space,b[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)
arrow1d <- 2*cos(pi/4)*arrow1
arrow2d <- 2*sin(pi/4)*arrow1
correction <- 1
lines(c(b[1]-space-arrow1d+correction,b[1]-space),c(b[2]-1/3*rib-.8,b[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)
lines(c(b[1]-space,b[1]-space-.8),c(b[2]-1/3*rib,b[2]-1/3*rib+arrow2d),lwd=lwd2,lty=1,col=kleur)

# q_24:
correction <- 2.5
lines(c(c[1]+1/2*rib,c[1]+1/2*rib),c(c[2]-rib-space,b[2]+space),lwd=lwd2,lty=1,col=kleur)
lines(c(c[1]+1/2*rib-arrow1d+correction,c[1]+1/2*rib),c(b[2]+1/8*rib+space,b[2]+space),lwd=lwd2,lty=1,col=kleur)
lines(c(c[1]+1/2*rib,c[1]+1/2*rib+arrow1d-correction),c(b[2]+space,b[2]+1/8*rib+space),lwd=lwd2,lty=1,col=kleur)

# q_23:
lines(c(c[1]+rib+space,d[1]-space),c(c[2]-1/3*rib,d[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)
lines(c(d[1]-space-arrow1,d[1]-space),c(d[2]-1/3*rib-arrow2,d[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)
lines(c(d[1]-space-arrow1,d[1]-space),c(d[2]-1/3*rib+arrow2,d[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)

# q_32:
lines(c(c[1]+rib+space,d[1]-space),c(c[2]-2/3*rib,d[2]-2/3*rib),lwd=lwd2,lty=1,col=kleur)
lines(c(c[1]+rib+space,c[1]+rib+space+arrow1),c(c[2]-2/3*rib,c[2]-2/3*rib+arrow2),lwd=lwd2,lty=1,col=kleur)
lines(c(c[1]+rib+space,c[1]+rib+space+arrow1),c(c[2]-2/3*rib,c[2]-2/3*rib-arrow2),lwd=lwd2,lty=1,col=kleur)

# q_34
lines(c(d[1]+1/2*rib,b[1]+rib+space),c(a[2]-rib-space,b[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)
arrow1d <- 2*cos(pi/4)*arrow1
arrow2d <- 2*sin(pi/4)*arrow1
correction<-2
lines(c(b[1]+rib+space+arrow1d,b[1]+rib+space),c(b[2]-1/3*rib-.8,b[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)
lines(c(b[1]+rib+space+arrow1d-2,b[1]+rib+space),c(b[2]-1/3*rib+arrow2d,b[2]-1/3*rib),lwd=lwd2,lty=1,col=kleur)

