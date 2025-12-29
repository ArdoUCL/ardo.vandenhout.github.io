# Functions for P matrices in Gompertz and Weibull models

# Ardo, UCL 2016


####################################
# Gompertz model:                  #
####################################
Pmatrix.Gomp<-function(t1,t2,lambda,xi,nnodes,h){
    # xi cannot be zero:
    if( any(abs(xi)<eps.xi) ){ for(i in 1:3){ if(abs(xi[i])<eps.xi){ xi[i]<- eps.xi } }  }
    # Construct P-matrix:
    P<-diag(3)
    P[1,1]<-exp( - (lambda[1]/xi[1]*(exp(xi[1]*t2)-exp(xi[1]*t1)) + lambda[2]/xi[2]*(exp(xi[2]*t2)-exp(xi[2]*t1))) )
   
    # Integrand for P[1,2]:
    integrand<-function(u){
      p11<-exp(-(lambda[1]/xi[1]*(exp(xi[1]*u)-exp(xi[1]*t1))+lambda[2]/xi[2]*(exp(xi[2]*u)-exp(xi[2]*t1))))
      q12<-lambda[1]*exp(xi[1]*u)
      p22<- exp(-lambda[3]/xi[3]*(exp(xi[3]*t2)-exp(xi[3]*u)))
      p11*q12*p22
    }
    # Approx using composite Simpson's rule (code from Wikipedia):
    S <- integrand(t1) 
    for(i in seq(1,nnodes,by=2)){ x <- t1+ h*i; S = S+ 4 * integrand(x)} 
    for(i in seq(2,nnodes-1,by=2)){x <- t1 + h*i; S = S+ 2 * integrand(x)} 
    S<-S+integrand(t2)
    P[1,2]<- h * S / 3
    # Check approximation:
    #if(is.na(P[1,2])){P[1,2]<-0}else{if((P[1,1]+P[1,2])>1){P[1,2]<-0}}  #P[1,1:2]<-P[1,1:2]/som}
    
    # Rest of P-matrix:
    P[1,3]<-1-P[1,1]-P[1,2]
    P[2,2]<- exp(-lambda[3]/xi[3]*(exp(xi[3]*t2)-exp(xi[3]*t1)))
    P[2,3]<-1-P[2,2]
    P
}


####################################
# Weibull model:                  #
####################################
Pmatrix.Weib<-function(t1,t2,lambda,tau,nnodes,h){
    # Construct P-matrix:
    P<-diag(3)
    P[1,1]<-exp(-lambda[1]*(t2^tau[1]-t1^tau[1])-lambda[2]*(t2^tau[2]-t1^tau[2]))
   
    # Integrand for P[1,2]:
    integrand<-function(u){
      p11<-exp(-lambda[1]*(u^tau[1]-t1^tau[1])-lambda[2]*(u^tau[2]-t1^tau[2]))
      q12<-lambda[1]*tau[1]*u^(tau[1]-1)
      p22<-exp(-lambda[3]*(t2^tau[3]-u^tau[3]))
      p11*q12*p22
    }
    # Approx using composite Simpson's rule (code from Wikipedia):
    S <- integrand(t1) 
    for(i in seq(1,nnodes,by=2)){ x <- t1+ h*i; S = S+ 4 * integrand(x)} 
    for(i in seq(2,nnodes-1,by=2)){x <- t1 + h*i; S = S+ 2 * integrand(x)} 
    S<-S+integrand(t2)
    P[1,2]<- h * S / 3
    # Check approximation:
    #if(is.na(P[1,2])){P[1,2]<-0}else{if((P[1,1]+P[1,2])>1){P[1,2]<-0}}  #P[1,1:2]<-P[1,1:2]/som}
    
    # Rest of P-matrix:
    P[1,3]<-1-P[1,1]-P[1,2]
    P[2,2]<- exp(-lambda[3]*(t2^tau[3]-t1^tau[3]))
    P[2,3]<-1-P[2,2]
    P
}

