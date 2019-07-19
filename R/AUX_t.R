
#location-scale student-T distribution pdf
dent<-function(x,mu,sigma2,nu){
  z<-(x-mu)/sqrt(sigma2)
  return(1/sqrt(sigma2)*dt(x = z,df = nu))
}

#########################################################################

#location-scale student-T distribution pcf
pent<-function(x,mu,sigma2,nu){
  return(pt((x-mu)/sqrt(sigma2),nu))
}

#########################################################################
