
#-------------------------------------------------------------------------------
#----------------- Simple Method: Inter quartile range interval, same n --------
#-------------------------------------------------------------------------------


unlink(".RData")
unlink(".Rhistory")

setwd("G:/Dell - Krepel")

#-------------------------------------------------------------------------------

### Data generation 

betbin <- function(anzhist, nhist, exprop, phi) 
{
  asum <- (phi-nhist)/(1-phi) 
  a <- exprop*asum 
  b <- asum-a
  pis <- rbeta(n=anzhist, shape1=a, shape2=b) 
  y <- rbinom(n=anzhist, size=nhist, prob=pis) 
  x <- nhist-y 
  dat <- data.frame(tumor=y, notumor=x)
  return(dat)
}

#-------------------------------------------------------------------------------

# Interval:
# Lower limit: 25 % Quartile
# Upper limit: 75 % Quartile

IQRprop <- function(dat){
  # Proportion of tumors
  prop<-dat$tumor/(dat$tumor+dat$notumor)
  
  IQR<-quantile(prop, probs=c(0.25, 0.75))
  
  df<-data.frame(lower=unname(IQR[1]), upper=unname(IQR[2]))
  
  return(df)
}


#-------------------------------------------------------------------------------

# Simulation

predsim <- function(nsim, anzhist, nhist, exprop, phi) {
  
  IQRpropOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin(anzhist=anzhist+1, nhist=nhist, exprop=exprop, phi=phi)
    dathist1 <- dathist[1:anzhist,]
    datakt <- dathist[anzhist+1, , drop=FALSE]
    propakt<-datakt$tumor/(datakt$tumor+datakt$notumor)
    
    # IQRprop
    PIIQRprop <- try(IQRprop(dat=dathist1))
    
    
    if(class(PIIQRprop)=="try-error")
    {
      IQRpropOK[i] <- NA
    }#
    
    else
    {
      IQRpropOK[i] <- (PIIQRprop[1]<= propakt & 
                         propakt<=PIIQRprop[2])
    }#
  }
  #
  
  # Coverage Probabilities
  cov_IQRprop <- (sum(IQRpropOK, na.rm=TRUE)/nsim)
  
  return(c(IQRprop=cov_IQRprop,
           nsim=nsim, 
           anzhist=anzhist, 
           nhist=nhist, 
           exprop=exprop, 
           givenphi=phi))
}


#-------------------------------------------------------------------------------

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  nhist=c(10, 30, 50, 100),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  phi=c(1.01, 1.5, 2, 3))

system.time(
  sim_IQRprop <- apply(simdat, MARGIN=1, 
                       function(x){
                         predsim(nsim=5000, anzhist=unname(x[1]), 
                                 nhist=unname(x[2]), exprop=unname(x[3]), 
                                 phi=unname(x[4]))
                       }) 
)


simulation_IQRprop <- data.frame(t(sim_IQRprop))

write.csv2(simulation_IQRprop, "coverage_IQRprop_50.csv")



#-------------------------------------------------------------------------------
