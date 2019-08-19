
#-------------------------------------------------------------------------------
#----------------------------- Simple Method: Whisker interval, same n ---------
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

# Whisker:
# Lower limit is either lower IQR limit minus 1.5*IQR or
# 0 or minimum of data, depending on which is the highest value.
# Upper limit is either upper IQR limit plus 1.5*IQR or
# maximum of data, depending on which is the lowest value.

whiskerprop <- function(dat){
  # Proportion of tumors
  prop<-dat$tumor/(dat$tumor+dat$notumor)
  
  limitsIQR<-quantile(prop, probs=c(0.25, 0.75))
  
  IQR<-IQR(prop)
  
  lower<-max(max(0, unname(limitsIQR[1])-1.5*IQR), min(prop))
  upper<-max(unname(limitsIQR[2])+1.5*IQR, max(prop))
  
  df<-data.frame(lower, upper)
  
  return(df)
}


#-------------------------------------------------------------------------------

# Simulation

predsim <- function(nsim, anzhist, nhist, exprop, phi) {
  
  whiskerpropOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin(anzhist=anzhist+1, nhist=nhist, exprop=exprop, phi=phi)
    dathist1 <- dathist[1:anzhist,]
    datakt <- dathist[anzhist+1, , drop=FALSE]
    propakt<-datakt$tumor/(datakt$tumor+datakt$notumor)
    
    # whiskerprop1
    PIwhiskerprop <- try(whiskerprop(dat=dathist1))
    
    
    if(class(PIwhiskerprop)=="try-error")
    {
      whiskerpropOK[i] <- NA
    }#
    
    else
    {
      whiskerpropOK[i] <- (PIwhiskerprop[1]<= propakt & 
                             propakt<=PIwhiskerprop[2])
    }#
  }
  #
  
  # Coverage Probabilities
  cov_whiskerprop <- (sum(whiskerpropOK, na.rm=TRUE)/nsim)
  
  return(c(whiskerprop=cov_whiskerprop,
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
  sim_whiskerprop <- apply(simdat, MARGIN=1, 
                           function(x){
                             predsim(nsim=5000, anzhist=unname(x[1]), 
                                     nhist=unname(x[2]), exprop=unname(x[3]), 
                                     phi=unname(x[4]))
                           }) 
)


simulation_whiskerprop <- data.frame(t(sim_whiskerprop))

write.csv2(simulation_whiskerprop, "coverage_whiskerprop_neu.csv")


#-------------------------------------------------------------------------------

