
#-------------------------------------------------------------------------------
#----------------------------- Simple Method: Mean sd interval, diff n ---------
#-------------------------------------------------------------------------------

unlink(".RData")
unlink(".Rhistory")


setwd("G:/Dell - Krepel")

#-------------------------------------------------------------------------------

### Data generation future study

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
#

#-------------------------------------------------------------------------------

### Data generation historical studies

betbin_diffn <- function(anzhist, nhistmean, exprop, phi) 
{
  nhist1<-rpois(n=anzhist, lambda=nhistmean)
  nhist<-rep(NA, times=length(nhist1))
  nlim <- ceiling(phi)+1
  for (x in 1:length(nhist1)) {nhist[x]<-max(nlim, nhist1[x])}
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
# Mean plus/minus standard deviation of the mean 

meansd <- function(dat){
  # Proportion of tumors
  prop<-dat$tumor/(dat$tumor+dat$notumor)
  
  # Expected number of tumors in the the future sample
  mean<-mean(prop)
  
  # expected standard deviation
  sd<-sd(prop) 
  
  # Lower boundary of the interval
  lower <- mean - sd
  
  # Upper boundary of the interval
  upper <- mean + sd
  
  # Lower limited at 0 
  if(lower < 0){
    lower <- 0
  }
  
  msd <- data.frame(lower=lower, upper=upper)
  
  return(msd)
}

#-------------------------------------------------------------------------------

# Simulation

predsim <- function(nsim, anzhist, nhistmean, m, exprop, phi) {
  
  meansdOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin_diffn(anzhist=anzhist, nhistmean=nhistmean, exprop=exprop, phi=phi)
    datakt <- betbin(anzhist=1, nhist=m, exprop=exprop, phi=phi)
    propakt<-datakt$tumor/(datakt$tumor+datakt$notumor)
    
    # meansd1
    PImeansd <- try(meansd(dat=dathist))
    
    
    if(class(PImeansd)=="try-error")
    {
      meansdOK[i] <- NA
    }#
    
    else
    {
      meansdOK[i] <- (PImeansd[1]<= propakt & 
                        propakt<=PImeansd[2])
    }#
  }
  #
  
  # Coverage Probabilities
  cov_meansd <- (sum(meansdOK, na.rm=TRUE)/nsim)
  
  return(c(meansd=cov_meansd,
           nsim=nsim, 
           anzhist=anzhist, 
           nhistmean=nhistmean, 
           exprop=exprop,
           m=m,
           givenphi=phi))
}


#-------------------------------------------------------------------------------

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(10, 30, 50, 100),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  phi=c(1.01, 1.5, 2, 3))


system.time(
  sim_meansd_df <- apply(simdat, MARGIN=1, 
                      function(x){
                        predsim(nsim=5000, anzhist=unname(x[1]), m=unname(x[2]),
                                nhistmean=unname(x[3]), exprop=unname(x[4]), 
                                phi=unname(x[5]))
                      }) 
)


simulation_meansd_df <- data.frame(t(sim_meansd_df))

write.csv2(simulation_meansd_df, "coverage_meansd_prop_diffn.csv")

