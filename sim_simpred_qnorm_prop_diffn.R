
#-------------------------------------------------------------------------------
#------------------- Simple Method: Simple prediction interval, diff n ---------
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

# Interval: mean plus/minus circa 1.96*sd
# Normal distribution: 95 % of data are in this interval 

predint <- function(dat, alpha=0.05){
  # Proportion of tumors
  prop<-dat$tumor/(dat$tumor+dat$notumor)
  
  mean<-mean(prop)
  sd<-sd(prop)
  
  # Lower boundary of the interval
  lower <- mean - qnorm(p=1-alpha/2, mean=0, sd=1)*sd
  
  # Upper boundary of the interval
  upper <- mean + qnorm(p=1-alpha/2, mean=0, sd=1)*sd
  
  # Lower auf null begrenzen
  if(lower < 0){
    lower <- 0
  }
  
  int <- data.frame(lower=lower, upper=upper)
  
  return(int)
}


#-------------------------------------------------------------------------------

# Simulation
# since the interval is on scale proportion 
# the success count in the actual study is also converted to proportion

predsim <- function(nsim, anzhist, nhistmean, m, exprop, phi) {
  
  predintOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin_diffn(anzhist=anzhist, nhistmean=nhistmean, exprop=exprop, phi=phi)
    datakt <- betbin(anzhist=1, nhist=m, exprop=exprop, phi=phi)
    propakt<-datakt$tumor/(datakt$tumor+datakt$notumor) 
    
    # predint
    PIpredint <- try(predint(dat=dathist))
    
    
    if(class(PIpredint)=="try-error")
    {
      predintOK[i] <- NA
    }#
    
    else
    {
      predintOK[i] <- (PIpredint[1]<= propakt & 
                      propakt<=PIpredint[2])
    }#
  }
  #
  
  # Coverage Probabilities
  cov_predint <- (sum(predintOK, na.rm=TRUE)/nsim)
  
  return(c(predint=cov_predint,
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
  sim_simpred_diffn <- apply(simdat, MARGIN=1, 
                       function(x){
                         predsim(nsim=5000, anzhist=unname(x[1]), m=unname(x[2]),
                                 nhistmean=unname(x[3]), exprop=unname(x[4]), 
                                 phi=unname(x[5]))
                       }) 
)


simulation_simpred_diffn <- data.frame(t(sim_simpred_diffn))

write.csv2(simulation_simpred_diffn, "coverage_simpred_qnorm_prop_diffn.csv")


#-------------------------------------------------------------------------------
