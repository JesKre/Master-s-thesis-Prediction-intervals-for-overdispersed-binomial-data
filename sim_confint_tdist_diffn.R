
#-------------------------------------------------------------------------------
#-------------------- Simple Method: Confidence interval, diff n ---------------
#-------------------------------------------------------------------------------


# If both the mean and the standard deviation of the population are unknown and 
# they can only be estimated from a sample of the population, the distribution 
# of the sample mean is not a normal distribution but a students t distribution

# Rees (1989) and Hahn & Meeker (1991) state that a normaldistribution should be
# used for large sample sizes (>30) and the t distribution should be used for
# small sample sizes

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

# Confidence interval of the mean:
# mean +/- (quantile of t distribution)*sd / square root(n)
# according to Hahn & Meeker (1991), page 54 

cint <- function(dat, alpha=0.05){
  
  n<-nrow(dat)
  
  # Proportion of tumors
  prop<-dat$tumor/(dat$tumor+dat$notumor)
  
  mean<-mean(prop)
  sd<-sd(prop)
  
  # Lower boundary of the interval
  lower <- mean - qt(p=1-alpha/2, df=n-1)*sd/sqrt(n)
  
  # Upper boundary of the interval
  upper <- mean + qt(p=1-alpha/2, df=n-1)*sd/sqrt(n) 
  
  
  # Lower limited at 0
  if(lower < 0){
    lower <- 0
  }
  
  int <- data.frame(lower=lower, upper=upper)
  
  return(int)
}


#-------------------------------------------------------------------------------

# Simulation

predsim <- function(nsim, anzhist, nhistmean, m, exprop, phi) {
  
  cintOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin_diffn(anzhist=anzhist, nhistmean=nhistmean, exprop=exprop, phi=phi)
    datakt <- betbin(anzhist=1, nhist=m, exprop=exprop, phi=phi)
    propakt<-datakt$tumor/(datakt$tumor+datakt$notumor) 
    
    # cint
    PIcint <- try(cint(dat=dathist))
    
    
    if(class(PIcint)=="try-error")
    {
      cintOK[i] <- NA
    }#
    
    else
    {
      cintOK[i] <- (PIcint[1]<= propakt & 
                      propakt<=PIcint[2])
    }#
  }
  #
  
  # Coverage Probabilities
  cov_cint <- (sum(cintOK, na.rm=TRUE)/nsim)
  
  return(c(cint=cov_cint,
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

str(simdat)

system.time(
  sim_cint_df <- apply(simdat, MARGIN=1, 
                             function(x){
                               predsim(nsim=5000, anzhist=unname(x[1]), m=unname(x[2]),
                                       nhistmean=unname(x[3]), exprop=unname(x[4]), 
                                       phi=unname(x[5]))
                             }) 
)

simulation_cint_df <- data.frame(t(sim_cint_df))

write.csv2(simulation_cint_df, "coverage_cint_diffn.csv")


#-------------------------------------------------------------------------------
