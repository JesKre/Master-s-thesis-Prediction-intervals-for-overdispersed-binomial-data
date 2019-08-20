
#-------------------------------------------------------------------------------
#-------------------- Simple Method: Confidence interval, same n ---------------
#-------------------------------------------------------------------------------


# If both the mean and the standard deviation of the population are unknown and 
# they can only be estimated from a sample of the population, the distribution
# of the sample mean is not a normal distribution but a students t distribution

?qt

qt(0.975, 9)
# [1] 2.262157
qt(0.975, 9)/sqrt(10)
# [1] 0.7153569

qnorm(p=0.975, mean=0, sd=1)
# [1] 1.959964
qnorm(p=0.975, mean=0, sd=1)/sqrt(10) 
# [1] 0.619795
 
# Rees (1989) and Hahn & Meeker (1991) state that a normaldistribution should be
# used for large sample sizes (>30) and the t distribution should be used for
# small sample sizes

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

predsim <- function(nsim, anzhist, nhist, exprop, phi) {
  
  cintOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin(anzhist=anzhist+1, nhist=nhist, exprop=exprop, phi=phi)
    dathist1 <- dathist[1:anzhist,]
    datakt <- dathist[anzhist+1, , drop=FALSE]
    propakt<-datakt$tumor/(datakt$tumor+datakt$notumor) 
    
    # cint1vvvv
    PIcint <- try(cint(dat=dathist1))
    
    
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
  sim_cint <- apply(simdat, MARGIN=1, 
                    function(x){
                      predsim(nsim=5000, anzhist=unname(x[1]), 
                              nhist=unname(x[2]), exprop=unname(x[3]), 
                              phi=unname(x[4]))
                    }) 
)


simulation_cint <- data.frame(t(sim_cint))

write.csv2(simulation_cint, "coverage_cint_t.csv")


#-------------------------------------------------------------------------------
