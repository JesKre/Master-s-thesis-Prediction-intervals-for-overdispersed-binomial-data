
#-------------------------------------------------------------------------------
#------------------- Simple Method: Simple prediction interval, same n ---------
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

predsim <- function(nsim, anzhist, nhist, exprop, phi) {
  
  predintOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin(anzhist=anzhist+1, nhist=nhist, exprop=exprop, phi=phi)
    dathist1 <- dathist[1:anzhist,]
    datakt <- dathist[anzhist+1, , drop=FALSE]
    propakt<-datakt$tumor/(datakt$tumor+datakt$notumor) 
    
    # predint
    PIpredint <- try(predint(dat=dathist1))
    
    
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
  sim_simpred <- apply(simdat, MARGIN=1, 
                    function(x){
                      predsim(nsim=5000, anzhist=unname(x[1]), 
                              nhist=unname(x[2]), exprop=unname(x[3]), 
                              phi=unname(x[4]))
                    }) 
)


simulation_simpred <- data.frame(t(sim_simpred))

write.csv2(simulation_simpred, "coverage_simpred_qnorm_prop.csv")


#-------------------------------------------------------------------------------
