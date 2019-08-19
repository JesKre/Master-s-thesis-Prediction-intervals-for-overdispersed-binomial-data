
#-------------------------------------------------------------------------------
#------------------- Simple Method: Range, same n ------------------------------
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
# Lower limit: minimum of proportion of success in hist. data
# Upper limit: maximum of proportion of success in hist. data

rangeprop <- function(dat){
  # Proportion of tumors
  prop<-dat$tumor/(dat$tumor+dat$notumor)
  
  range<-range(prop)
  
  # Lower boundary of the interval
  lower <- range[1]
  
  # Upper boundary of the interval
  upper <- range[2] 
  
  # Lower auf null begrenzen
  if(lower < 0){
    lower <- 0
  }
  
  df<-data.frame(lower=lower, upper=upper)
  
  return(df)
}


#-------------------------------------------------------------------------------

# Simulation
# since the interval is on scale proportion 
# the success count in the actual study is also converted to proportion

predsim <- function(nsim, anzhist, nhist, exprop, phi) {
  
  rangepropOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin(anzhist=anzhist+1, nhist=nhist, exprop=exprop, phi=phi)
    dathist1 <- dathist[1:anzhist,]
    datakt <- dathist[anzhist+1, , drop=FALSE]
    propakt<-datakt$tumor/(datakt$tumor+datakt$notumor)
    
    # rangeprop1
    PIrangeprop <- try(rangeprop(dat=dathist1))
    
    
    if(class(PIrangeprop)=="try-error")
    {
      rangepropOK[i] <- NA
    }#
    
    else
    {
      rangepropOK[i] <- (PIrangeprop[1]<= propakt & 
                         propakt<=PIrangeprop[2])
    }#
  }
  #
  
  # Coverage Probabilities
  cov_rangeprop <- (sum(rangepropOK, na.rm=TRUE)/nsim)
  
  return(c(rangeprop=cov_rangeprop,
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
  sim_rangeprop <- apply(simdat, MARGIN=1, 
                         function(x){
                           predsim(nsim=5000, anzhist=unname(x[1]), 
                                   nhist=unname(x[2]), exprop=unname(x[3]), 
                                   phi=unname(x[4]))
                         }) 
)


simulation_rangeprop <- data.frame(t(sim_rangeprop))

write.csv2(simulation_rangeprop, "coverage_rangeprop_50.csv")



#-------------------------------------------------------------------------------

