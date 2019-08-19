
#-------------------------------------------------------------------------------
#------------------------ Simple qBB Interval ----------------------------------
#-------------------------------------------------------------------------------


unlink(".RData")
unlink(".Rhistory")

setwd("G:/Dell - Krepel/Ergebnisse")


### Datenerzeugung

# Generation of future studies with same m 
# (a + b) are calculated from rho instead of phi and nk
# same a and b for all studies in data set 

betbin_rho <- function(anzhist, nhist, exprop, rho) 
{
  asum <- (-rho+1)/rho
  
  a <- exprop*asum 
  b <- asum-a
  pis <- rbeta(n=anzhist, shape1=a, shape2=b) 
  
  y <- rbinom(n=anzhist, size=nhist, prob=pis) 
  x <- nhist-y 
  
  dat <- data.frame(succ=y, fail=x)
  return(dat)
}

#-------------------------------------------------------------------------------

# Generation of historical studies 
# with n_k drawn from Poisson distribution with defined mean (nhistmean) 
# (a + b) are calculated from rho instead of phi and nk.
# same a and b for all studies in data set 

betbin_diffn_rho <- function(anzhist, nhistmean, exprop, rho) 
{
  nhist1<-rpois(n=anzhist, lambda=nhistmean)
  nhist<-rep(NA, times=length(nhist1))
  for (x in 1:length(nhist1)) {nhist[x]<-max(2, nhist1[x])}
  asum <- 1/rho-1
  a <- exprop*asum 
  b <- asum-a
  pis <- rbeta(n=anzhist, shape1=a, shape2=b) 
  y <- rbinom(n=anzhist, size=nhist, prob=pis)
  x <- nhist-y 
  dat <- data.frame(succ=y, fail=x)
  return(dat)
}

data<-betbin_diffn_rho(anzhist = 5, nhistmean = 20, rho=0.003, exprop = 0.1)



#-------------------------------------------------------------------------------

##############################################
# Same sample size:
# adjusted for the cases
# all(succ==0); all(fail==0)
# all rows contain either 0 successes or 0 failures
# in either case,
# the first row is chosen
# if all(succ==0) 0.5 successes are added, 0.5 failures are subtracted
# if all(fail==0) 0.5 failures are added, 0.5 successes are subtracted

# Different sample size:
# adjusted for the case that 
# in each row either succ==0 or fail==0
# the row with the largest occuring nk is chosen
# and 0.5 successes and 0.5 failures are added to this row.


phiBBLui <- function(data){
  
  xk <- data$succ
  nk <- data$succ+data$fail
  
  sf0 <- pmin(data$succ, data$fail)
  
  if(all(sf0==0)){
    wmnk <- which.max(nk)
    xk[wmnk] <- xk[wmnk]+0.5
    nk[wmnk] <- nk[wmnk]+1
  }
  
  nmean<-sum(nk)/length(nk)
  k <- nrow(data)
  
  #dcalc <- data.frame(xk=xk, n.xk=nk-xk, nk=nk, sf0=sf0)
  #print(dcalc)
  
  ### Intra class correlation nach Lui et al 2000
  
  # Between mean squared error
  part1B <- sum(xk^2/nk)
  part2B <- (sum(xk)^2/sum(nk))
  part3B <- k-1
  
  BMS <- (part1B-part2B)/part3B
  
  #print(BMS)
  
  # Within mean squared error
  part1W <- sum(xk)
  part2W <- sum(xk^2/nk)
  part3W <- sum(nk-1)
  
  WMS <- (part1W-part2W)/part3W
  
  #print(WMS)
  
  # mstar
  part1m <- sum(nk)^2
  part2m <- sum(nk^2)
  part3m <- (k-1)*sum(nk)
  
  mstar <- (part1m-part2m)/part3m
  
  
  # Estimated intrasclass correlation
  rhohatLui <- (BMS-WMS)/(BMS+(mstar-1)*WMS)
  rholim <- (1.001-1)/(nmean-1)
  rhohatlim <- max(rholim, rhohatLui)
  
  sigmacalc <- rhohatlim/(-rhohatlim+1)
  
  pi <- sum(xk)/sum(nk)
  
  return(c(pi=pi, rho=rhohatlim, sigma=sigmacalc, rhohat=rhohatLui))
}

pd<-phiBBLui(data)

#-------------------------------------------------------------------------------

# Calculation of qBB interval:
# Changes: take sample size of actual study m instead of n of hist. studies

# install.packages("gamlss.dist")
library(gamlss.dist)

qBB_noboot <- function(dat, m, alpha){
  
  # Phi and Pi estimates
  pd <- phiBBLui(data=dat)
  
  # Mu and Sigma
  mu <- as.numeric(pd[1])
  sigma <- as.numeric(pd[3])
  
  # Interval
  l <- qBB(alpha/2, mu=mu, sigma=sigma, bd=m)
  u <- qBB(1-alpha/2, mu=mu, sigma=sigma, bd=m)
  
  
  # Interval as vector
  intcalib <- c(lower=l, upper=u)
  
  return(intcalib)
}

# bd=vector of binomial denominators


#-------------------------------------------------------------------------------

# Simulation function

# Changes: 
# Instead of calculating both historic and actual data together
# with betbin and then seperating it into dathist1 and datakt,
# now dathist is calculated with betbin_diffn_rho
# and datakt is calculated with betbin_rho.

# add argument m to predsim and qBB_noboot function 

predsim <- function(nsim, nboot, anzhist, m, nhistmean, exprop, rho, alpha) {
  
  qBBOK <- numeric(length=nsim)
  
  for(i in 1:nsim) { 
    dathist <- betbin_diffn_rho(anzhist=anzhist, nhistmean=nhistmean, exprop=exprop, rho=rho)
    datakt <- betbin_rho(anzhist=1, nhist=m, exprop=exprop, rho=rho)
    
    # qBB
    PIqBB <- try(qBB_noboot(dat=dathist, m=m, alpha = alpha))
    
    
    if(class(PIqBB)=="try-error")
    {
      qBBOK[i] <- NA
    }#
    
    else
    {
      qBBOK[i] <- (PIqBB[1]<= datakt$succ & 
                     datakt$succ<=PIqBB[2])
    }#
  }
  #
  
  # Coverage Probabilities
  cov_qBB <- (sum(qBBOK, na.rm=TRUE)/nsim)
  
  return(c(qBB=cov_qBB,
           
           nsim=nsim, 
           anzhist=anzhist, 
           nhistmean=nhistmean, 
           m=m,
           exprop=exprop, 
           rho=rho))
}
#

predsim(nsim=5, nboot=10, anzhist=10, nhistmean=50, m=40, exprop=0.1, rho=0.01, alpha=0.05)



#-------------------------------------------------------------------------------
# Simulation n=10

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(10),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  rho=c(0.001111111, 0.055555556,0.111111111, 0.22222222),
  alpha=c(0.05))

system.time(
  sim_qBB_10 <- apply(simdat, MARGIN=1, 
                                  function(x){
                                    predsim(nsim=5000, anzhist=x["anzhist"], 
                                            nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                            rho=x["rho"], alpha=x["alpha"])
                                  }) 
)
#

# ca 11 Tage auf Sim12b
simulation_qBB_10 <- data.frame(t(sim_qBB_10))

write.csv2(simulation_qBB_10, "qBB_10.csv")

#-------------------------------------------------------------------------------
# Simulation n=30

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(30),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  rho=c(0.0003448276,0.0172413793,0.0344827586, 0.06896552),
  alpha=c(0.05))

system.time(
  sim_qBB_30 <- apply(simdat, MARGIN=1, 
                                  function(x){
                                    predsim(nsim=5000, anzhist=x["anzhist"], 
                                            nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                            rho=x["rho"], alpha=x["alpha"])
                                  }) 
)
#

# ca 11 Tage auf Sim12b
simulation_qBB_30 <- data.frame(t(sim_qBB_30))

write.csv2(simulation_qBB_30, "qBB_30.csv")

#-------------------------------------------------------------------------------
# Simulation n=50

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(50),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  rho=c(0.0002040816,0.0102040816,0.0204081633, 0.04081633),
  alpha=c(0.05))
# str(simdat)
# str(simdat)

system.time(
  sim_qBB_50 <- apply(simdat, MARGIN=1, 
                                  function(x){
                                    predsim(nsim=5000, anzhist=x["anzhist"], 
                                            nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                            rho=x["rho"], alpha=x["alpha"])
                                  }) 
)
#

# ca 11 Tage auf Sim12b
simulation_qBB_50 <- data.frame(t(sim_qBB_50))

write.csv2(simulation_qBB_50, "qBB_50.csv")

#-------------------------------------------------------------------------------
# Simulation n=100

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(100),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  rho=c(0.0001010101,0.0050505051,0.0101010101, 0.02020202),
  alpha=c(0.05))
# str(simdat)
# str(simdat)

system.time(
  sim_qBB_100 <- apply(simdat, MARGIN=1, 
                                   function(x){
                                     predsim(nsim=5000, anzhist=x["anzhist"], 
                                             nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                             rho=x["rho"], alpha=x["alpha"])
                                   }) 
)
#

# ca 11 Tage auf Sim12b
simulation_qBB_100 <- data.frame(t(sim_qBB_100))

write.csv2(simulation_qBB_100, "qBB_100.csv")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

