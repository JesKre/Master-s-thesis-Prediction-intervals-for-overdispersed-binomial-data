

#-------------------------------------------------------------------------------
#------------------------ calib qBB Interval ----------------------------------
#-------------------------------------------------------------------------------

# Version: Generation of data with rho estimated by phiBBlui() function
#          (according to Lui et al. 2000)


unlink(".RData")
unlink(".Rhistory")

setwd("G:/Dell - Krepel/Ergebnisse")


library(gamlss.dist)


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
  
  ### Intra class correlation according to Lui et al 2000
  
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
  
  
  # Estimated intra class correlation
  rhohatLui <- (BMS-WMS)/(BMS+(mstar-1)*WMS)
  rholim <- (1.001-1)/(nmean-1)
  rhohatlim <- max(rholim, rhohatLui)
  
  sigmacalc <- rhohatlim/(-rhohatlim+1)
  
  pi <- sum(xk)/sum(nk)
  
  return(c(pi=pi, rho=rhohatlim, sigma=sigmacalc, rhohat=rhohatLui))
}


#-------------------------------------------------------------------------------

# Calculation of qBB interval:
# Changes: take sample size of actual study m instead of n of hist. studies

qBB_boot <- function(dat, m, alpha){
  
  # Phi and Pi estimates
  pd <- phiBBLui(data=dat)
  
  # Mu and Sigma
  mu <- as.numeric(pd[1])
  sigma <- as.numeric(pd[3])
  
  # Interval
  l <- qBB(alpha/2, mu=mu, sigma=sigma, bd=m)
  u <- qBB(1-alpha/2, mu=mu, sigma=sigma, bd=m)
  
  
  # Interval as vector
  intcalib <- c(lower=l, upper=u, alpha=alpha)
  
  return(intcalib)
}

qBB_boot(dat=data, m=40, alpha=0.05)

#-------------------------------------------------------------------------------

# install.packages("magrittr")
# install.packages("plyr")
# install.packages("tidyverse")


library(magrittr)
library(plyr)
library(tidyverse)


#-------------------------------------------------------------------------------

# Generation of future studies with same m 
# (a + b) are calculated from rho instead of phi and nk
# same a and b for all studies in data set 

betbin_rho <- function(anzhist, m, exprop, rho) 
{
  asum <- (-rho+1)/rho
  
  a <- exprop*asum 
  b <- asum-a
  pis <- rbeta(n=anzhist, shape1=a, shape2=b) 
  
  y <- rbinom(n=anzhist, size=m, prob=pis) 
  x <- nhist-y 
  
  dat <- data.frame(succ=y, fail=x)
  return(dat)
}

dat<-betbin_rho(anzhist=10, m = 50, exprop=0.1, rho=0.001) 

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

# Generation of bootstrap historical data set with vector of n_k
# from original historical data 

betbin_nvec_rho <- function(anzhist, nhist, exprop, rho) 
{
  asum <- 1/rho-1
  
  a <- exprop*asum 
  b <- asum-a
  pis <- rbeta(n=anzhist, shape1=a, shape2=b) 
  
  y <- rbinom(n=anzhist, size=nhist, prob=pis)
  x <- nhist-y 
  
  dat <- data.frame(succ=y, fail=x)
  return(dat)
}

betbin_nvec_rho(anzhist = 5, nhist= c(50,80,30,55,60), rho=0.003, exprop = 0.1)


#-------------------------------------------------------------------------------

### Parametric Bootstrap function
# Changes: 
# - bootstrap hist. data are generated with bet_bin_nvec_rho and defined rho (instead Phi)
# - bootstrap future study is generated with betbin_rho and defined sample size m and rho

histdat <- function(data, m, k=nrow(data), pi1, rho, nboot=1000){
  
  datastar <- replicate(n=nboot, simplify=FALSE, 
                        expr=betbin_nvec_rho(anzhist=k, exprop=pi1, nhist = data$succ+data$fail,
                                         rho=rho))
  
  futvalue<-unlist(replicate(n=nboot, simplify=FALSE, 
                             expr=betbin_rho(anzhist=1, m=m, exprop=pi1,
                                         rho=rho)[1])) 
  
  # hd contains historical and future BS data
  hd <- list(datastar, futvalue)
  
  return(hd)
  
}

histdat(data, m=40, k=nrow(data), pi1=0.1, rho=0.003, nboot=10)

#-------------------------------------------------------------------------------

### Bisection
# initial search interval (minval, maxval) 
# each step:
# 1) Calculate midpoint of interval
# 2) update upper or lower limit with midpoint 
#    until the distance between the observed cov.prob. and the nominal cov.prob. 
#    is minimized to a satisfactory level 

# Changes: add m=m to function bisec and f 

bisec <- function(f, a1=0.00001, b1=0.3, steps = 15, m=m, tolerance = 0.003, 
                  alpha=0.05, traceplot=TRUE) {
  # nominal coverage probability
  nomcov <- 1-alpha
  
  
  # Startpoints
  # estimated with function defined in bisection(f=..., ) 
  minval <- f(lambda=a1, alpha=alpha, m=m)
  maxval <- f(lambda=b1, alpha=alpha, m=m)
  # Both coverages are higher than 0.95: Return b1 as lambda
  # if difference between calculated cov.prob. and nominal 
  # cov.prob. is positive, b1 stays and a1 is replaced
  if ((minval > 0) && (maxval > 0)) {
    
    val <- c(minval, maxval)
    ab <- c(a1, b1)
    
    # Traceplot
    if(traceplot==TRUE){
      plot(x=ab, y=val, type="p", pch=20, 
           xlab="lambda", ylab=paste("coverage-", nomcov), 
           main="Both higher than 0")
      lines(x=ab, y=val, type="s", col="red")
      abline(a=0, b=0, lty="dashed")
    }
    
    
    return(b1)
    
    
    # Both coverages are smaller than 0.95: Return a1 as lambda       
  } else if ((minval < 0) && (maxval < 0)) {
    
    val <- c(minval, maxval)
    ab <- c(a1, b1)
    
    # Traceplot
    if(traceplot==TRUE){
      plot(x=ab, y=val, type="p", pch=20, 
           xlab="lambda", ylab=paste("coverage-", nomcov),
           main="Both smaller than 0")
      lines(x=ab, y=val, type="s", col="red")
      abline(a=0, b=0, lty="dashed")
    }
    
    
    return(a1)
  }
  
  
  # Open vectors for the estimated lambda and coverage probabilities
  lam <- vector()
  cover <- vector()
  
  for (i in 1:steps) {
    
    # Calculate midpoint
    c1 <- (a1 + b1) / 2 
    
    # Calculation of coverdiff based on c1
    runval <- f(lambda=c1, alpha=alpha, m=m)
    
    # Assigning c1 and runval into the vectors
    lam[i] <- c1
    cover[i] <- runval
    
    # If runval is 0 or
    # If runval is smaller than the tolerance and positive
    # return the value of c1
    if ((runval == 0)  || near(runval, tolerance) || (runval < tolerance) && sign(runval)==1) {
      
      # Traceplot
      if(traceplot==TRUE){
        plot(x=lam, y=cover, type="p", pch=20, 
             xlab="lambda", ylab=paste("coverage-", nomcov), 
             main=paste("Trace with", i, "iterations"))
        lines(x=lam, y=cover, type="s", col="red")
        abline(a=0, b=0, lty="dashed")
      }
      
      
      return(c1)
    }
    
    
    
    # If another iteration is required,
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    if(sign(runval)==1){
      a1 <- c1}
    
    else if(sign(runval)==-1){
      b1 <- c1}
  }
  
  # If the full n iterations are done:
  
  # Open a data.frame
  lamdf <- data.frame("lambda"=lam, "coverdiff"=cover)
  
  # print(cover)
  # 
  # print(lamdf)
  
  # Traceplot
  if(traceplot==TRUE){
    plot(x=lamdf$lambda, y=lamdf$coverdiff, type="p", pch=20, 
         xlab="lambda", ylab=paste("coverage-", nomcov), 
         main=paste("Trace with", i, "iterations"))
    lines(x=lam, y=cover, type="s", col="red")
    abline(a=0, b=0, lty="dashed")
  }
  
  # Take the lambda for which the coverage difference lies between 0 and -tolerance
  if(all(lamdf$coverdiff < 0) & (any(lamdf$coverdiff >= -tolerance) | 
                                 any(near(lamdf$coverdiff, -tolerance)))){
    
    # Order lamdf by increasing lambda
    ldfo <- lamdf[order(lamdf$lambda, decreasing = FALSE),]
    
    # Take the lambda which is closest to 0
    c1 <- ldfo[which.max(ldfo$coverdiff),][1,1]
    return(c1)
  }
  
  # All coverdiff are smaller than -tolerance, return a1
  else if(all(lamdf$coverdiff < -tolerance)){
    return(a1)
  }
  
  # Take the lambda with the smallest positive coverage difference
  else{
    ldf <- filter(lamdf, coverdiff > 0) %>% 
      arrange(desc(lambda))
    c1 <- ldf$lambda[which.min(ldf$coverdiff)]
    return(c1)
  }
}

#-------------------------------------------------------------------------------


### Calibrated interval

# Changes: 
# m = m instead of m = rowSums(data)[1]
# add m=m to covfun and map function

qBB_bisec <- function(data, alpha, k=nrow(data), nboot=1000, m, 
                          a1=0.00001, b1=0.3, steps = 15, tolerance = 0.003, traceplot=FALSE){
  
  # Parameter Estimates
  
  pd <-  phiBBLui(data)
  pi1 <- as.numeric(pd[1])
  rho <- as.numeric(pd[2])
  
  # Bootstrapped historical data and future observation
  hd <- histdat(data=data, k=k, pi1=pi1, m=m, rho=rho, nboot=nboot)
  
  # Coverage Probability (function is needed in the next step)
  covfun <- function(input=hd, lambda, alpha, m=m){
    
    dint <- map(.x=hd[[1]], .f=qBB_boot, m=m, alpha=lambda)  %>% 
      ldply() 
    
    lambdaint <- split(dint, f=factor(dint$alpha)) 
    
    #-----------------------------------------------------------------------
    
    # Is futvalue an element of each interval?
    futinint <- as.data.frame(map(lambdaint, function(x){x$lower<=hd[[2]]&
        x$upper>=hd[[2]]}))
    
    # Coverage probabilities for each lambda
    covprob <- unname(apply(futinint, MARGIN=2, mean))
    
    # Difference between observed coverage and 0.95
    covdiff <- covprob-(1-alpha)
    
    return(covdiff)
  }
  
  # print(hd)
  
  # Generation of the calibrated alpha
  alphacalib <- bisec(f=covfun, a1=a1, m=m, b1=b1, steps=steps, tolerance=tolerance, 
                      alpha=alpha, traceplot=traceplot)
  
  
  
  # Calculation of the interval with alpha calib
  int <- qBB_boot(dat=data, m=m, alpha=alphacalib)
  
  return(int)
 }


#-------------------------------------------------------------------------------

# Simulation function

# Changes: 
# Instead of calculating both historic and actual data together
# with betbin and then seperating it into dathist1 and datakt,
# now dathist is calculated with betbin_diffn_rho
# and datakt is calculated with betbin_rho.

# add argument m to predsim and qBB_bisec function 

predsim <- function(nsim, nboot, anzhist, nhistmean, m, exprop, rho, alpha=0.05) {
  
  qBB_bootOK <- numeric(length=nsim)
  
  
  for(i in 1:nsim) { 
    dathist <- betbin_diffn_rho(anzhist=anzhist, nhistmean=nhistmean, exprop=exprop, rho=rho)
    datakt <- betbin_rho(anzhist=1, m=m, exprop=exprop, rho=rho)
    
    # cat("i", i, "\n")
    # 
    # cat("datakt$tumor", datakt$tumor, "\n")
    
    
    # Alpha calibrated qBBphi1
    PIbootstrap <- try(qBB_bisec(data=dathist, alpha=alpha, m=m, nboot=nboot))
    
    # cat("PIbootstrap")
    # print(PIbootstrap)
    
    
    #---------------------------------------------------------------
    # Logische Abfrage geht die Intervalberechnung?
    
    if(class(PIbootstrap)=="try-error"){
      
      qBB_bootOK[i] <- NA
    }#
    
    else{
      try(datinint <- (PIbootstrap[1]<=  datakt$succ & 
                         datakt$succ<=PIbootstrap[2]))
      
      if(class(datinint)=="try-error"){
        qBB_bootOK[i] <- NA   
        
      }
      
      else{ 
        qBB_bootOK[i] <- datinint
      }
    }#
    
    
    
    
    
    # cat("qBB_bootOK[i]", qBB_bootOK[i])
  }
  #
  
  # Coverage Probabilities
  cov_qBB_boot <- (sum(qBB_bootOK, na.rm=TRUE)/nsim)
  na <- (sum(is.na(qBB_bootOK))/nsim)
  
  outp <- c( qBB_boot=cov_qBB_boot,
             nsim=nsim, 
             anzhist=anzhist, 
             nhistmean=nhistmean, 
             m=m,
             exprop=exprop, 
             rho=rho)
  
  # print(outp)
  # print(str(outp))
  
  return(outp)
}

predsim(nsim=5, nboot=10, anzhist=10, nhistmean=50, m=40, exprop=0.1, rho=0.05, alpha=0.05)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# n=10

# Changes: Instead of defining Phi = 1.01, 1.5, 2, 3 
# Rho is defined depending on Phi and nhistmean,
# see file 'Estimation_rho_diffn.R' 

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(10),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  rho=c(0.001111111, 0.055555556, 0.111111111, 0.22222222))
  
system.time(
  sim_alpha_calib_qBB_10 <- apply(simdat, MARGIN=1, 
                                  function(x){
                                    predsim(nsim=2500, nboot=1000, anzhist=x["anzhist"], 
                                            nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                            rho=x["rho"])
                                  }) 
)
#

# ca 11 Tage auf Sim12b
simulation_alpha_calib_qBB_10 <- data.frame(t(sim_alpha_calib_qBB_10))

write.csv2(simulation_alpha_calib_qBB_10, "qBB_bisec_rho_adapt_n10.csv")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# n=30

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(30),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  rho=c(0.0003448276, 0.0172413793, 0.0344827586, 0.06896552))
  
system.time(
  sim_alpha_calib_qBB_30 <- apply(simdat, MARGIN=1, 
                               function(x){
                                 predsim(nsim=2500, nboot=1000, anzhist=x["anzhist"], 
                                         nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                         rho=x["rho"])
                               }) 
)
#

# ca 11 Tage auf Sim12b
simulation_alpha_calib_qBB_30 <- data.frame(t(sim_alpha_calib_qBB_30))

write.csv2(simulation_alpha_calib_qBB_30, "qBB_bisec_rho_adapt_n30.csv")


#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# n=50

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(50),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  rho=c(0.0002040816,0.0102040816,0.0204081633, 0.04081633))
  # str(simdat)
# str(simdat)

system.time(
  sim_alpha_calib_qBB_50 <- apply(simdat, MARGIN=1, 
                               function(x){
                                 predsim(nsim=2500, nboot=1000, anzhist=x["anzhist"], 
                                         nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                         rho=x["rho"])
                               }) 
)
#

# ca 11 Tage auf Sim12b
simulation_alpha_calib_qBB_50 <- data.frame(t(sim_alpha_calib_qBB_50))

write.csv2(simulation_alpha_calib_qBB_50, "qBB_bisec_rho_adapt_n50.csv")

#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# n=100

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(100),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  rho=c(0.0001010101,0.0050505051,0.0101010101, 0.02020202))
  # str(simdat)
# str(simdat)

system.time(
  sim_alpha_calib_qBB_100 <- apply(simdat, MARGIN=1, 
                               function(x){
                                 predsim(nsim=2500, nboot=1000, anzhist=x["anzhist"], 
                                         nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                         rho=x["rho"])
                               }) 
)
#

# ca 11 Tage auf Sim12b
simulation_alpha_calib_qBB_100 <- data.frame(t(sim_alpha_calib_qBB_100))


write.csv2(simulation_alpha_calib_qBB_100, "qBB_bisec_rho_adapt_n100.csv")

