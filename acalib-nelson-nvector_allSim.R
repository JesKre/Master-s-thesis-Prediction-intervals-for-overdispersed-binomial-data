

#-------------------------------------------------------------------------------
#------------------------ calib Nelsonphi1 Interval ----------------------------
#-------------------------------------------------------------------------------


# Fallzahl für zukünftige Studie m wird festgelegt, 
# verschiedene Fallzahlen n_k für historische Daten

unlink(".RData")
unlink(".Rhistory")

setwd("G:/Dell - Krepel")

#-------------------------------------------------------------------------------

# install.packages("magrittr")
# install.packages("plyr")
# install.packages("tidyverse")


library(magrittr)
library(plyr)
library(tidyverse)


#-------------------------------------------------------------------------------

# Generation of future studies with same m 

betbin <- function(anzhist, nhist, exprop, phi) 
{
  asum <- (phi-nhist)/(1-phi) 
  
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
# and limited to a lower limit of Phi+1 to avoid negative a and b.
# different a and b and therefore different beta distributions and pis 
# for each historical study 

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
  
  dat <- data.frame(succ=y, fail=x, nhist=nhist)
  return(dat)
}


#-------------------------------------------------------------------------------

# Generation of bootstrap historical data set with vector of n_k
# from original historical data 

betbin_nvec <- function(anzhist, nhist, exprop, phi)   
{
  nhist1<-nhist
  
  nhist<-rep(NA, times=length(nhist1))
  nlim <- ceiling(phi)+1
  for (x in 1:length(nhist1)) {nhist[x]<-max(nlim, nhist1[x])}
  
  asum <- (phi-nhist)/(1-phi)
  
  a <- exprop*asum 
  b <- asum-a
  pis <- rbeta(n=anzhist, shape1=a, shape2=b) 
  
  y <- rbinom(n=anzhist, size=nhist, prob=pis)
  x <- nhist-y 
  
  dat <- data.frame(succ=y, fail=x)
  return(dat)
}


#-------------------------------------------------------------------------------


# The dispersion parameter (quasibin) and Pi (probability of success) 
# are estimated from a historical data set (made by betbin function) 

# adjusted for the cases
# all(succ==0); all(fail==0)
# all rows contain either 0 successes or 0 failures

# Same sample size:
# in either case, the first row is chosen
# if all(succ==0) 0.5 successes are added, 0.5 failures are subtracted
# if all(fail==0) 0.5 failures are added, 0.5 successes are subtracted

# Different sample size:
# the row with the largest occuring nk is chosen
# and 0.5 successes and 0.5 failures are added to this row.

pidisp<- function(dat){
  
  # Adjustment
  if(all(dat$succ==0))
  {
    wmnk <- which.max(dat$succ+dat$fail)
    dat$succ[wmnk] <- dat$succ[wmnk]+0.5
    dat$fail[wmnk] <- dat$fail[wmnk]+0.5
  }
  #
  
  if(all(dat$fail==0))
  {
    wmnk <- which.max(dat$succ+dat$fail)
    dat$succ[wmnk] <- dat$succ[wmnk]+0.5
    dat$fail[wmnk] <- dat$fail[wmnk]+0.5
  }
  # 
  
  # dispersionparameter estimated by glm() 
  fit <- glm(cbind(succ, fail)~1, dat,
             family=quasibinomial())
  
  disp <- max(1.001, summary(fit)$dispersion)
  # eta ist das Verhältnis Erfolg/Misserfolg auf der Logarithmus-Skala
  # eta steht im glm-Modell im Intercept
  # Um in die normale Skala zu kommen wird eta eponiert: exp(eta)
  # pi ist die Erfolgswahrscheinlichkeit, diese wird berechnet, indem
  # das Verhältnis durch (1 + das Verhältnis) dividiert wird
  eta <- as.numeric(coef(fit))
  pi <- exp(eta)/(1+exp(eta))
  return(c(disp=disp, pi=pi))
}

#-------------------------------------------------------------------------------

### Parametric Bootstrap function
# Changes: 
# - bootstrap hist. data are generated with bet_bin_nvec_rho and defined rho (instead Phi)
# - bootstrap future study is generated with betbin_rho and defined sample size m and rho


histdat <- function(data, m, k=nrow(data), pi1, disp1, nboot){
  
  datastar <- replicate(n=nboot, simplify=FALSE, 
                        expr=betbin_nvec(anzhist=k, exprop=pi1, nhist = data$nhist,
                                         phi=disp1))
  
  futvalue<-unlist(replicate(n=nboot, simplify=FALSE, 
                             expr=betbin(anzhist=1, nhist=m, exprop=pi1,
                                         phi=disp1)[1])) 
  
  # hd contains historical and future BS data
  hd <- list(datastar, futvalue)
  
  return(hd)
  
}


#-------------------------------------------------------------------------------

### Nelsonphi1 interval

# Adjustment like described for pidisp()

nelsphi1b <- function(dat, alpha, m){

  # Adjustment
  if(all(dat$succ==0))
  {
    wmnk <- which.max(dat$succ+dat$fail)
    dat$succ[wmnk] <- dat$succ[wmnk]+0.5
    dat$fail[wmnk] <- dat$fail[wmnk]+0.5
  }
  #
  
  if(all(dat$fail==0))
  {
    wmnk <- which.max(dat$succ+dat$fail)
    dat$succ[wmnk] <- dat$succ[wmnk]+0.5
    dat$fail[wmnk] <- dat$fail[wmnk]+0.5
  }
  # 
  
  fit <- glm(cbind(succ, fail)~1, dat, 
             family=quasibinomial(link="identity"))
  
  
  # Estimated dispersion parameter
  estphi <- max(1, summary(fit)$dispersion)
  
  # Estimated proportion
  estprop <- coef(fit)[1] 
  names(estprop) <- NULL
  
  # Expected proportion of no tumors in the future sample
  estq <- 1-estprop
  
  # Variance of the estimated proportion
  varestprop <- estphi*estprop*estq*(1/(sum(rowSums(dat))))
  
  # Expected number of tumors in the future study
  esty <- m*estprop
  
  # Standard normal quantile
  z <- qnorm(alpha, 0, 1)
  
  # Lower boundary of the interval
  lower <- esty + z*sqrt(estphi*m*estprop*estq + 
                           (m^2)*varestprop)
  
  # Upper boundary of the interval
  upper <- esty - z*sqrt(estphi*m*estprop*estq +
                           (m^2)*varestprop)
  
  # Auf bzw abrunden 
  lower <- ceiling(lower)
  upper <- floor(upper)
  
  # Lower auf null begrenzen
  if(lower < 0){
    lower <- 0
  }
  
  
  nelsonphi1 <- data.frame(lower=lower, upper=upper, alpha=alpha)
  
  return(nelsonphi1)
}


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
    # If runval is smaller than the tolerance and positve
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
  
  # Take the lambda for wich the coveragedifference lies between 0 and -tolerance
  if(all(lamdf$coverdiff < 0) & (any(lamdf$coverdiff >= -tolerance) | 
                                 any(near(lamdf$coverdiff, -tolerance)))){
    
    # Order lamdf by increasing lambda
    ldfo <- lamdf[order(lamdf$lambda, decreasing = FALSE),]
    
    # Take the lamda which is closest to 0
    c1 <- ldfo[which.max(ldfo$coverdiff),][1,1]
    return(c1)
  }
  
  # All coverdiff are smaller than -tolerance, return a1
  else if(all(lamdf$coverdiff < -tolerance)){
    return(a1)
  }
  
  # Take the lambda with the smallest positive coveragedifference
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

nelson_bisec <- function (data, alpha, k=nrow(data), nboot=1000, m, 
                          a1=0.00001, b1=0.3, steps = 15, tolerance = 0.003, traceplot=FALSE){
  
  # Parameter Estimates
  
  pd <- pidisp(data)
  pi1 <- pd[2]
  disp1 <- pd[1]
  
  # Bootstrapped historical data and future observation
  hd <- histdat(data=data, k=k, pi1=pi1, m=m, disp1=disp1, nboot=nboot)
  
  # Coverage Probability (function is needed in the next step)
  covfun <- function(input=hd, alpha, m=m, lambda){
    
    dint <- map(.x=hd[[1]], .f=nelsphi1b, m=m, alpha=lambda)  %>% 
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
  int <- nelsphi1b(dat=data, m=m, alpha=alphacalib)
  
  return(int)
}


#-------------------------------------------------------------------------------


# Simulation function

# Changes: 
# Instead of calculating both historic and actual data together
# with betbin and then seperating it into dathist1 and datakt,
# now dathist is calculated with betbin_diffn
# and datakt is calculated with betbin.

# add argument m to predsim and nelson_bisec function 

predsim <- function(nsim, nboot, anzhist, nhistmean, m, exprop, phi, alpha=0.05) {
  
  nelson_bootOK <- numeric(length=nsim)
  
  
  for(i in 1:nsim) { 
    dathist <- betbin_diffn(anzhist=anzhist, nhistmean=nhistmean, exprop=exprop, phi=phi)
    datakt <- betbin(anzhist=1, nhist=m, exprop=exprop, phi=phi)
    
    # cat("i", i, "\n")
    # 
    # cat("datakt$tumor", datakt$tumor, "\n")
    
    
    # Alpha calibrated Nelsonphi1
    PIbootstrap <- try(nelson_bisec(data=dathist, alpha=alpha, m=m, nboot=nboot))
    
    # cat("PIbootstrap")
    # print(PIbootstrap)
    
    
    #---------------------------------------------------------------
    # Logische Abfrage geht die Intervalberechnung?
    
    if(class(PIbootstrap)=="try-error"){
      
      nelson_bootOK[i] <- NA
    }#
    
    else{
      try(datinint <- (PIbootstrap[1]<=  datakt$succ & 
                         datakt$succ<=PIbootstrap[2]))
      
      if(class(datinint)=="try-error"){
        nelson_bootOK[i] <- NA   
        
      }
      
      else{ 
        nelson_bootOK[i] <- datinint
      }
    }#
    
    
    
    
    
    # cat("nelson_bootOK[i]", nelson_bootOK[i])
  }
  #
  
  # Coverage Probabilities
  cov_nelson_boot <- (sum(nelson_bootOK, na.rm=TRUE)/nsim)
  na <- (sum(is.na(nelson_bootOK))/nsim)
  
  outp <- c( nelson_boot=cov_nelson_boot,
             nsim=nsim, 
             anzhist=anzhist, 
             nhistmean=nhistmean, 
             m=m,
             exprop=exprop, 
             givenphi=phi)
  
  # print(outp)
  # print(str(outp))
  
  return(outp)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

simdat <- expand.grid(
  anzhist=c(5, 10, 20, 100),
  m=c(10, 30, 50, 100),
  nhistmean=c(10, 30, 50, 100),
  exprop=c(0.01, 0.05, 0.1, 0.2),
  phi=c(1.01, 1.5, 2, 3))
# str(simdat)

system.time(
  sim_alpha_calib_nelson <- apply(simdat, MARGIN=1, 
                                  function(x){
                                    predsim(nsim=2500, nboot=1000, anzhist=x["anzhist"], 
                                            nhistmean=x["nhistmean"], exprop=x["exprop"], m=x["m"], 
                                            phi=x["phi"])
                                  }) 
)
#

# ca 11 Tage auf Sim12b
simulation_alpha_calib_nelson <- data.frame(t(sim_alpha_calib_nelson))

write.csv2(simulation_alpha_calib_nelson, "nelson_bisec_dn_adapt.csv")

