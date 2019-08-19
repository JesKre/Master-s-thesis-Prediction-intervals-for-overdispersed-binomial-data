#-------------------------------------------------------------------------------
#------------------------------- Nelsonphi1 ------------------------------------
#-------------------------------------------------------------------------------

unlink(".RData")
unlink(".Rhistory")

setwd("G:/Dell - Krepel/Ergebnisse")


# Phi bigger/equal 1 (no underdispersion) 

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
  dat <- data.frame(tumor=y, notumor=x)
  return(dat)
}

dat<-betbin_diffn(anzhist = 5, nhistmean = 20, phi=3, exprop = 0.1)
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
  dat <- data.frame(tumor=y, notumor=x)
  return(dat)
}


#-------------------------------------------------------------------------------

### Nelsonphi1 interval

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

nelsonphi1dn <- function(dat, m, alpha=0.025){
        
        # Adjustment
        if(all(dat$tumor==0))
        {
          wmnk <- which.max(dat$tumor+dat$notumor)
          dat$tumor[wmnk] <- dat$tumor[wmnk]+0.5
          dat$notumor[wmnk] <- dat$notumor[wmnk]+0.5
        }
        #
        
        if(all(dat$notumor==0))
        {
          wmnk <- which.max(dat$tumor+dat$notumor)
          dat$tumor[wmnk] <- dat$tumor[wmnk]+0.5
          dat$notumor[wmnk] <- dat$notumor[wmnk]+0.5
        }
        # 
        
        
        fit <- glm(cbind(tumor, notumor)~1, dat, 
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
        
        
        nelsonphi1 <- data.frame(lower=lower, upper=upper)
        
        return(nelsonphi1)
}

nelsonphi1dn(dat=hd[[1]][[1]], alpha=0.05, m=40)
#-------------------------------------------------------------------------------

# Simulation function

# Changes:
# Instead of calculating both historic and actual data together
# with betbin and then seperating it into dathist1 and datakt,
# now dathist is calculated with betbin_diffn
# and datakt is calculated with betbin.

# argument m in predsim function
# m=m instead of m=rowSums(dathist1)[1] in nelsonphi1dn()

predsim <- function(nsim, anzhist, m, nhistmean, exprop, phi) {
        
        nelsonphi1OK <- numeric(length=nsim)
        
        for(i in 1:nsim) { 
                dathist <- betbin_diffn(anzhist=anzhist, nhistmean=nhistmean, exprop=exprop, phi=phi)
                datakt <- betbin(anzhist=1, nhist=m, exprop=exprop, phi=phi)
                
                # Nelsonphi1
                PInelson <- try(nelsonphi1dn(dat=dathist, m=m))
                
                
                if(class(PInelson)=="try-error")
                {
                        nelsonphi1OK[i] <- NA
                }#
                
                else
                {
                        nelsonphi1OK[i] <- (PInelson[1]<= datakt$tumor & 
                                                    datakt$tumor<=PInelson[2])
                }#
        }
        #
        
        # Coverage Probabilities
        cov_nelsonphi1 <- (sum(nelsonphi1OK, na.rm=TRUE)/nsim)
        
        return(c(nelsonphi1=cov_nelsonphi1,
                  
                  nsim=nsim, 
                  anzhist=anzhist, 
                  nhistmean=nhistmean, 
                  m=m,
                  exprop=exprop, 
                  givenphi=phi))
}
#

#-------------------------------------------------------------------------------


simdat <- expand.grid(
        anzhist=c(5, 10, 20, 100),
        m=c(10, 30, 50, 100),
        nhistmean=c(10, 30, 50, 100),
        exprop=c(0.01, 0.05, 0.1, 0.2),
        phi=c(1.01, 1.5, 2, 3))

str(simdat)

system.time(
        sim_nelsonphi1_df <- apply(simdat, MARGIN=1, 
                                   function(x){
                                           predsim(nsim=5000, anzhist=x["anzhist"], m=x["m"],
                                                   nhistmean=x["nhistmean"], exprop=x["exprop"], 
                                                   phi=x["phi"])
                                   }) 
)
#

simulation_nelsonphi1_df <- data.frame(t(sim_nelsonphi1_df))

write.csv2(simulation_nelsonphi1_df, "coverage_nelsonphi1_df.csv")



