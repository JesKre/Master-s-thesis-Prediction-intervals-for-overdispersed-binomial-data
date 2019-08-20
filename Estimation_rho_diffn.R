
rho<-function(phi, n){
  
  rho<-(phi-1)/(n-1)
  
  outp <- c(rho=rho, 
            n=n,
            phi=phi)
  return(outp)
}

rho(phi=2, n=50)

simdat <- expand.grid(
  n=c(10,30,50,100),
  phi=c(1.01, 1.5, 2, 3))

sim_rho <- apply(simdat, MARGIN=1, 
                               function(x){
                                 rho(n=x["n"], phi=x["phi"])
                               })
rownames(sim_rho)<-c("rho","n","phi")

simulation_rho <- data.frame(t(sim_rho))

simulation_rho

subset(simulation_rho, subset = n==10)
#           rho  n  phi
# 1 0.001111111 10 1.01
# 5 0.055555556 10 1.50
# 9 0.111111111 10 2.00
#13 0.222222222 10 3.00

subset(simulation_rho, subset = n==30)
#             rho  n  phi
# 2  0.0003448276 30 1.01
# 6  0.0172413793 30 1.50
# 10 0.0344827586 30 2.00
# 14 0.0689655172 30 3.00

subset(simulation_rho, subset = n==50)
#             rho  n  phi
# 3  0.0002040816 50 1.01
# 7  0.0102040816 50 1.50
# 11 0.0204081633 50 2.00


subset(simulation_rho, subset = n==100)
#             rho   n  phi
# 4  0.0001010101 100 1.01
# 8  0.0050505051 100 1.50
# 12 0.0101010101 100 2.00
# 16 0.0202020202 100 3.00

