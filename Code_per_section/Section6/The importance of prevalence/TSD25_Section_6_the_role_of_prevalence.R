################################################################################
### code for hypothetical decision tree #######################################
################################################################################
#load required libraries
library(MASS)
library(R2WinBUGS)
require(mcmcplots)
################################################################################

##generate fake 2x2 data
ns <- 10

sens<-0.83; fpf<-1-0.78
theta<-rep(NA,2)
theta[1]<-log(sens/(1-sens))
theta[2]<-log(fpf/(1-fpf))

rho<-0.53
var1<-1;var2<-1.7

inmat<-c(var1,rho*sqrt(var1)*sqrt(var2),
         rho*sqrt(var1)*sqrt(var2),var2)

covmat <- matrix(inmat, ncol = 2,nrow=2,byrow = TRUE) 

set.seed(989) 
delta <- mvrnorm(n = ns, mu = theta,  Sigma = covmat)

se<-rep(NA,ns)
fp<-rep(NA,ns)

for (i in 1:ns){
  se[i]<-exp(delta[i,1])/(1+exp(delta[i,1]))
  fp[i]<-exp(delta[i,2])/(1+exp(delta[i,2]))
}

set.seed(344)
N1<-sample(25:50,ns,replace = TRUE)
set.seed(676)
TP<-rbinom(ns, N1, se)

set.seed(665)
N2<-sample(80:500,ns,replace = TRUE)
set.seed(677)
FP<-rbinom(ns, N2, fp)

##fit bivariate model through winbugs
r<-cbind(TP,FP)
# matrix with No of healthy in 1st column and No of diseased in 2nd
N<-cbind(N1,N2)

dat<-list(r=r,N=N)
data.names<-c(names(dat) ,'ns')

#create vector with names of parameters to monitor
parameter.names <- c( 'theta','rho','sd','sumtpr','sumfpr','spec','predtpr','predfpf','pred_tpr','predspec')

# defining the directory of WinBUGS, this should be replaced with the full path to the directory WinBUGS is located
winbugs.dir <- "C:/My_winbugs_directory/WinBUGS14/"

# select No of chains, iterations and burnin
n.chain=3
n.iters=30000
n.burn=5000
n.th=3
#create list containing initial value for 3 chains
inits1<-list(  
  list(theta =c(0,-1) , sd=c(1,1) ),
  list(theta =c(1,1) , sd=c(0.6,0.6) ),
  list(theta =c(0.5,0.5) , sd=c(1.5,1.5) )
)

#run model through WinBUGS and save output
#model1.txt contains the WinBUGS code for the bivariate model and should be saved in your current working directory
model1.sim <- bugs( data.names, inits1, model.file = "model1.txt", parameters = parameter.names,
                    n.chains = n.chain, n.iter = n.iters, n.burnin=n.burn, n.thin=n.th,  bugs.directory = winbugs.dir, debug=FALSE)

# print model summary rounded to 2 decimals
print(model1.sim,2)
#################################################################
##decision tree #################################################
#################################################################

#save needed posterior chains
theta<-model1.sim$sims.list$theta
thetapred1<-model1.sim$sims.list$predtpr
thetapred2<-model1.sim$sims.list$predfpf
thetapred<-cbind(thetapred1,thetapred2)

prev <- seq(0.01, 0.99, length.out =200) # Prevalence
EDT <- -100  # Early detected and treated
LDT <- -200 # Late detected and treated
UFI <- -50 # Unnecessary further investigations
C <- 10 # Cost of screening

nb_noscreen <- prev*LDT   # Net Benefit for no screening

nsims<-nrow(theta)#length of posterior sample
tpf <- fpf <- prev_ce <-rep(NA, nsims)
nb_screen <- inb <-matrix(NA, nsims, 200)
tpfpred <- fpfpred <-  rep(NA, nsims)
nb_screenpred <- inbpred <- matrix(NA, nsims, 200)
for(i in 1:nsims){
  tpf[i] <- exp(theta[i,1])/(1+exp(theta[i,1]))  # TPF 
  fpf[i] <- exp(theta[i,2])/(1+exp(theta[i,2]))  # FPF
  tpfpred[i] <- exp(thetapred[i,1])/(1+exp(thetapred[i,1]))  # predictive TPF 
  fpfpred[i] <- exp(thetapred[i,2])/(1+exp(thetapred[i,2]))  # predictive FPF
  prev_ce[i]<-(UFI*(fpf[i])-C)/(UFI*fpf[i]-(EDT-LDT)*tpf[i]) #prevalence above which test is ce
  for(j in 1:200){
    nb_screen[i,j] <- prev[j]*( tpf[i]*EDT + (1-tpf[i])*LDT) + (1-prev[j])*(fpf[i]*UFI) - C  # Net Benefit of screening
    inb[i,j] <- nb_screen[i,j] - nb_noscreen[j]  # Incremental net benefit (INB) relative to no screening
    ###predictive INB
    nb_screenpred[i,j] <- prev[j]*( tpfpred[i]*EDT + (1-tpfpred[i])*LDT) + (1-prev[j])*(fpfpred[i]*UFI) - C  
    inbpred[i,j] <- nb_screenpred[i,j] - nb_noscreen[j]  
  }
}

# Summaries at each prevalence:
summary_inb <-summary_inbpred<- matrix(NA, 200, 3)
for(j in 1:200){
  #INB (posterior median with 95% CrI)
  summary_inb[j, 1:3] <- quantile(inb[,j], c(0.5, 0.025, 0.975))
  #INB predictive interval
  summary_inbpred[j, 1:3] <- quantile(inbpred[,j], c(0.5, 0.025, 0.975))
}
#median prev_ce and 95% CrI
quantile(prev_ce, c(0.5, 0.025, 0.975))

#### INB vs prevalence plot
# Plot INB (posterior medians):
svg("inb_vs_prev2.svg", width = 8, height = 8, pointsize = 12) # open file to output to
plot(prev, summary_inb[,1], type = "l", ylim=c(-40,80),lwd = 3, col = "firebrick", xlab = "Prevalence", ylab = "Incremental net benefit (INB)", frame.plot =FALSE,cex.lab=1.2, cex.axis=1)

### add 95% cis on plot
coltp<-"red"
  
polygon( c(prev, rev(prev)),
         c(summary_inb[,2], rev(summary_inb[,3])),
         col=adjustcolor(coltp, alpha.f = 0.3),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))


# Add dashed line at INB = 0:

lines(prev, rep(0, 200), type = "l", lty = 3)  

#add predictive intervals

coltp<-"red"
  
polygon( c(prev, rev(prev)),
         c(summary_inbpred[,2], rev(summary_inbpred[,3])),
         col=adjustcolor(coltp, alpha.f = 0.15),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))

legend(x = "topleft",          # Position
       legend = c("Summary INB"),  # Legend texts
       lty = c(1),           # Line types
       col = c("firebrick"), # Line colors
       pch = c(NA),
       lwd = c(3),bty = "n") 
dev.off()

####solve INB eq wrt prevalence ################################################
prevroot<-function(prev,a){
  (quantile(prev*( tpf*EDT + (1-tpf)*LDT) + (1-prev)*(fpf*UFI) - C,c(0.025))-(prev*LDT))-a
}
xmin <- uniroot(prevroot, c(-1E6, 1E6), tol = 0.0001, a = 0)

