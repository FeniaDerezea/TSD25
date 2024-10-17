################################################################################
### TSD25 Section 6.5 ##########################################################
################################################################################
#load required libraries
library(MASS)
library(R2WinBUGS)
require(mcmcplots)
library(MCMCpack)
library(faraway)
################################################################################
##generate fake 2x2 data
ns <- 20
meanh<-rep(NA,4)

meanh[1]<-0.81
meanh[2]<-5.2
meanh[3]<-0.19
meanh[4]<-0.43

var1<-0.5;var2<-0.7;var3<-0.47;var4<-0.36
rhom<-0.38
rhos<-0.23
inmat<-c(var1,rhom*sqrt(var1)*sqrt(var2),rhos*sqrt(var1)*sqrt(var3),0,
         rhom*sqrt(var1)*sqrt(var2),var2,0,rhos*sqrt(var2)*sqrt(var4),
         rhos*sqrt(var1)*sqrt(var3),0,var3,0,
         0,rhos*sqrt(var2)*sqrt(var4),0,var4)

covmat <- matrix(inmat, ncol = 4,nrow=4,byrow = TRUE)

set.seed(989) 
delta <- mvrnorm(n = ns, mu = meanh,  Sigma = covmat)

set.seed(344)
N2<-sample(15:50,ns,replace = TRUE)
set.seed(665)
N1<-sample(210:700,ns,replace = TRUE)

set.seed(899)
Tc<-sample(1:10,ns,replace=TRUE)

C<-matrix(NA,ncol=10,nrow=ns)
for (i in 1:ns){
  set.seed(676+i)
  C[i,1:Tc[i]]<-sort(sample(seq(5,150,by=5),Tc[i],replace=FALSE))
}

inprobfp<-matrix(NA,ncol=10,nrow=ns)
inprobse<-matrix(NA,ncol=10,nrow=ns)

for (i in 1:ns){
  for (t in 1:Tc[i]){
    inprobfp[i,t]<-(delta[i,1]-log(C[i,t]))/exp(delta[i,3])
    inprobse[i,t]<-(delta[i,2]-log(C[i,t]))/exp(delta[i,4])
  }
}


FP<-TP<-matrix(NA,ncol=10,nrow=ns)
for (i in 1:ns){
  set.seed(456+i)
  TP[i,1]<-rbinom(1,N2[i],ilogit(inprobse[i,1]))
  set.seed(780+i)
  FP[i,1]<-rbinom(1,N1[i],ilogit(inprobfp[i,1]))
  for (t in 2:Tc[i]){
    set.seed(677+t)
    TP[i,t]<-rbinom(1,TP[i,t-1],ilogit(inprobse[i,t])/ilogit(inprobse[i,t-1]))
    set.seed(782+t)
    FP[i,t]<-rbinom(1,FP[i,t-1],ilogit(inprobfp[i,t])/ilogit(inprobfp[i,t-1]))
  }
}

TP[10,1]<-15
FP[10,1]<-101

FP[2,3]<-NA
FP[3,3:6]<-NA
FP[7,3]<-NA
FP[16,7:9]<-NA

x <- array(NA, c(ns, 2, 10))
for(i in 1:ns){
  for(t in 1:10){
    x[i, 2, t] <- FP[i,t]
    x[i, 1, t] <- TP[i,t]
  }
}

Tc<-matrix(NA,ncol=2,nrow=ns)
for (i in 1:ns){
  Tc[i,2]<-length(na.omit(FP[i,]))
  Tc[i,1]<-length(na.omit(TP[i,]))
}

N<-cbind(N2,N1)

###fit model in winbugs
dataList = list(ns = ns, Tc = Tc, N = N, C = C, x = x)

inits1<-list(  
  list(mean = c(3, 5,log(1), log(1)), 
       sd = c(1, 1, 1, 1),
       rho_mu = 0.8, 
       rho_mu_sigma = 0.2 ) ,
  list(mean = c(0, 6,log(2), log(2)), 
       sd = c(0.1, 0.1, 0.1, 0.1),
       rho_mu = 0.6, 
       rho_mu_sigma = 0.1 ) ,
  list(mean = c(2, 4,log(5), log(5)), 
       sd = c(0.1, 0.2, 0.1, 0.2),
       rho_mu = 0.1, 
       rho_mu_sigma = 0.1 )
)  



mymonitoredparamslist <- c( "mean","resdev","sd","mupred","logspred")#
data.names<-c(names(dataList) )
# defining the directory of WinBUGS, this should be replaced with the full path to the directory WinBUGS is located
winbugs.dir <- "C:/My_winbugs_directory/WinBUGS14/"

##Run the model
#model_log.txt contains the WinBUGS code for the multiple thresholds model and should be saved in your current working directory
set.seed(213)
model1.sim <- bugs(data.names, inits1, model.file = "model_log.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 90000, n.burnin=30000, n.thin=3,  bugs.directory = winbugs.dir, debug=FALSE)

######## save posterior samples ##################################################
nsims<-length(model1.sim$sims.list$resdev)

m_mu <- m_sigma <- matrix(NA, nsims, 2)
m_mu_pred <- m_sigma_pred <- matrix(NA, nsims, 2)
m_mu[,1] <- model1.sim$sims.list$mean[,1]
m_mu[,2] <- model1.sim$sims.list$mean[,2]
m_sigma[,1] <- model1.sim$sims.list$mean[,3]
m_sigma[,2] <- model1.sim$sims.list$mean[,4]

m_mu_pred[,1] <- model1.sim$sims.list$mupred[,1]
m_mu_pred[,2] <- model1.sim$sims.list$mupred[,2]
m_sigma_pred[,1] <- model1.sim$sims.list$logspred[,1]
m_sigma_pred[,2] <- model1.sim$sims.list$logspred[,2]

#### decission tree ############################################################
prev<-0.28
EDT <- -100  # Early detected and treated
LDT <- -200 # Late detected and treated
UFI <- -50 # Unnecessary further investigations
C <- 10 # Cost of screening

nb_noscreen <- prev*LDT   # Net Benefit for no screening


thres<-seq(5,140,length.out=500)
m<-length(thres)

logit_pr <- pr <- array(NA, c(nsims, 2, m )) 
logit_pr_pred <- pr_pred <- array(NA, c(nsims, 2, m )) 
for(k in 1:nsims){
  for(j in 1:2){
    for(t in 1:m){
      # Evaluate logit(TPR) and logit(FPR) at each iteration, k
      logit_pr[k,j,t]<-(m_mu[k,j]-log(thres[t]))/exp(m_sigma[k,j])
      # Evaluate predictive logit(TPR) and logit(FPR) at each iteration, k
      logit_pr_pred[k,j,t]<-(m_mu_pred[k,j]-log(thres[t]))/exp(m_sigma_pred[k,j])
      # Undo the logit transformation:
      pr[k,j,t] <- plogis(logit_pr[k,j,t])
      pr_pred[k,j,t] <- plogis(logit_pr_pred[k,j,t])
    }
  }
}


nb_screen <- inb <- matrix(NA, nsims, m)
nb_screen_pred <- inb_pred <- matrix(NA, nsims, m)

for(i in 1:nsims){
  for(j in 1:m){
    nb_screen[i,j] <- prev*( pr[i,2,j]*EDT + (1-pr[i,2,j])*LDT) + (1-prev)*(pr[i,1,j]*UFI) - C  # Net Benefit of screening
    inb[i,j] <- nb_screen[i,j] - nb_noscreen  # Incremental net benefit (INB) relative to no screening
    #predictive incremental net benefit
    nb_screen_pred[i,j] <- prev*( pr_pred[i,2,j]*EDT + (1-pr_pred[i,2,j])*LDT) + (1-prev)*(pr_pred[i,1,j]*UFI) - C  # Net Benefit of screening
    inb_pred[i,j] <- nb_screen_pred[i,j] - nb_noscreen
    
  }
}

# Summary: INB (posterior median with 95% CrI) for each possible threshold:

summary_inb <- matrix(NA, m, 3)
summary_inb_pred <- matrix(NA, m, 3)
for(j in 1:m){
  summary_inb[j, 1:3] <- quantile(inb[,j], c(0.5, 0.025, 0.975))
  summary_inb_pred[j, 1:3] <- quantile(inb_pred[,j], c(0.5, 0.025, 0.975))
}
## calculating threshold that maximizes the inb with 95% CRi
optthres<-rep(NA,nsims)
for (i in 1:nsims){
  optthres[i]<-thres[which.max(inb[i,])]
}
#summaries
sumopt<-quantile( optthres,c(0.5,0.025,0.975))

## probability of haning highest inb
countinb<-matrix(0, nsims, m)
for(i in 1:nsims){
  for(j in 1:m){
    if (inb[i,j]==max(inb[i,])){
      countinb[i,j]<-1
    }
  }
}

probinb<-rep(NA,m)
for(j in 1:m){
  probinb[j]<-sum(countinb[,j])/nsims
}


#### INB vs threshold plot
# Plot INB (posterior medians):
svg("inb_vs_thres.svg", width = 16, height = 8, pointsize = 14) # open file to output to
par(mfrow=c(1,2))
plot(thres, summary_inb[,1], type = "l", ylim=c(-20,20),xlim=c(5,140),lwd = 3, col = "darkcyan", xlab = "Threshold", ylab = "Incremental net benefit (INB)", frame.plot =FALSE,cex.lab=1.2, cex.axis=1)
#zm<-which(prev>0.1&prev<0.6)
#plot(prev[zm], summary_inb[zm,1], type = "l", ylim=c(-25,35),lwd = 3, col = "firebrick", xlab = "Prevalence", ylab = "Incremental net benefit (INB)", frame.plot =FALSE,cex.lab=1.2, cex.axis=1)

### add 95% cis on plot
coltp<-"darkcyan"
  
polygon( c(thres, rev(thres)),
         c(summary_inb[,2], rev(summary_inb[,3])),
         col=adjustcolor(coltp, alpha.f = 0.3),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))

### add 95% predictive interval

polygon( c(thres, rev(thres)),
         c(summary_inb_pred[,2], rev(summary_inb_pred[,3])),
         col=adjustcolor(coltp, alpha.f = 0.15),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))

#add line INB=0
lines(thres, rep(0, m), type = "l", lty = 3)  
#optimal threshold
points(sumopt[1],summary_inb[which(thres==sumopt[1]),1],pch=18,cex=1.5)
#points(thres[which.max(summary_inb[,1])],max(summary_inb[,1]),pch=18,cex=1.5)

legend(x = "topright",          # Position
       legend = c("Summary INB","Optimal threshold"),  # Legend texts
       lty = c(1,NA),           # Line types
       col = c("darkcyan","black"), # Line colors
       pch = c(NA,18),
       lwd = c(3,NA),bty = "n")    

plot(thres,probinb,type="l",lwd=2,xlab="Threshold",ylab="P(max INB)",frame.plot =FALSE,cex.lab=1.2, cex.axis=1)
lines( rep(thres[which.max(probinb)], m),probinb,col="red",lty = 3,lwd=2)
legend(x = "topright",          # Position
       legend = c("Threshold=25.29"),  # Legend texts
       lty = c(3),           # Line types
       col = c("red"), # Line colors
       pch = c(NA),
       lwd = c(2),bty = "n") 
dev.off()

### dependence on prevalence ########
prev<-seq(0,1,by=0.015)
nb_screen2 <- inb2 <- array(NA, c(nsims, m,length(prev)))

nb_noscreen <- prev*LDT
for(i in 1:nsims){
  for(j in 1:m){
    for (k in 1:length(prev)){
      nb_screen2[i,j,k] <- prev[k]*( pr[i,2,j]*EDT + (1-pr[i,2,j])*LDT) + (1-prev[k])*(pr[i,1,j]*UFI) - C  # Net Benefit of screening
      inb2[i,j,k] <- nb_screen2[i,j,k] - nb_noscreen[k]  # Incremental net benefit (INB) relative to no screening
    }
    
  }
}

##optimal threshold for each prevalence
optthresprev<-matrix(NA,nrow=nsims,ncol=length(prev))
for (i in 1:nsims){
  for (k in 1:length(prev)){
    optthresprev[i,k]<-thres[which.max(inb2[i,,k])]
  }
}

sumthresprev<-matrix(NA,nrow=length(prev),ncol=3)
for (k in 1:length(prev)){
  sumthresprev[k,]<-quantile(optthresprev[,k],c(0.5,0.025,0.975))
}

###optimal threshold plot #########################
svg("optthres_v_prev.svg", width = 8, height = 8, pointsize = 12) # open file to output to
plot(prev,sumthresprev[,1], type = "l", ylim=c(5,140),xlim=c(0,1),lwd = 3, col = "sienna", xlab = "Prevalence", ylab = "Threshold that maximizes INB", frame.plot =FALSE,cex.lab=1.2, cex.axis=1)
coltp<-"chocolate"
  
polygon( c(prev, rev(prev)),
         c(sumthresprev[,2], rev(sumthresprev[,3])),
         col=adjustcolor(coltp, alpha.f = 0.3),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
legend(x = "topright",          # Position
       legend = c("Median optimal threshold"),  # Legend texts
       lty = c(1),           # Line types
       col = c("sienna"), # Line colors
       pch = c(NA),
       lwd = c(3),bty = "n")  
dev.off()

### prevalence estimate below which no thres is cost effective
allneg<-matrix(NA,nrow=nsims,ncol=length(prev))
prevnothres<-rep(NA,nsims)
for (i in 1:nsims){
  for (k in 1:length(prev)){
    allneg[i,k]<-any(inb2[i,,k]>=0)
  
  }
  prevnothres[i]<-prev[length(which(allneg[i,]==FALSE))]
}

#summaries
quantile(prevnothres,c(0.5,0.025,0.975))
