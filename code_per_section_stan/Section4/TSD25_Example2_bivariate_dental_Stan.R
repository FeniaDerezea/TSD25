################################################################################
#### TSD25 - Example 2 - Stan ##################################################
################################################################################
#load required libraries
library(R2WinBUGS)
require(mcmcplots)
library(DTAplots)
require(rstan)
library(ellipse)
#read data
data<-read.table("data.txt",sep=",",header = TRUE)
#select only the Radio test from the dataset
data<-data[data$Test=="Radio",1:7]

#calculate number of diseased
data$nA<-data$TP+data$FN
#calculate number of healthy
data$nB<-data$FP+data$TN
# matrix with TP counts in first column and FP counts in second
r<-cbind(data$TP,data$FP)
# matrix with No of healthy in 1st column and No of diseased in 2nd
N<-cbind(data$nA,data$nB)
#calculate No of rows (which is equal to No of studies)
ns<-nrow(r)

#create list for Stan that contains the data
dat<-list(r=r,N=N)
data.names<-c(names(dat) ,'ns')

#create vector with names of parameters to monitor
mymonitoredparamslist <- c( 'predlambda','theta','rho','sd','sumtpf','sumfpf','delta','spec','Theta','beta','predtpr','predfpf','Lambda','vartheta','varalpha','sdtheta','sdalpha')
#and a separate list used for plots, with varaibles separated by index
mymonitoredparamslist2 <- c( 'predlambda','theta[1]','theta[2]','rho','sd[1]','sd[2]',
                             'sumtpf','sumfpf','spec','Theta','beta','predtpr','predfpf',
                             'Lambda','vartheta','varalpha','sdalpha','sdtheta')#,'delta')


# select No of chains, iterations and burnin
n.chain=3
n.iters=3000
n.burn=1000
n.th=1

#create list containing initial value for 3 chains
inits1<-list(  
  list(theta =c(0,-1) , sd=c(1,1), rho=0.5 ),
  list(theta =c(1,1) , sd=c(0.6,0.6), rho=0.1 ),
  list(theta =c(0.5,0.5) , sd=c(1.5,1.5),rho=0.8 )
)

set.seed(455)
#run model through Stan and save output
#model_bivariate.stan contains the Stan code for the bivariate model and should be saved in your current working directory
start.time <- Sys.time()
set.seed(675)
stan_output <- stan(
  file = "model_bivariate.stan",   # Stan model file
  data = data.names,               # Named list of data
  init = inits1,                   # List of initial values (optional)
  pars = mymonitoredparamslist,    # Parameters to monitor
  chains = n.chain,                # Number of chains
  iter = n.iters,                  # Total iterations per chain
  warmup = n.burn,                 # Burn-in (warmup) iterations
  thin = n.th,                     # Thinning interval
  refresh = 100,                   # Progress output frequency
  control = list(adapt_delta=0.99)                
)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

posterior <- as.array(stan_output)

# print model summary rounded to 3 decimals
print(stan_output,digits=3)



##################################################################################
#### sroc plot Figure 5 #######################################################

#calculate TPF and FPF of each study
senss<-data$TP/data$nA
fpfs<-data$FP/data$nB


#calculate the posterior correlation in the MCMC chains between logit sens and fpf
mcmc1<-rstan::extract(stan_output)$theta[,1]
mcmc2<-rstan::extract(stan_output)$theta[,2]
r<-cor(mcmc1,mcmc2)
#save summary sensitivity and FPF
sens<-median(rstan::extract(stan_output)$sumtpf)
fpf<-median(rstan::extract(stan_output)$sumfpf)
#collect required values for the plot in the same object
results_biv <- data.frame(
  mean0 = c( median(mcmc1) ), #logit sens
  mean1 =c( median(mcmc2)),    #logit spec
  s0 =c(sd(mcmc1) ),          # mcmc chain logit sens sd
  s1 =c( sd(mcmc2) ),         # mcmc chain logit sens sd
  cor = c(r), #the posterior correlation in the chains between the parameters mean0 and mean 1
  summtpr = c(sens  ), #sens
  summfpr = c(fpf)   #fpf
  
) 
#the plot will be saved in the working directory unless otherwise specified within the svg command
svg("roc_prob_dent2.svg", width = 7, height = 7, pointsize = 14) 
#par(mfrow=c(1,2))
cex <- 1.2
plot(fpfs, senss,col="white",  pch = 20, xlim = c(0,1), ylim = c(0,1),
     cex = cex, cex.lab = cex, xlab = "1-Specificity (FPF)", ylab = "Sensitivity",
     axes = T,frame.plot =FALSE)
#add observed data with varying size according to sample size
mycol <- rgb(0, 0, 0, max = 255, alpha = 55, names = "myblack")
wno<-data$nA+data$nB
wnos<-wno^(1/1)
symbols(x=fpfs, y=senss, circles=wnos, inches=1/7,
        ann=F, bg=mycol, fg=NULL,add=TRUE)
# add summary TPF and FPF point
points(results_biv$summfpr[1], results_biv$summtpr[1], col ="red", pch = 17,cex=cex)
#create covariance matrix used for the ellipses 
varcov <- array(NA, dim = c(2, 2, 1))

varcov[1,1,1] <- (results_biv$s0[1])^2
varcov[2,2,1] <- (results_biv$s1[1])^2
varcov[1,2,1] <- varcov[2,1,1] <- results_biv$s1[1]*results_biv$s0[1]*results_biv$cor[1]
#create ellipse on the logit scale
logitellipse1 <- ellipse(varcov[,,1], centre = c(results_biv$mean0[1], results_biv$mean1[1]), level = 0.95)
#back-transform ellipse on the probability scale
ROCellipse1 <- matrix(0, ncol = 2, nrow = nrow(logitellipse1))
ROCellipse1[,1] <- exp(logitellipse1[,2]) / ( 1 + exp(logitellipse1[,2]))
ROCellipse1[,2] <- exp(logitellipse1[,1]) / (1 + exp(logitellipse1[,1]))
lines(ROCellipse1, lty = 2, lwd = 2, col = "red") 

###code for prediction ellipse
mcmc1_pred<-rstan::extract(stan_output)$predtpr
mcmc2_pred<-rstan::extract(stan_output)$predfpf
r_pred<-cor(mcmc1_pred,mcmc2_pred)

results_biv_pred <- data.frame(
  mean0 = c( median(mcmc1_pred) ), #logit sens
  mean1 =c( median(mcmc2_pred)),    #logit spec
  s0 =c(sd(mcmc1_pred) ),          # mcmc chain logit sens sd
  s1 =c( sd(mcmc2_pred) ),         # mcmc chain logit sens sd
  cor = c(r_pred) #the posterior correlation in the chains between the parameters mean0 and mean 1
  
)

varcov_pred <- array(NA, dim = c(2, 2, 1))

varcov_pred[1,1,1] <- (results_biv_pred$s0[1])^2
varcov_pred[2,2,1] <- (results_biv_pred$s1[1])^2
varcov_pred[1,2,1] <- varcov_pred[2,1,1] <- results_biv_pred$s1[1]*results_biv_pred$s0[1]*results_biv_pred$cor[1]

logitellipse2 <- ellipse(varcov_pred[,,1], centre = c(results_biv$mean0[1], results_biv$mean1[1]), level = 0.95)
ROCellipse2 <- matrix(0, ncol = 2, nrow = nrow(logitellipse2))
ROCellipse2[,1] <- exp(logitellipse2[,2]) / ( 1 + exp(logitellipse2[,2]))
ROCellipse2[,2] <- exp(logitellipse2[,1]) / (1 + exp(logitellipse2[,1]))
lines(ROCellipse2, lty = 2, lwd = 2, col = "blue") 

### HSROC line ############## 
#create sequence of FPFs within the observed range
fpfsr<-seq(min(fpfs),max(fpfs),by=0.01)[1:100]
#length of posterior sample
nsims<-length(rstan::extract(stan_output)$sumtpf)
#calculate the HSROC line for each mcmc iteration along with predictive values
beta<-meanalpha<-meanalphapred<-rep(NA,nsims)
q<-y<-matrix(NA,ncol=length(fpfsr),nrow = nsims)
qpred<-ypred<-matrix(NA,ncol=length(fpfsr),nrow = nsims)
for (i in 1:nsims) {
  beta[i]<-rstan::extract(stan_output)$beta[i]
  meanalpha[i]<-rstan::extract(stan_output)$Lambda[i]
  meanalphapred[i]<-rstan::extract(stan_output)$predlambda[i]
  for (j in 1:length(fpfsr)){
    q[i,j]<-exp(-beta[i])*logit(fpfsr[j])+exp(-beta[i]/2)*meanalpha[i]
    y[i,j]<-exp(q[i,j])/(exp(q[i,j])+1)
    qpred[i,j]<-exp(-beta[i])*logit(fpfsr[j])+exp(-beta[i]/2)*meanalphapred[i]
    ypred[i,j]<-exp(qpred[i,j])/(exp(qpred[i,j])+1)
  }
}
yl<-matrix(NA,ncol=3,nrow=length(fpfsr))
ylpred<-matrix(NA,ncol=3,nrow=length(fpfsr))
ql<-matrix(NA,ncol=3,nrow=length(fpfsr))
qlpred<-matrix(NA,ncol=3,nrow=length(fpfsr))
for (j in 1:length(fpfsr)){
  #ql[j,]<-quantile(q[,j], c(0.5, 0.025, 0.975))
  yl[j,]<-quantile(y[,j], c(0.5, 0.025, 0.975))
  ylpred[j,]<-quantile(ypred[,j], c(0.5, 0.025, 0.975))
  ql[j,]<-quantile(q[,j], c(0.5, 0.025, 0.975))
  qlpred[j,]<-quantile(qpred[,j], c(0.5, 0.025, 0.975))
}
hsr<-as.data.frame(cbind(fpfsr,yl[,1],yl[,2],yl[,3]))
hsr<-hsr[order(hsr$fpfsr),]
hsrpred<-as.data.frame(cbind(fpfsr,ylpred[,1],ylpred[,2],ylpred[,3]))
hsrpred<-hsrpred[order(hsrpred$fpfsr),]

#add hsroc curve on plot
lines(hsr$fpfsr,hsr$V2,type="l",col="cornflowerblue",lwd=2)
# add 95% cri around hsroc curve
colrc<-"blue"  
polygon( c(hsr$fpfsr, rev(hsr$fpfsr)),
         c(hsr$V3, rev(hsr$V4)),
         col=adjustcolor(colrc, alpha.f = 0.15),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
#add 95% predictive interval around hsroc curve 
polygon( c(hsrpred$fpfsr, rev(hsrpred$fpfsr)),
         c(hsrpred$V3, rev(hsrpred$V4)),
         col=adjustcolor(colrc, alpha.f = 0.07),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1)) 

points(results_biv$summfpr[1], results_biv$summtpr[1], col ="red", pch = 17,cex=cex)

legend(x = "bottomright", cex = 1,          # Position
       legend = c("Observed (TPF, FPF)","Summary (TPF,FPF)", "95% credible region","95% prediction region","SROC curve"),  # Legend texts
       lty = c(NA,NA,2,2,1),           # Line types
       col = c(mycol, "red","red","blue","cornflowerblue"), # Line colors
       pch = c(20,17, NA,NA,NA),
       lwd = c(NA,NA,2,2,2),bty = "n")   
dev.off()

####same plot on logit scale (not necessary for reporting)#####
#calculate logit-sens and logit-spec for each study
svg("roc_logit_dent2.svg", width = 7, height = 7, pointsize = 14)

logitsenss<-log(senss/(1-senss))
logitfpfs<-log(fpfs/(1-fpfs))

plot(logitfpfs, logitsenss,col="white",  pch = 20,xlim=c(-7,5),ylim=c(-7,5),
     cex = cex, cex.lab = cex, xlab = "logit-FPF", ylab = "logit-Sensitivity",
     axes = T,frame.plot =FALSE)


symbols(x=logitfpfs, y=logitsenss, circles=wnos, inches=1/7,
        ann=F, bg=mycol, fg=NULL,add=TRUE)

points(results_biv$mean1[1], results_biv$mean0[1], col ="red", pch = 17,cex=cex)
lines(logitellipse1[,2],logitellipse1[,1], lty = 2, lwd = 2, col = "red") 
lines(logitellipse2[,2],logitellipse2[,1], lty = 2, lwd = 2, col = "blue") 

logitfpfsr<-logit(fpfsr)
lhsr<-as.data.frame(cbind(logitfpfsr,ql[,1],ql[,2],ql[,3]))
lhsr<-lhsr[order(lhsr$logitfpfsr),]
lhsrpred<-as.data.frame(cbind(logitfpfsr,qlpred[,1],qlpred[,2],qlpred[,3]))
lhsrpred<-lhsrpred[order(lhsrpred$logitfpfsr),]
lines(lhsr$logitfpfsr[1:55],lhsr$V2[1:55],type="l",col="cornflowerblue",lwd=2)

polygon( c(lhsr$logitfpfsr[1:55], rev(lhsr$logitfpfsr[1:55])),
         c(lhsr$V3[1:55], rev(lhsr$V4[1:55])),
         col=adjustcolor(colrc, alpha.f = 0.15),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))

polygon( c(lhsrpred$logitfpfsr[1:55], rev(lhsrpred$logitfpfsr[1:55])),
         c(lhsrpred$V3[1:55], rev(lhsrpred$V4[1:55])),
         col=adjustcolor(colrc, alpha.f = 0.07),
         border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1)) 
points(results_biv$mean1[1], results_biv$mean0[1], col ="red", pch = 17,cex=cex)
lines(-7:5,-7:5,lty=3,type="l")
legend(x = "bottomright", cex = 1,          # Position
       legend = c("Observed (logit-TPF, logit-FPF)","Summary (m1,m2)", "95% credible ellipse","95% prediction ellipse","SROC line"),  # Legend texts
       lty = c(NA,NA,2,2,1),           # Line types
       col = c(mycol, "red","red","blue","cornflowerblue"), # Line colors
       pch = c(20,17, NA,NA,NA),
       lwd = c(NA,NA,2,2,2),bty = "n")   

dev.off()

## program ends 