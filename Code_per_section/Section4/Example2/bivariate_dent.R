################################################################################
#### enamel carries ############################################################
################################################################################
library(R2WinBUGS)
require(mcmcplots)
library(DTAplots)
#read data
data<-read.table("data.txt",sep=",",header = TRUE)
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

#create list for WinBUGS that contains the data
dat<-list(r=r,N=N)
data.names<-c(names(dat) ,'ns')

#create vector with names of parameters to monitor
#parameter.names <- c( 'theta','rho','sd','sumtpr','sumfpr','spec','Theta','delta','beta','predtpr','predfpf','logitsp','deltasp','Lambda','vartheta','varalpha')
parameter.names <- c( 'predlambda','theta','rho','sd','sumtpf','sumfpf','delta','spec','Theta','beta','predtpr','predfpf','Lambda','vartheta','varalpha','sdtheta','sdalpha')

# defining the directory of WinBUGS
winbugs.dir <- "C:/Users/cb22323/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14/"

# select No of chains, iterations and burnin
n.chain=3
n.iters=70000
n.burn=5000
n.th=1

#create list containing initial value for 3 chains
inits1<-list(  
  list(theta =c(0,-1) , sd=c(1,1) ),
  list(theta =c(1,1) , sd=c(0.6,0.6) ),
  list(theta =c(0.5,0.5) , sd=c(1.5,1.5) )
)
start_time <- Sys.time()
set.seed(455)
#run model through WinBUGS and safe output
model2.sim <- bugs( data.names, inits1, model.file = "model1.txt", parameters = parameter.names,
                    n.chains = n.chain, n.iter = n.iters, n.burnin=n.burn, n.thin=n.th,  bugs.directory = winbugs.dir, debug=FALSE)
end_time <- Sys.time()
end_time - start_time
# print model summary rounded to 3 decimals
print(model2.sim,2)

#### sroc_plot DTA-plots version #################################################################
mu<-cbind(model2.sim$sims.list$theta[,1],-model2.sim$sims.list$theta[,2])
tau<-cbind(model2.sim$sims.list$sd[,1],model2.sim$sims.list$sd[,2])
rho<-matrix(-model2.sim$sims.list$rho,ncol=1)
Summary_Se<-matrix(model2.sim$sims.list$sumtpf,ncol=1)
Summary_Sp<-matrix(model2.sim$sims.list$spec,ncol=1)
se<-plogis(model2.sim$sims.list$delta[,,1])
sp<-plogis(-model2.sim$sims.list$delta[,,2])
Xsroc<-cbind(mu,tau,rho,Summary_Se,Summary_Sp,se,sp)

colse<-colsp<-rep(NA,ns)
for (i in 1:ns){
  colse[i]<-paste("se","[",i,"]",sep="")
  colsp[i]<-paste("sp","[",i,"]",sep="")
}

colnames(Xsroc)<-c("mu[1]","mu[2]","tau[1]","tau[2]","rho","Summary_Se","Summary_Sp",colse,colsp)

svg("sroc_dent.svg", width = 8, height = 8, pointsize = 12)
SROC_rjags(X=Xsroc,n=ns,ref_std = TRUE,cred_region = TRUE,dataset=data[,4:7], predict_region = TRUE,SROC_curve = TRUE)
dev.off()

##################################################################################
#### sroc plot our version #######################################################
senss<-data$TP/data$nA
fpfs<-data$FP/data$nB


#calculate the posterior correlation in the MCMC chains between logit sens and fpf
mcmc1<-model2.sim$sims.list$theta[,1]
mcmc2<-model2.sim$sims.list$theta[,2]
r<-cor(mcmc1,mcmc2)

sens<-median(model2.sim$sims.list$sumtpf)
fpf<-median(model2.sim$sims.list$sumfpf)

results_biv <- data.frame(
  mean0 = c( median(mcmc1) ), #logit sens
  mean1 =c( median(mcmc2)),    #logit spec
  s0 =c(sd(mcmc1) ),          # mcmc chain logit sens sd
  s1 =c( sd(mcmc2) ),         # mcmc chain logit sens sd
  cor = c(r), # NB in WinBUGS you have to use the 'correlation' tool to get this, i.e.
  # this is not a parameter. It's the posterior correlation in the chains between
  # the parameters mean0 and mean 1
  summtpr = c(sens  ), #sens
  summfpr = c(fpf)   #fpf
  
)   
svg("roc_logit_dent.svg", width = 14, height = 7, pointsize = 14) 
par(mfrow=c(1,2))
cex <- 1.2
plot(fpfs, senss,col="white",  pch = 20, xlim = c(0,1), ylim = c(0,1),
     cex = cex, cex.lab = cex, xlab = "1-Specificity (FPF)", ylab = "Sensitivity",
     axes = T,frame.plot =FALSE)

mycol <- rgb(0, 0, 0, max = 255, alpha = 55, names = "myblack")
wno<-data$nA+data$nB
wnos<-wno^(1/1)#sqrt(wno)
symbols(x=fpfs, y=senss, circles=wnos, inches=1/7,
        ann=F, bg=mycol, fg=NULL,add=TRUE)

points(results_biv$summfpr[1], results_biv$summtpr[1], col ="red", pch = 17,cex=cex)

varcov <- array(NA, dim = c(2, 2, 1))

varcov[1,1,1] <- (results_biv$s0[1])^2
varcov[2,2,1] <- (results_biv$s1[1])^2
varcov[1,2,1] <- varcov[2,1,1] <- results_biv$s1[1]*results_biv$s0[1]*results_biv$cor[1]

logitellipse1 <- ellipse(varcov[,,1], centre = c(results_biv$mean0[1], results_biv$mean1[1]), level = 0.95)
ROCellipse1 <- matrix(0, ncol = 2, nrow = nrow(logitellipse1))
ROCellipse1[,1] <- exp(logitellipse1[,2]) / ( 1 + exp(logitellipse1[,2]))
ROCellipse1[,2] <- exp(logitellipse1[,1]) / (1 + exp(logitellipse1[,1]))
lines(ROCellipse1, lty = 2, lwd = 2, col = "red") 

###code for prediction ellipse
mcmc1_pred<-model2.sim$sims.list$predtpr
mcmc2_pred<-model2.sim$sims.list$predfpf
r_pred<-cor(mcmc1_pred,mcmc2_pred)

results_biv_pred <- data.frame(
  mean0 = c( median(mcmc1_pred) ), #logit sens
  mean1 =c( median(mcmc2_pred)),    #logit spec
  s0 =c(sd(mcmc1_pred) ),          # mcmc chain logit sens sd
  s1 =c( sd(mcmc2_pred) ),         # mcmc chain logit sens sd
  cor = c(r_pred) # NB in WinBUGS you have to use the 'correlation' tool to get this, i.e.
  # this is not a parameter. It's the posterior correlation in the chains between
  # the parameters mean0 and mean 1
  
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

###hsroc
fpfsr<-seq(min(fpfs),max(fpfs),by=0.01)[1:100]
nsims<-length(model2.sim$sims.list$sumtpf)
beta<-meanalpha<-meanalphapred<-rep(NA,nsims)
q<-y<-matrix(NA,ncol=length(fpfsr),nrow = nsims)
qpred<-ypred<-matrix(NA,ncol=length(fpfsr),nrow = nsims)
for (i in 1:nsims) {
  beta[i]<-model2.sim$sims.list$beta[i]
  meanalpha[i]<-model2.sim$sims.list$Lambda[i]
  meanalphapred[i]<-model2.sim$sims.list$predlambda[i]
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


lines(hsr$fpfsr,hsr$V2,type="l",col="cornflowerblue",lwd=2)
colrc<-"blue"  
  polygon( c(hsr$fpfsr, rev(hsr$fpfsr)),
           c(hsr$V3, rev(hsr$V4)),
           col=adjustcolor(colrc, alpha.f = 0.15),
           border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
  
  polygon( c(hsrpred$fpfsr, rev(hsrpred$fpfsr)),
           c(hsrpred$V3, rev(hsrpred$V4)),
           col=adjustcolor(colrc, alpha.f = 0.07),
           border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1)) 
  
  points(results_biv$summfpr[1], results_biv$summtpr[1], col ="red", pch = 17,cex=cex)
  
  legend(x = "bottomright", cex = 1,          # Position
         legend = c("Observed (TPF, FPF)","Summary (TPF,FPF)", "95% credible ellipse","95% prediction ellipse","SROC curve"),  # Legend texts
         lty = c(NA,NA,2,2,1),           # Line types
         col = c(mycol, "red","red","blue","cornflowerblue"), # Line colors
         pch = c(20,17, NA,NA,NA),
         lwd = c(NA,NA,2,2,2),bty = "n")   

####logit version 
  #calculate logit-sens and logit-spec for each study
  
  logitsenss<-log(senss/(1-senss))
  logitfpfs<-log(fpfs/(1-fpfs))
  
  plot(logitfpfs, logitsenss,col="white",  pch = 20,xlim=c(-7,3),ylim=c(-6,4),
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
  
  legend(x = "bottomright", cex = 1,          # Position
         legend = c("Observed (logit-TPF, logit-FPF)","Summary (m1,m2)", "95% credible ellipse","95% prediction ellipse","SROC line"),  # Legend texts
         lty = c(NA,NA,2,2,1),           # Line types
         col = c(mycol, "red","red","blue","cornflowerblue"), # Line colors
         pch = c(20,17, NA,NA,NA),
         lwd = c(NA,NA,2,2,2),bty = "n")   
  
  dev.off()