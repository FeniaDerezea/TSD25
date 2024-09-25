################################################################################
###### code for multiple threshold model #######################################
################################################################################
library(R2WinBUGS);library(MASS);library(MCMCpack)
require(mcmcplots)

### OC sensor
data <- read.csv("ACD response_All FOB Gold_20072023.csv",na="NA",header = TRUE)
## bring data in form needed for the model
ns <- length(data$ID)
Tc <- data$T
Tc<-cbind(Tc,Tc)
N <- as.matrix(data[, grepl("N", names(data))])
C <- as.matrix(data[, grepl("C", names(data))])
tp <- as.data.frame(data[, grepl("tp", names(data))])
fp <- as.data.frame(data[, grepl("fp", names(data))])

x <- array(NA, c(ns, 2, max(Tc[,1])))
for(i in 1:ns){
  for(t in 1:max(Tc[,1])){
    x[i, 1, t] <- fp[i,t]
    x[i, 2, t] <- tp[i,t]
  }
}

dataList = list(ns = ns, Tc = Tc, N = N, C = C, x = x)

inits1<-list(  
  list(mean = c(3, 5,log(1), log(1)), 
       sd = c(1, 1, 1, 1) ) ,
  list(mean = c(0, 6,log(2), log(2)), 
       sd = c(0.1, 0.1, 0.1, 0.1) ) ,
  list(mean = c(2, 4,log(5), log(5)), 
       sd = c(0.1, 0.2, 0.1, 0.2))
)  



mymonitoredparamslist <- c( "mean","resdev","sd","mupred","logspred")#
data.names<-c(names(dataList) )
# usual directory
winbugs.dir <- "C:/Users/cb22323/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14/"

##Run the model 
set.seed(213)
model1.sim <- bugs(data.names, inits1, model.file = "model_log_ind.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 150000, n.burnin=50000, n.thin=5,  bugs.directory = winbugs.dir, debug=FALSE)

print(model1.sim,3)

################################################################################
### check linearity
sens<-matrix(NA,ncol=4,nrow=ns)
fpr<-matrix(NA,ncol=4,nrow=ns)
for (i in 1:ns){
  sens[i,1:4]<-(as.matrix(tp[i,])/N[i,2])
  fpr[i,1:4]<-(as.matrix(fp[i,])/N[i,1])
}

logitsens<-matrix(NA,ncol=4,nrow=ns)
logitfpr<-matrix(NA,ncol=4,nrow=ns)
for (i in 1:ns){
  logitsens[i,1:4]<-log(as.matrix(sens[i,])/(1-as.matrix(sens[i,])))
  logitfpr[i,1:4]<-log(as.matrix(fpr[i,])/(1-as.matrix(fpr[i,])))
}
svg("linear.svg", width = 8, height = 7, pointsize = 14) # open file to output to
plot(log(C[,1:4]),logitsens[,1:4],col="white",type="l",ylim=c(-5,5),xlab="log(Threshold)",ylab="logit probability of a positive test result",xlim=c(0,5.4),xaxt="n",frame.plot =FALSE,cex.lab=1.2, cex.axis=1.1,cex.main=1) 
aa<-c(0,1,2,3,4,5)
#axis(side = 1,at=seq(log(thres)[2],max(log(thres)),length.out=8),labels=round(exp(aa),1))
axis(side = 1,at=aa,labels=aa,cex.axis=1.1)
#grid(nx = NULL, ny = NULL,
#     lty = 6,      # Grid line type
#     col = "gray", # Grid line color
#     lwd = 1) 
for(i in 1:ns){
  lines(log(C[i,1:4]),logitsens[i,1:4],col="pink",type="l")
  lines(log(C[i,1:4]),logitsens[i,1:4],col=2,type="p",pch=20)
  lines(log(C[i,1:4]),logitfpr[i,1:4],col="lightblue",type="l")
  lines(log(C[i,1:4]),logitfpr[i,1:4],col="blue",type="p",pch=20)
}

legend(x = "topright",          # Position
       legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)"),  # Legend texts
       lty = c(NA,NA),           # Line types
       col = c("red", "blue"), # Line colors
       pch = c(20, 20),
       lwd = c(NA,NA),bty = "n")  
dev.off()

################################################################################
### summary sens and fpf for different thresholds ##############################
nsims = length(model1.sim$sims.list$resdev)# total number of simulations
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

# sequence of thresholds
threshold<-na.omit(as.vector(C))
thres<-sort(threshold)

thres<-seq(min(thres),max(thres),length.out=800)
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


summ_fpr <- summ_tpr <- matrix(NA, m, 3)
summ_fpr_pred <- summ_tpr_pred <- matrix(NA, m, 3)
for(t in 1:m){
  summ_fpr[t,] <- quantile(pr[,1,t], c(0.5, 0.025, 0.975))
  summ_tpr[t,] <- quantile(pr[,2,t], c(0.5, 0.025, 0.975))
  summ_fpr_pred[t,] <- quantile(pr_pred[,1,t], c(0.5, 0.025, 0.975))
  summ_tpr_pred[t,] <- quantile(pr_pred[,2,t], c(0.5, 0.025, 0.975))
}

#### summary sense spec for each thres #########################################
sumtprfpr<-round(data.frame(thres, summ_tpr, summ_fpr), 2)
write.table(sumtprfpr, "Summary_sens_spec_log_ind.txt", sep=",")

sumtprfpr_pred<-round(data.frame(thres, summ_tpr_pred, summ_fpr_pred), 2)
write.table(sumtprfpr_pred, "Summary_sens_spec_pred_log.txt", sep=",")


#### fit fixed effects version ##################################
inits1<-list(  
  list(mean = c(3, 5,log(1), log(1)) ) ,
  list(mean = c(0, 6,log(2), log(2)) ) ,
  list(mean = c(2, 4,log(5), log(5)))
)  



mymonitoredparamslist <- c( "mean","resdev")#
data.names<-c(names(dataList) )
# usual directory
winbugs.dir <- "C:/Users/cb22323/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14/"

##Run the model 
set.seed(213)
model2.sim <- bugs(data.names, inits1, model.file = "model_log_fe.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 150000, n.burnin=50000, n.thin=5,  bugs.directory = winbugs.dir, debug=FALSE)

print(model2.sim,3)

### summary sens and fpf for different thresholds ##############################
nsims = length(model2.sim$sims.list$resdev)# total number of simulations
m_mu <- m_sigma <- matrix(NA, nsims, 2)
m_mu_pred <- m_sigma_pred <- matrix(NA, nsims, 2)
m_mu[,1] <- model2.sim$sims.list$mean[,1]
m_mu[,2] <- model2.sim$sims.list$mean[,2]
m_sigma[,1] <- model2.sim$sims.list$mean[,3]
m_sigma[,2] <- model2.sim$sims.list$mean[,4]




logit_pr2 <- pr2 <- array(NA, c(nsims, 2, m )) 

for(k in 1:nsims){
  for(j in 1:2){
    for(t in 1:m){
      # Evaluate logit(TPR) and logit(FPR) at each iteration, k
      logit_pr2[k,j,t]<-(m_mu[k,j]-log(thres[t]))/exp(m_sigma[k,j])
      
      # Undo the logit transformation:
      pr2[k,j,t] <- plogis(logit_pr2[k,j,t])
      
    }
  }
}


summ_fpr2 <- summ_tpr2 <- matrix(NA, m, 3)

for(t in 1:m){
  summ_fpr2[t,] <- quantile(pr2[,1,t], c(0.5, 0.025, 0.975))
  summ_tpr2[t,] <- quantile(pr2[,2,t], c(0.5, 0.025, 0.975))
 
}

#### summary sense spec for each thres #########################################
sumtprfpr2<-round(data.frame(thres, summ_tpr2, summ_fpr2), 2)
write.table(sumtprfpr2, "Summary_sens_spec_log_fe.txt", sep=",")

#### plot #######################################
svg("thresplot_fob1.svg", width = 16, height = 7, pointsize = 14) # open file to output to
par(mfrow=c(1,2))
#plot(log(C[1,1:48]),sens[1,1:48],col=2,type="l",xlab="log-threshold",ylab="probability of a positive test result",ylim=c(0,1),xlim = c(0,10),frame.plot =FALSE )
plot(log(C[,1:4]),sens[,1:4],col="white",type="l",xlab=expression("Threshold ("*mu *"g Hb/g)"),ylab="probability of a positive test result",main = "Full model with vague priors",ylim=c(0,1),xlim=c(min(log(thres)),max(log(thres))),frame.plot =FALSE,xaxt="n",cex.lab=1.2, cex.axis=1,cex.main=1)
#aa<-seq(log(thres)[2],max(log(thres)),length.out=8)
aa<-c(2,5,10,20,50,100,150)
#axis(side = 1,at=seq(log(thres)[2],max(log(thres)),length.out=8),labels=round(exp(aa),1))
axis(side = 1,at=log(aa),labels=aa)
#grid(nx = NULL, ny = NULL,
#     lty = 6,      # Grid line type
#     col = "gray", # Grid line color
#     lwd = 1) 



mycol <- rgb(255, 0, 0, max = 255, alpha = 55, names = "myred")
mycol2 <- rgb(0, 0, 255, max = 255, alpha = 55, names = "myblue")



xno<-NULL
yno<-NULL
wno<-NULL
yno2<-NULL
wno2<-NULL

for(i in 1:ns){
  lines(log(C[i,1:4]),sens[i,1:4],col="darkgray",type="l")
  newx<-log(C[i,1:4])
  newy<-sens[i,1:4]
  newy2<-fpr[i,1:4]
  neww<-rep(N[i,2],4)
  neww2<-rep(N[i,1],4)
  xno<-c(xno,newx)
  yno<-c(yno,newy)
  wno<-c(wno,neww)
  yno2<-c(yno2,newy2)
  wno2<-c(wno2,neww2)
  #lines(log(C[i,1:10]),sens[i,1:10],col=2,type="p",pch=20)
  lines(log(C[i,1:4]),fpr[i,1:4],col="darkgray",type="l")
  #lines(log(C[i,1:10]),fpr[i,1:10],col="blue",type="p",pch=20)
}

wnos<-wno^(1/1)#sqrt(wno)
symbols(x=xno, y=yno, circles=wnos, inches=1/7,
        ann=F, bg=mycol, fg=NULL,add=TRUE)
wnos2<-wno^(1/1)#sqrt(wno2)
symbols(x=xno, y=yno2, circles=wnos2, inches=1/7,
        ann=F, bg=mycol2, fg=NULL,add=TRUE)

lines(log(thres), summ_fpr[,1], type = "l",col="blue",lwd=3)
lines(log(thres), summ_tpr[,1], type = "l",col="red",lwd=3)

### add 95% cis on plot
coltp<-"red"
  colfp<-"blue"  
    polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
             c(summ_tpr[,2], rev(summ_tpr[,3])),
             col=adjustcolor(coltp, alpha.f = 0.3),
             border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
    
    polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
             c(summ_fpr[,2], rev(summ_fpr[,3])),
             col=adjustcolor(colfp, alpha.f = 0.2),
             border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
    
    ### predictive intervals
    polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
             c(summ_tpr_pred[,2], rev(summ_tpr_pred[,3])),
             col=adjustcolor(coltp, alpha.f = 0.1),
             border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
    
    polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
             c(summ_fpr_pred[,2], rev(summ_fpr_pred[,3])),
             col=adjustcolor(colfp, alpha.f = 0.05),
             border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))   
    
    legend(x = "topright",          # Position
           legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)","Summary TPF", "Summary FPF"),  # Legend texts
           lty = c(NA,NA,1, 1),           # Line types
           col = c("red", "blue","red","blue"), # Line colors
           pch = c(20, 20, NA, NA),
           lwd = c(NA,NA,2,2),bty = "n")    
    
    #### fixed effects
    plot(log(C[,1:4]),sens[,1:4],col="white",type="l",xlab=expression("Threshold ("*mu *"g Hb/g)"),ylab="probability of a positive test result",main = "Fixed effects version",ylim=c(0,1),xlim=c(min(log(thres)),max(log(thres))),frame.plot =FALSE,xaxt="n",cex.lab=1.2, cex.axis=1,cex.main=1)
    #aa<-seq(log(thres)[2],max(log(thres)),length.out=8)
    aa<-c(2,5,10,20,50,100,150)
    #axis(side = 1,at=seq(log(thres)[2],max(log(thres)),length.out=8),labels=round(exp(aa),1))
    axis(side = 1,at=log(aa),labels=aa)
    #grid(nx = NULL, ny = NULL,
    #     lty = 6,      # Grid line type
    #     col = "gray", # Grid line color
    #     lwd = 1) 
    xno<-NULL
    yno<-NULL
    wno<-NULL
    yno2<-NULL
    wno2<-NULL
    for(i in 1:ns){
      lines(log(C[i,1:4]),sens[i,1:4],col="darkgray",type="l")
      newx<-log(C[i,1:4])
      newy<-sens[i,1:4]
      newy2<-fpr[i,1:4]
      neww<-rep(N[i,2],4)
      neww2<-rep(N[i,1],4)
      xno<-c(xno,newx)
      yno<-c(yno,newy)
      wno<-c(wno,neww)
      yno2<-c(yno2,newy2)
      wno2<-c(wno2,neww2)
      #lines(log(C[i,1:10]),sens[i,1:10],col=2,type="p",pch=20)
      lines(log(C[i,1:4]),fpr[i,1:4],col="darkgray",type="l")
      #lines(log(C[i,1:10]),fpr[i,1:10],col="blue",type="p",pch=20)
    }
    wnos<-wno^(1/1)#sqrt(wno)
    symbols(x=xno, y=yno, circles=wnos, inches=1/7,
            ann=F, bg=mycol, fg=NULL,add=TRUE)
    wnos2<-wno2^(1/1)#sqrt(wno2)
    symbols(x=xno, y=yno2, circles=wnos2, inches=1/7,
            ann=F, bg=mycol2, fg=NULL,add=TRUE)
    
    lines(log(thres), summ_fpr2[,1], type = "l",col="blue",lwd=3)
    lines(log(thres), summ_tpr2[,1], type = "l",col="red",lwd=3)
    
    coltp<-"red"
      colfp<-"blue"  
        polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
                 c(summ_tpr2[,2], rev(summ_tpr2[,3])),
                 col=adjustcolor(coltp, alpha.f = 0.3),
                 border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
        
        polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
                 c(summ_fpr2[,2], rev(summ_fpr2[,3])),
                 col=adjustcolor(colfp, alpha.f = 0.2),
                 border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
        
      
        
        
        legend(x = "topright",          # Position
               legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)","Summary TPF", "Summary FPF"),  # Legend texts
               lty = c(NA,NA,1, 1),           # Line types
               col = c("red", "blue","red","blue"), # Line colors
               pch = c(20, 20, NA, NA),
               lwd = c(NA,NA,2,2),bty = "n")    
        
        dev.off() # turn off device   
        
 #### ##########################################################################
#### informative priors  version ###############################################
        inits1<-list(  
          list(mean = c(3, 5,log(1), log(1)), 
               sd = c(1, 1, 1, 1) ) ,
          list(mean = c(0, 6,log(2), log(2)), 
               sd = c(0.1, 0.1, 0.1, 0.1) ) ,
          list(mean = c(2, 4,log(5), log(5)), 
               sd = c(0.1, 0.2, 0.1, 0.2))
        )  
        
        
        
        mymonitoredparamslist <- c( "mean","resdev","sd","mupred","logspred")#
        data.names<-c(names(dataList) )
        # usual directory
        winbugs.dir <- "C:/Users/cb22323/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14/"
        
        ##Run the model 
        set.seed(213)
        model3.sim <- bugs(data.names, inits1, model.file = "model_log_inf.txt", parameters = mymonitoredparamslist ,
                           n.chains = 3, n.iter = 3000000, n.burnin=1000000, n.thin=20,  bugs.directory = winbugs.dir, debug=FALSE)
        
        mcmcplot(model3.sim)
        print(model3.sim,3)
 
        ### summary sens and fpf for different thresholds ##############################
        nsims = length(model3.sim$sims.list$resdev)# total number of simulations
        m_mu <- m_sigma <- matrix(NA, nsims, 2)
        m_mu_pred <- m_sigma_pred <- matrix(NA, nsims, 2)
        m_mu[,1] <- model3.sim$sims.list$mean[,1]
        m_mu[,2] <- model3.sim$sims.list$mean[,2]
        m_sigma[,1] <- model3.sim$sims.list$mean[,3]
        m_sigma[,2] <- model3.sim$sims.list$mean[,4]
        
        m_mu_pred[,1] <- model3.sim$sims.list$mupred[,1]
        m_mu_pred[,2] <- model3.sim$sims.list$mupred[,2]
        m_sigma_pred[,1] <- model3.sim$sims.list$logspred[,1]
        m_sigma_pred[,2] <- model3.sim$sims.list$logspred[,2]
        
      
        
        logit_pr3 <- pr3 <- array(NA, c(nsims, 2, m )) 
        logit_pr_pred3 <- pr_pred3 <- array(NA, c(nsims, 2, m )) 
        for(k in 1:nsims){
          for(j in 1:2){
            for(t in 1:m){
              # Evaluate logit(TPR) and logit(FPR) at each iteration, k
              logit_pr3[k,j,t]<-(m_mu[k,j]-log(thres[t]))/exp(m_sigma[k,j])
              # Evaluate predictive logit(TPR) and logit(FPR) at each iteration, k
              logit_pr_pred3[k,j,t]<-(m_mu_pred[k,j]-log(thres[t]))/exp(m_sigma_pred[k,j])
              # Undo the logit transformation:
              pr3[k,j,t] <- plogis(logit_pr3[k,j,t])
              pr_pred3[k,j,t] <- plogis(logit_pr_pred3[k,j,t])
            }
          }
        }
        
        
        summ_fpr3 <- summ_tpr3 <- matrix(NA, m, 3)
        summ_fpr_pred3 <- summ_tpr_pred3 <- matrix(NA, m, 3)
        for(t in 1:m){
          summ_fpr3[t,] <- quantile(pr3[,1,t], c(0.5, 0.025, 0.975))
          summ_tpr3[t,] <- quantile(pr3[,2,t], c(0.5, 0.025, 0.975))
          summ_fpr_pred3[t,] <- quantile(pr_pred3[,1,t], c(0.5, 0.025, 0.975))
          summ_tpr_pred3[t,] <- quantile(pr_pred3[,2,t], c(0.5, 0.025, 0.975))
        }
        
        #### summary sense spec for each thres #########################################
        sumtprfpr3<-round(data.frame(thres, summ_tpr3, summ_fpr3), 2)
        write.table(sumtprfpr3, "Summary_sens_spec_log_inf.txt", sep=",")
        
        sumtprfpr_pred3<-round(data.frame(thres, summ_tpr_pred3, summ_fpr_pred3), 2)
        write.table(sumtprfpr_pred3, "Summary_sens_spec_pred_log.txt", sep=",")
 
        
        #### plot #######################################
        svg("thresplot_fob2.svg", width = 9, height = 7, pointsize = 12) # open file to output to
        par(mfrow=c(1,1))
        #plot(log(C[1,1:48]),sens[1,1:48],col=2,type="l",xlab="log-threshold",ylab="probability of a positive test result",ylim=c(0,1),xlim = c(0,10),frame.plot =FALSE )
        plot(log(C[,1:4]),sens[,1:4],col="white",type="l",xlab=expression("Threshold ("*mu *"g Hb/g)"),ylab="probability of a positive test result",main = "Full model with informative priors",ylim=c(0,1),xlim=c(min(log(thres)),max(log(thres))),frame.plot =FALSE,xaxt="n",cex.lab=1.2, cex.axis=1,cex.main=1)
        #aa<-seq(log(thres)[2],max(log(thres)),length.out=8)
        aa<-c(2,5,10,20,50,100,150)
        #axis(side = 1,at=seq(log(thres)[2],max(log(thres)),length.out=8),labels=round(exp(aa),1))
        axis(side = 1,at=log(aa),labels=aa)
        #grid(nx = NULL, ny = NULL,
        #     lty = 6,      # Grid line type
        #     col = "gray", # Grid line color
        #     lwd = 1) 
        
        
        
        mycol <- rgb(255, 0, 0, max = 255, alpha = 55, names = "myred")
        mycol2 <- rgb(0, 0, 255, max = 255, alpha = 55, names = "myblue")
        
        
        
        xno<-NULL
        yno<-NULL
        wno<-NULL
        yno2<-NULL
        wno2<-NULL
        
        for(i in 1:ns){
          lines(log(C[i,1:4]),sens[i,1:4],col="darkgray",type="l")
          newx<-log(C[i,1:4])
          newy<-sens[i,1:4]
          newy2<-fpr[i,1:4]
          neww<-rep(N[i,2],4)
          neww2<-rep(N[i,1],4)
          xno<-c(xno,newx)
          yno<-c(yno,newy)
          wno<-c(wno,neww)
          yno2<-c(yno2,newy2)
          wno2<-c(wno2,neww2)
          #lines(log(C[i,1:10]),sens[i,1:10],col=2,type="p",pch=20)
          lines(log(C[i,1:4]),fpr[i,1:4],col="darkgray",type="l")
          #lines(log(C[i,1:10]),fpr[i,1:10],col="blue",type="p",pch=20)
        }
        
        wnos<-wno^(1/1)#sqrt(wno)
        symbols(x=xno, y=yno, circles=wnos, inches=1/7,
                ann=F, bg=mycol, fg=NULL,add=TRUE)
        wnos2<-wno^(1/1)#sqrt(wno2)
        symbols(x=xno, y=yno2, circles=wnos2, inches=1/7,
                ann=F, bg=mycol2, fg=NULL,add=TRUE)
        
        lines(log(thres), summ_fpr3[,1], type = "l",col="blue",lwd=3)
        lines(log(thres), summ_tpr3[,1], type = "l",col="red",lwd=3)
        
        ### add 95% cis on plot
        coltp<-"red"
          colfp<-"blue"  
            polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
                     c(summ_tpr3[,2], rev(summ_tpr3[,3])),
                     col=adjustcolor(coltp, alpha.f = 0.3),
                     border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
            
            polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
                     c(summ_fpr3[,2], rev(summ_fpr3[,3])),
                     col=adjustcolor(colfp, alpha.f = 0.2),
                     border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
            
            ### predictive intervals
            polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
                     c(summ_tpr_pred3[,2], rev(summ_tpr_pred3[,3])),
                     col=adjustcolor(coltp, alpha.f = 0.1),
                     border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
            
            polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
                     c(summ_fpr_pred3[,2], rev(summ_fpr_pred3[,3])),
                     col=adjustcolor(colfp, alpha.f = 0.05),
                     border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))   
            
            legend(x = "topright",          # Position
                   legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)","Summary TPF", "Summary FPF"),  # Legend texts
                   lty = c(NA,NA,1, 1),           # Line types
                   col = c("red", "blue","red","blue"), # Line colors
                   pch = c(20, 20, NA, NA),
                   lwd = c(NA,NA,2,2),bty = "n")    
            dev.off()