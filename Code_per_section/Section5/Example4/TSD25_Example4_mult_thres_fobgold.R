################################################################################
###### TSD25 - Example 4  #######################################
################################################################################
#load required libraries
library(R2WinBUGS);library(MASS);library(MCMCpack)
require(mcmcplots)
################################################################################
#read FobGold data
data <- read.csv("FOBGold.csv",na="NA",header = TRUE)
## bring data in form needed for the model
ns <- length(data$ID) # No of studies

#table with number of patients in the diseased and disease-free group
N <- as.matrix(data[, grepl("N", names(data))])
N<-cbind(N[,2],N[,1])
#table with threshold vales per study
C <- as.matrix(data[, grepl("C", names(data))])
tp <- as.data.frame(data[, grepl("tp", names(data))])
fp <- as.data.frame(data[, grepl("fp", names(data))])

# table with number of thresholds per study
Tc <- data$T
# remove any consequite zero FP counts
T1 <- Tc
for(i in 1:ns){
  if(Tc[i] > 1){
    for(j in 1:(Tc[i]-1)){
      if(fp[i,j] == 0 | is.na(fp[i,j])){
        fp[i,j+1] <- NA
        T1[i] <- T1[i] - 1
      }
    }
  }
}
Tc<-cbind(Tc,T1)

#create array with counts as needed for the multiple thresholds model
x <- array(NA, c(ns, 2, max(Tc[,1],Tc[,2])))
for(i in 1:ns){
  for(t in 1:max(Tc[,1])){
    x[i, 2, t] <- fp[i,t]
    x[i, 1, t] <- tp[i,t]
  }
}


#list with data neded for winbugs
dataList = list(ns = ns, Tc = Tc, N = N, C = C, x = x)
data.names<-c(names(dataList) )

##### plot to check the linearity assumption ###################################
################################################################################

maxthr<-max(Tc[,1],Tc[,2]) #no of maximum thresholds in a study
### calculate observed sens and fpf of each study
sens<-matrix(NA,ncol=maxthr,nrow=ns)
fpr<-matrix(NA,ncol=maxthr,nrow=ns)
for (i in 1:ns){
  sens[i,1:maxthr]<-(as.matrix(tp[i,])/N[i,1])
  fpr[i,1:maxthr]<-(as.matrix(fp[i,])/N[i,2])
}
# calculate logit version of the above
logitsens<-matrix(NA,ncol=maxthr,nrow=ns)
logitfpr<-matrix(NA,ncol=maxthr,nrow=ns)
for (i in 1:ns){
  logitsens[i,1:maxthr]<-log(as.matrix(sens[i,])/(1-as.matrix(sens[i,])))
  logitfpr[i,1:maxthr]<-log(as.matrix(fpr[i,])/(1-as.matrix(fpr[i,])))
}

#this plot will be saved at the current working directory
svg("linear.svg", width = 8, height = 7, pointsize = 14) # open file to output to
# create base plot
plot(log(C[,1:4]),logitsens[,1:4],col="white",type="l",ylim=c(-5,5),xlab="log(Threshold)",ylab="logit probability of a positive test result",xlim=c(0,5.4),xaxt="n",frame.plot =FALSE,cex.lab=1.2, cex.axis=1.1,cex.main=1) 
# x-axis points
aa<-c(0,1,2,3,4,5)
# add x-axis
axis(side = 1,at=aa,labels=aa,cex.axis=1.1)
#add logit sens, fpf points connected with a line if they come from the same study
for(i in 1:ns){
  lines(log(C[i,1:4]),logitsens[i,1:4],col="pink",type="l")
  lines(log(C[i,1:4]),logitsens[i,1:4],col=2,type="p",pch=20)
  lines(log(C[i,1:4]),logitfpr[i,1:4],col="lightblue",type="l")
  lines(log(C[i,1:4]),logitfpr[i,1:4],col="blue",type="p",pch=20)
}
# add a legend
legend(x = "topright",          # Position
       legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)"),  # Legend texts
       lty = c(NA,NA),           # Line types
       col = c("red", "blue"), # Line colors
       pch = c(20, 20),
       lwd = c(NA,NA),bty = "n")  
dev.off() # to close graphical device and comple the plot

################################################################################
###### fit vague prior version of the model with indipendent r.e.s #############
###############################################################################
#initial values
inits1<-list(  
  list(mean = c(5, -3,log(1), log(1)), 
       sd = c(1, 1, 1, 1) ) ,
  list(mean = c(6, 0,log(2), log(2)), 
       sd = c(0.1, 0.1, 0.1, 0.1) ) ,
  list(mean = c(4, -2,log(5), log(5)), 
       sd = c(0.2, 0.1, 0.2, 0.1))
)  

#vector with parameters to monitor in winbugs
mymonitoredparamslist <- c( "mean","resdev","sd","mupred","logspred")#

# defining the directory of WinBUGS, this should be replaced with the full path to the directory WinBUGS is located
winbugs.dir <- "C:/My_winbugs_directory/WinBUGS14/"
##Run the model 
#model_log_ind.txt contains the WinBUGS code for the multiple thresholds model and should be saved in your current working directory
set.seed(213)
model1.sim <- bugs(data.names, inits1, model.file = "model_log_ind.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 300000, n.burnin=100000, n.thin=5,  bugs.directory = winbugs.dir, debug=FALSE)

################################################################################
##### run the fixed-effects version ############################################
################################################################################

#initial values
inits1<-list(  
  list(mean = c(5, 3,log(1), log(1)) ) ,
  list(mean = c(6, 0,log(2), log(2)) ) ,
  list(mean = c(4, 2,log(5), log(5)))
)  

#parameters to monitor
mymonitoredparamslist <- c( "mean","resdev")#

##Run the model 
set.seed(213)
model2.sim <- bugs(data.names, inits1, model.file = "model_log_fe.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 150000, n.burnin=50000, n.thin=5,  bugs.directory = winbugs.dir, debug=FALSE)

################################################################################
### version with HN(0,1) priors ################################################
################################################################################

#initial values
inits1<-list(  
  list(mean = c(5, -3,log(1), log(1)), 
       sd = c(1, 1, 1, 1) ) ,
  list(mean = c(6, 0,log(2), log(2)), 
       sd = c(0.1, 0.1, 0.1, 0.1) ) ,
  list(mean = c(4, -2,log(5), log(5)), 
       sd = c(0.2, 0.1, 0.2, 0.1))
)  

#vector with parameters to monitor in winbugs
mymonitoredparamslist <- c( "mean","resdev","sd","mupred","logspred")#

#run  the model in winbugs
set.seed(213)
model3.sim <- bugs(data.names, inits1, model.file = "model_log_ind_hn.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 300000, n.burnin=100000, n.thin=5,  bugs.directory = winbugs.dir, debug=FALSE)

################################################################################
#### code for threshold plot (Figure 11) #######################################
################################################################################
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

# sequence of thresholds
threshold<-na.omit(as.vector(C))
thres<-sort(threshold)

thres<-seq(min(thres),max(thres),length.out=800)
m<-length(thres)

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
  summ_fpr3[t,] <- quantile(pr3[,2,t], c(0.5, 0.025, 0.975))
  summ_tpr3[t,] <- quantile(pr3[,1,t], c(0.5, 0.025, 0.975))
  summ_fpr_pred3[t,] <- quantile(pr_pred3[,2,t], c(0.5, 0.025, 0.975))
  summ_tpr_pred3[t,] <- quantile(pr_pred3[,1,t], c(0.5, 0.025, 0.975))
}

#### summary sense spec for each thres #########################################
sumtprfpr3<-round(data.frame(thres, summ_tpr3, summ_fpr3), 2)
write.table(sumtprfpr3, "Summary_sens_spec_log_infhn.txt", sep=",")

sumtprfpr_pred3<-round(data.frame(thres, summ_tpr_pred3, summ_fpr_pred3), 2)
write.table(sumtprfpr_pred3, "Summary_sens_spec_pred_loghn.txt", sep=",")


#### plot #######################################
svg("thresplot_fobhn.svg", width = 9, height = 7, pointsize = 12) # open file to output to
par(mfrow=c(1,1),mai = c(1.4, 1, 0.3, 0.2))
#plot(log(C[1,1:48]),sens[1,1:48],col=2,type="l",xlab="log-threshold",ylab="probability of a positive test result",ylim=c(0,1),xlim = c(0,10),frame.plot =FALSE )
plot(log(C[,1:maxthr]),sens[,1:maxthr],col="white",type="l",xlab=expression("Threshold ("*mu *"g Hb/g)"),ylab="probability of a positive test result",main = "Full model with HN(0,1) priors",ylim=c(0,1),xlim=c(min(log(thres)),max(log(thres))),frame.plot =FALSE,xaxt="n",cex.lab=1.2, cex.axis=1,cex.main=1)
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
  lines(log(C[i,1:maxthr]),sens[i,1:maxthr],col="darkgray",type="l")
  newx<-log(C[i,1:maxthr])
  newy<-sens[i,1:maxthr]
  newy2<-fpr[i,1:maxthr]
  neww<-rep(N[i,2],maxthr)
  neww2<-rep(N[i,1],maxthr)
  xno<-c(xno,newx)
  yno<-c(yno,newy)
  wno<-c(wno,neww)
  yno2<-c(yno2,newy2)
  wno2<-c(wno2,neww2)
  #lines(log(C[i,1:10]),sens[i,1:10],col=2,type="p",pch=20)
  lines(log(C[i,1:maxthr]),fpr[i,1:maxthr],col="darkgray",type="l")
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
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    
    legend(x = "bottom",          # Position
           legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)","Summary TPF", "Summary FPF"),  # Legend texts
           lty = c(NA,NA,1, 1),           # Line types
           col = c("red", "blue","red","blue"), # Line colors
           pch = c(20, 20, NA, NA),
           lwd = c(NA,NA,2,2),bty = "n",ncol=2)       
    dev.off() 
### Program ends    