################################################################################
###### TSD25 - Example 3  #######################################
################################################################################
#load required libraries
library(R2WinBUGS);library(MASS);library(MCMCpack)
require(mcmcplots)

### read OC sensor data
data <- read.csv("OC-Sensor.csv",na="NA",header = TRUE)
## bring data in form needed for the model
ns <- length(data$ID) # No of studies
#table with number of thresholds per study
Tc <- data$T
Tc<-cbind(Tc,Tc)
#table with number of patients in the diseased and disease-free group
N <- as.matrix(data[, grepl("N", names(data))])
N<-cbind(N[,2],N[,1])
#table with threshold vales per study
C <- as.matrix(data[, grepl("C", names(data))])
tp <- as.data.frame(data[, grepl("tp", names(data))])
fp <- as.data.frame(data[, grepl("fp", names(data))])

#create array with counts as needed for the multiple thresholds model
x <- array(NA, c(ns, 2, max(Tc[,1])))
for(i in 1:ns){
  for(t in 1:max(Tc[,1])){
    x[i, 2, t] <- fp[i,t]
    x[i, 1, t] <- tp[i,t]
  }
}
#list with data neded for winbugs
dataList = list(ns = ns, Tc = Tc, N = N, C = C, x = x)
data.names<-c(names(dataList) )
### box-cox version #########################################################

#initial values for three chains
inits1<-list(  
  list(mean = c(3, 5,log(1), log(1)), 
       sd = c(1, 1, 1, 1),
       rho_mu = 0.8, 
       rho_mu_sigma = 0.2 , lambda=-2) ,
  list(mean = c(0, 6,log(2), log(2)), 
       sd = c(0.1, 0.1, 0.1, 0.1),
       rho_mu = 0.6, 
       rho_mu_sigma = 0.1 , lambda=-2.1) ,
  list(mean = c(2, 4,log(5), log(5)), 
       sd = c(0.1, 0.2, 0.1, 0.2),
       rho_mu = 0.1, 
       rho_mu_sigma = 0.1 , lambda=-2.5)
)  

#vector with parameters to monitor in winbugs
mymonitoredparamslist <- c( "mean","resdev","sd","lambda","mupred","logspred","s")#

# defining the directory of WinBUGS, this should be replaced with the full path to the directory WinBUGS is located
winbugs.dir <- "C:/My_winbugs_directory/WinBUGS14/"
##Run the model 
#model_boxcox.txt contains the WinBUGS code for the multiple thresholds model and should be saved in your current working directory
set.seed(213)
model2.sim <- bugs(data.names, inits1, model.file = "model_boxcox.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 400000, n.burnin=150000, n.thin=10,  bugs.directory = winbugs.dir, debug=FALSE)
#print summary od results
print(model2.sim,3)



########  log version   ########################################################

#initian values for three chains
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


#parameters to monitor in winbugs
mymonitoredparamslist <- c( "mean","resdev","sd","mupred","logspred")#


##Run the model 
#model_log.txt contains the WinBUGS code for the multiple thresholds model and should be saved in your current working directory
set.seed(213)
model1.sim <- bugs(data.names, inits1, model.file = "model_log.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 150000, n.burnin=50000, n.thin=5,  bugs.directory = winbugs.dir, debug=FALSE)

print(model1.sim,3)


################################################################################
##### code for threshold plots (Figure 6) ######################################
################################################################################

##log

nsims = length(model1.sim$sims.list$resdev)# total number of simulations
#save required posterior samples
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

#create vector with threshold values for which we wish to calculate TPF and FPF for
threshold<-na.omit(as.vector(C))
thres<-sort(threshold)

thres<-seq(min(thres),max(thres),length.out=800)
m<-length(thres) #No of total thresholds

#for each mcmc iteration and each threshold calculate TPF and FPF along with their predictive values
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

#calculate median TPF and FPF for each threshold along with their corresponding 95% CrI and PrI
summ_fpr <- summ_tpr <- matrix(NA, m, 3)
summ_fpr_pred <- summ_tpr_pred <- matrix(NA, m, 3)
for(t in 1:m){
  summ_fpr[t,] <- quantile(pr[,2,t], c(0.5, 0.025, 0.975))
  summ_tpr[t,] <- quantile(pr[,1,t], c(0.5, 0.025, 0.975))
  summ_fpr_pred[t,] <- quantile(pr_pred[,2,t], c(0.5, 0.025, 0.975))
  summ_tpr_pred[t,] <- quantile(pr_pred[,1,t], c(0.5, 0.025, 0.975))
}

#### summary sense spec for each thres #########################################
sumtprfpr<-round(data.frame(thres, summ_tpr, summ_fpr), 2)
write.table(sumtprfpr, "Summary_sens_spec_log.txt", sep=",")

sumtprfpr_pred<-round(data.frame(thres, summ_tpr_pred, summ_fpr_pred), 2)
write.table(sumtprfpr_pred, "Summary_sens_spec_pred_log.txt", sep=",")


#### reperate the same process for the Box-Coc version
nsims = length(model2.sim$sims.list$resdev)# total number of simulations
m_mu <- m_sigma <- matrix(NA, nsims, 2)
m_mu_pred <- m_sigma_pred <- matrix(NA, nsims, 2)
m_mu[,1] <- model2.sim$sims.list$mean[,1]
m_mu[,2] <- model2.sim$sims.list$mean[,2]
m_sigma[,1] <- model2.sim$sims.list$mean[,3]
m_sigma[,2] <- model2.sim$sims.list$mean[,4]
lambda<-model2.sim$sims.list$lambda

m_mu_pred[,1] <- model2.sim$sims.list$mupred[,1]
m_mu_pred[,2] <- model2.sim$sims.list$mupred[,2]
m_sigma_pred[,1] <- model2.sim$sims.list$logspred[,1]
m_sigma_pred[,2] <- model2.sim$sims.list$logspred[,2]

g <- matrix(NA, nsims, m)
logit_pr2 <- pr2 <- array(NA, c(nsims, 2, m )) 
logit_pr_pred2 <- pr_pred2 <- array(NA, c(nsims, 2, m )) 
for(k in 1:nsims){
  for(t in 1:m){
    if(lambda[k] == 1){g[k,t] <- log(thres[t])}
    if(lambda[k] != 1){g[k,t] <- ((thres[t])^lambda[k] - 1)/lambda[k] }
    for(j in 1:2){
      logit_pr2[k,j,t]<-(m_mu[k,j]-g[k,t])/exp(m_sigma[k,j])
      logit_pr_pred2[k,j,t]<-(m_mu_pred[k,j]-g[k,t])/exp(m_sigma_pred[k,j])
      
      # Undo the logit transformation:
      pr2[k,j,t] <- plogis(logit_pr2[k,j,t])
      pr_pred2[k,j,t] <- plogis(logit_pr_pred2[k,j,t])
    }
  }
}

summ_fpr2 <- summ_tpr2 <- matrix(NA, m, 3)
summ_fpr_pred2 <- summ_tpr_pred2 <- matrix(NA, m, 3)
for(t in 1:m){
  summ_fpr2[t,] <- quantile(pr2[,2,t], c(0.5, 0.025, 0.975))
  summ_tpr2[t,] <- quantile(pr2[,1,t], c(0.5, 0.025, 0.975))
  summ_fpr_pred2[t,] <- quantile(pr_pred2[,2,t], c(0.5, 0.025, 0.975))
  summ_tpr_pred2[t,] <- quantile(pr_pred2[,1,t], c(0.5, 0.025, 0.975))
}

#### summary sense spec for each thres #########################################
sumtprfpr2<-round(data.frame(thres, summ_tpr2, summ_fpr2), 2)
write.table(sumtprfpr2, "Summary_sens_spec_boxcox.txt", sep=",")

sumtprfpr_pred2<-round(data.frame(thres, summ_tpr_pred2, summ_fpr_pred2), 2)
write.table(sumtprfpr_pred2, "Summary_sens_spec_boxcox_pred.txt", sep=",")
#################################################################################

#calculate observed TPF and FPF of each study at each threshold
sens<-matrix(NA,ncol=10,nrow=ns)
fpr<-matrix(NA,ncol=10,nrow=ns)
for (i in 1:ns){
  sens[i,1:10]<-(as.matrix(tp[i,])/N[i,1])
  fpr[i,1:10]<-(as.matrix(fp[i,])/N[i,2])
}

#create the plot
#the plot will be saved in the working directory unless otherwise specified within the svg command
svg("thresplot_pred2.svg", width = 16, height = 7, pointsize = 12) # open file to output to
par(mfrow=c(1,2))
#create base of plot
plot(log(C[,1:10]),sens[,1:10],col="white",type="l",xlab=expression("Threshold ("*mu *"g/g)"),ylab="probability of a positive test result",main = "Log version",ylim=c(0,1),xlim=c(min(log(thres)),max(log(thres))),frame.plot =FALSE,xaxt="n",cex.lab=1.2, cex.axis=1,cex.main=1)
#points of x axis
aa<-c(5,10,20,50,100,200)
#add x axis
axis(side = 1,at=log(aa),labels=aa)

#add observed points with their size varying according to sample size
# add lines connecting points from the same study

mycol <- rgb(255, 0, 0, max = 255, alpha = 55, names = "myred")
mycol2 <- rgb(0, 0, 255, max = 255, alpha = 55, names = "myblue")



xno<-NULL
yno<-NULL
wno<-NULL
yno2<-NULL
wno2<-NULL

for(i in 1:ns){
  lines(log(C[i,1:10]),sens[i,1:10],col="darkgray",type="l")
  newx<-log(C[i,1:10])
  newy<-sens[i,1:10]
  newy2<-fpr[i,1:10]
  neww<-rep(N[i,1],10)
  neww2<-rep(N[i,2],10)
  xno<-c(xno,newx)
  yno<-c(yno,newy)
  wno<-c(wno,neww)
  yno2<-c(yno2,newy2)
  wno2<-c(wno2,neww2)
  lines(log(C[i,1:10]),fpr[i,1:10],col="darkgray",type="l")
}

wnos<-wno^(1/4)#rescale size
symbols(x=xno, y=yno, circles=wnos, inches=1/7,
        ann=F, bg=mycol, fg=NULL,add=TRUE)
wnos2<-wno^(1/4)#rescale size
symbols(x=xno, y=yno2, circles=wnos2, inches=1/7,
        ann=F, bg=mycol2, fg=NULL,add=TRUE)

#addd summary TPF and FPF lines
lines(log(thres), summ_fpr[,1], type = "l",col="blue",lwd=3)
lines(log(thres), summ_tpr[,1], type = "l",col="red",lwd=3)

### add 95% cris on plot
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
    
    ####box-cox
    plot(log(C[,1:10]),sens[,1:10],col="white",type="l",xlab=expression("Threshold ("*mu *"g/g)"),ylab="probability of a positive test result",main = "Box-Cox version",ylim=c(0,1),xlim=c(min(log(thres)),max(log(thres))),frame.plot =FALSE,xaxt="n",cex.lab=1.2, cex.axis=1,cex.main=1)
    aa<-c(5,10,20,50,100,200)
    axis(side = 1,at=log(aa),labels=aa)
    xno<-NULL
    yno<-NULL
    wno<-NULL
    yno2<-NULL
    wno2<-NULL
    for(i in 1:ns){
      lines(log(C[i,1:10]),sens[i,1:10],col="darkgray",type="l")
      newx<-log(C[i,1:10])
      newy<-sens[i,1:10]
      newy2<-fpr[i,1:10]
      neww<-rep(N[i,1],10)
      neww2<-rep(N[i,2],10)
      xno<-c(xno,newx)
      yno<-c(yno,newy)
      wno<-c(wno,neww)
      yno2<-c(yno2,newy2)
      wno2<-c(wno2,neww2)
      lines(log(C[i,1:10]),fpr[i,1:10],col="darkgray",type="l")
      
    }
    wnos<-wno^(1/4)#rescale size
    symbols(x=xno, y=yno, circles=wnos, inches=1/7,
            ann=F, bg=mycol, fg=NULL,add=TRUE)
    wnos2<-wno2^(1/4)#rescale size
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
        
        polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
                 c(summ_tpr_pred2[,2], rev(summ_tpr_pred2[,3])),
                 col=adjustcolor(coltp, alpha.f = 0.1),
                 border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
        
        polygon( c(log(thres[1:m]), rev(log(thres[1:m]))),
                 c(summ_fpr_pred2[,2], rev(summ_fpr_pred2[,3])),
                 col=adjustcolor(colfp, alpha.f = 0.05),
                 border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))   
        
        
        
        legend(x = "topright",          # Position
               legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)","Summary TPF", "Summary FPF"),  # Legend texts
               lty = c(NA,NA,1, 1),           # Line types
               col = c("red", "blue","red","blue"), # Line colors
               pch = c(20, 20, NA, NA),
               lwd = c(NA,NA,2,2),bty = "n")    
        
        dev.off() # turn off device   
        
################################################################################        
#### plot for checking the linearity assumption (Figure 7) #####################
        
  #calculate observed logit TPF and FPF for each study and threshold      
        logitsens<-matrix(NA,ncol=10,nrow=ns)
        logitfpr<-matrix(NA,ncol=10,nrow=ns)
        for (i in 1:ns){
          logitsens[i,1:10]<-log(as.matrix(sens[i,])/(1-as.matrix(sens[i,])))
          logitfpr[i,1:10]<-log(as.matrix(fpr[i,])/(1-as.matrix(fpr[i,])))
        }
  #create plot
  #the plot will be saved in the working directory unless otherwise specified within the svg command      
        svg("linear.svg", width = 8, height = 7, pointsize = 14) # open file to output to
        plot(log(C[,1:10]),logitsens[,1:10],col="white",type="l",ylim=c(-5,5),xlab="log(Threshold)",ylab="logit probability of a positive test result",xlim=c(1,5.4),xaxt="n",frame.plot =FALSE,cex.lab=1.2, cex.axis=1.1,cex.main=1) 
        #define x-axis points
        aa<-c(1,2,3,4,5)
        #add x-axis
        axis(side = 1,at=aa,labels=aa,cex.axis=1.1)
        #add logit TPF and FPF with lines connecting points from the same study
        for(i in 1:ns){
          lines(log(C[i,1:10]),logitsens[i,1:10],col="pink",type="l")
          lines(log(C[i,1:10]),logitsens[i,1:10],col=2,type="p",pch=20)
          lines(log(C[i,1:10]),logitfpr[i,1:10],col="lightblue",type="l")
          lines(log(C[i,1:10]),logitfpr[i,1:10],col="blue",type="p",pch=20)
        }
        
        legend(x = "topright",          # Position
               legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)"),  # Legend texts
               lty = c(NA,NA),           # Line types
               col = c("red", "blue"), # Line colors
               pch = c(20, 20),
               lwd = c(NA,NA),bty = "n")  
        dev.off()
        
     
   
        
        
        
 