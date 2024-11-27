################################################################################
#### TSD 25 -Example 3 #########################################################
################################################################################

#load required libraries
library(R2WinBUGS);library(MASS);library(MCMCpack)
require(mcmcplots)
library(coda)

### function for calculating pD outside winbugs

pD_calc<-function(bugsres,tp,fp,Tc,ns,N){ #bugsres is the simulation matrix form winbugs output

totaldatapoints <- sum(Tc[,1])+sum(Tc[,2])

# Get mean of fitted values ('rhat'):
tphat.mean <- apply(bugsres[,grepl("xhatt", names(bugsres))],2,mean)
fphat.mean <- apply(bugsres[,grepl("xhatf", names(bugsres))],2,mean)

# fitted logit TPR and logit FPR:
dt.mean <- apply(bugsres[,grepl("dt", names(bugsres))],2,mean)
df.mean <- apply(bugsres[,grepl("df", names(bugsres))],2,mean)

# Put observed values in same format:
tp.obs <- c(tp[1,1:Tc[1,1]])
fp.obs <- c(fp[1,1:Tc[1,2]])

NT <- max(Tc[,1],Tc[,2]) 
n0 <- n1 <- matrix(NA, ns, NT)

n0[,1] <- N[,2]
n1[,1] <- N[,1]
for(s in 1:ns){
  if(Tc[s,1] > 1){
    for(t in 2:Tc[s,1]){
      n0[s,t] <- fp[s,(t-1)]
      n1[s,t] <- tp[s,(t-1)]
    }
  }
}

n1.obs <- c(n1[1,1:Tc[1,1]])
n0.obs <- c(n0[1,1:Tc[1,1]])
C.obs <- c(C[1,1:Tc[1,1]])
for(s in 2:ns){
  tp.obs <- c(tp.obs, tp[s, 1:Tc[s,1]])
  fp.obs <- c(fp.obs, fp[s, 1:Tc[s,1]])
  n1.obs <- c(n1.obs, n1[s, 1:Tc[s,1]])
  n0.obs <- c(n0.obs, n0[s, 1:Tc[s,1]])
  C.obs <- c(C.obs, C[s, 1:Tc[s,1]])
}
# Get mean residual deviance (dev):
mean.dev.fp <- apply(bugsres[,grepl("devf", names(bugsres))],2,mean)
mean.dev.tp <- apply(bugsres[,grepl("devt", names(bugsres))],2,mean)

# Calculate deviance at posterior mean for each data point:
# FP:
fp.dev <- tp.dev <- rep(NA, sum(Tc[,1])) 
for(i in 1:sum(Tc[,1])){
  
  fp.dev[i] <-  2*( 
    as.numeric(fp.obs[i])*(log(as.numeric(fp.obs[i])) - log(fphat.mean[i]))
    + (n0.obs[i] - as.numeric(fp.obs[i]))*(log(n0.obs[i] - as.numeric(fp.obs[i])) - log(n0.obs[i] - fphat.mean[i]))
  ) 
  
  if(as.numeric(fp.obs[i]) ==0){
    fp.dev[i] <-  2*( 
      (n0.obs[i] - as.numeric(fp.obs[i]))*(log(n0.obs[i] - as.numeric(fp.obs[i])) - log(n0.obs[i] - fphat.mean[i]))
    ) 
  }
  
  if(as.numeric(fp.obs[i]) == n0.obs[i]){
    fp.dev[i] <-  2*( 
      as.numeric(fp.obs[i])*(log(as.numeric(fp.obs[i])) - log(fphat.mean[i]))
    )
  } 
  
  tp.dev[i] <-  2*( 
    as.numeric(tp.obs[i])*(log(as.numeric(tp.obs[i]) )- log(tphat.mean[i]))
    + (n1.obs[i] - as.numeric(tp.obs[i]))*(log(n1.obs[i] - as.numeric(tp.obs[i])) - log(n1.obs[i] - tphat.mean[i]))
  ) 
  
  if(as.numeric(tp.obs[i]) == n1.obs[i]){
    tp.dev[i] <-  2*( 
      as.numeric(tp.obs[i])*(log(as.numeric(tp.obs[i])) - log(tphat.mean[i]))
    ) 
  }
  
}

# Leverage (contribution to pD) for each data point:
fp.lev <- mean.dev.fp - fp.dev
tp.lev <- mean.dev.tp - tp.dev

# pD = sum of these:
pd <- sum(fp.lev[]) + sum(tp.lev[])
return(pd)
}

### Read OC-Sensor data
data <- read.csv("OC-Sensor.csv",na="NA",header = TRUE)
#bring data in form needed for the model
ns <- length(data$ID) # No of studies
#table with number of patients in the diseased and disease-free group
N <- as.matrix(data[, grepl("N", names(data))])
N <- cbind(N[,2],N[,1])
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
  for(t in 1:max(Tc[,1],Tc[,2])){
    x[i, 2, t] <- fp[i,t]
    x[i, 1, t] <- tp[i,t]
  }
}

#list with data neded for winbugs
dataList = list(ns = ns, Tc = Tc, N = N, C = C, x = x)
data.names<-c(names(dataList) )

### plot for checking the linearity assumption (Figure 6) ######################

maxthr<-max(Tc[,1],Tc[,2]) #no of maximum thresholds in a study
#calculate observed TPF,FPF and their logit transforms
sens<-matrix(NA,ncol=maxthr,nrow=ns)
fpr<-matrix(NA,ncol=maxthr,nrow=ns)
logitsens<-matrix(NA,ncol=maxthr,nrow=ns)
logitfpr<-matrix(NA,ncol=maxthr,nrow=ns)
for (i in 1:ns){
  sens[i,1:maxthr]<-(as.matrix(tp[i,])/N[i,1])
  fpr[i,1:maxthr]<-(as.matrix(fp[i,])/N[i,2])
  logitsens[i,1:maxthr]<-log(as.matrix(sens[i,])/(1-as.matrix(sens[i,])))
  logitfpr[i,1:maxthr]<-log(as.matrix(fpr[i,])/(1-as.matrix(fpr[i,])))
}
#create and save plot
svg("linear.svg", width = 8, height = 7, pointsize = 14) # open file to output to
plot(log(C[,1:maxthr]),logitsens[,1:maxthr],col="white",type="l",ylim=c(-5,5),xlab="log(Threshold)",ylab="logit probability of a positive test result",xlim=c(1,5.4),xaxt="n",frame.plot =FALSE,cex.lab=1.2, cex.axis=1.1,cex.main=1) 
aa<-c(1,2,3,4,5)

axis(side = 1,at=aa,labels=aa,cex.axis=1.1)

for(i in 1:ns){
  lines(log(C[i,1:maxthr]),logitsens[i,1:maxthr],col="pink",type="l")
  lines(log(C[i,1:maxthr]),logitsens[i,1:maxthr],col=2,type="p",pch=20)
  lines(log(C[i,1:maxthr]),logitfpr[i,1:maxthr],col="lightblue",type="l")
  lines(log(C[i,1:maxthr]),logitfpr[i,1:maxthr],col="blue",type="p",pch=20)
}

legend(x = "topright",          # Position
       legend = c("Observed diseased (TPF)", "Observed disease-free (FPF)"),  # Legend texts
       lty = c(NA,NA),           # Line types
       col = c("red", "blue"), # Line colors
       pch = c(20, 20),
       lwd = c(NA,NA),bty = "n")  
dev.off()

################################################################################
### Box-Cox + restricted cov structure #########################################
################################################################################

# initial values
inits1<-list(  
  list(mean = c(5, 3,log(1), log(1)), 
       sd = c(1, 1, 1, 1),
       rho_mu = 0.8, 
       rho_mu_sigma = 0.2,lambda=-0.8 ) ,
  list(mean = c(6, 0,log(2), log(2)), 
       sd = c(0.1, 0.1, 0.1, 0.1),
       rho_mu = 0.6, 
       rho_mu_sigma = 0.1,lambda=-0.5 ) ,
  list(mean = c(4, 2,log(5), log(5)), 
       sd = c(0.2, 0.1, 0.2, 0.1),
       rho_mu = 0.1, 
       rho_mu_sigma = 0.1,lambda=-0.3 )
)  

#vector with parameters to monitor in winbugs
mymonitoredparamslist <- c( "mean","resdev","sd","mupred","s","xhatt","xhatf","pr","dt","df","devt","devf","lambda")#

# defining the directory of WinBUGS, this should be replaced with the full path to the directory WinBUGS is located
winbugs.dir <- "C:/Users/cb22323/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14/"
##Run the model 
#model_boxcox.txt contains the WinBUGS code for the multiple thresholds model and should be saved in your current working directory
set.seed(213)
model1.sim <- bugs(data.names, inits1, model.file = "model_boxcox_restr.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 800000, n.burnin=480000, n.thin=15,  bugs.directory = winbugs.dir, debug=FALSE)

#print summary of results
print(model1.sim,3)

#save matrix with posterior samples
bugsres<-as.data.frame(model1.sim$sims.matrix)

# calculate number of effective parameters
pD_calc(bugsres=bugsres,tp=tp,fp=fp,Tc=Tc,ns=ns,N=N)

################################################################################
#### Box-cox + independence cov structure ######################################
################################################################################

#initial values
inits1<-list(  
  list(mean = c(5, 3,log(1), log(1)), 
       sd = c(1, 1, 1, 1),lambda=-0.2) ,
  list(mean = c(6, 0,log(0.5), log(0.5)), 
       sd = c(0.1, 0.1, 0.1, 0.1),lambda=-0.5) ,
  list(mean = c(4, 2,log(1.2), log(1.2)), 
       sd = c(0.2, 0.1, 0.2, 0.1),lambda=-0.3)
)  

#vector with parameters to monitor in winbugs
mymonitoredparamslist <- c( "mean","resdev","sd","mupred","s","xhatt","xhatf","pr","dt","df","devt","devf","lambda")#

##Run the model 
#model_boxcox.txt contains the WinBUGS code for the multiple thresholds model and should be saved in your current working directory
set.seed(213)

model2.sim <- bugs(data.names, inits1, model.file = "model_boxcox_ind.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 600000, n.burnin=300000, n.thin=10,  bugs.directory = winbugs.dir, debug=FALSE)

#print summary of results
print(model2.sim,3)

#save matrix with posterior samples
bugsres<-as.data.frame(model2.sim$sims.matrix)

# calculate number of effective parameters
pD_calc(bugsres=bugsres,tp=tp,fp=fp,Tc=Tc,ns=ns,N=N)

################################################################################
###### Box-Cox + unrestricted covariance structure #############################
################################################################################

# initial values
inits1<-list(  
  list(mean = c(5, 3,log(1), log(1)), 
       sd = c(1, 1, 1, 1),
       rho_mu = 0.8, 
       rho_mu0_sigma0 = 0.2,
       rho_mu1_sigma1=0.1,
       rho_mu0_sigma1=0.1,
       rho_mu1_sigma0=0.1,
       rho_sigma=0.3,lambda=-0.2) ,
  list(mean = c(6, 0.1,log(0.5), log(0.5)), 
       sd = c(0.1, 0.1, 0.1, 0.1),
       rho_mu = 0.6, 
       rho_mu0_sigma0 = 0.24,
       rho_mu1_sigma1=0.2,
       rho_mu0_sigma1=0.1,
       rho_mu1_sigma0=0.2,
       rho_sigma=0.6 ,lambda=-0.5) ,
  list(mean = c(4, 0.2,log(1.2), log(1.2)), 
       sd = c(0.2, 0.1, 0.2, 0.1),
       rho_mu = 0.1, 
       rho_mu0_sigma0 = 0.1,
       rho_mu1_sigma1=0.15,
       rho_mu0_sigma1=0.15,
       rho_mu1_sigma0=0.15,
       rho_sigma=0.4 ,lambda=-0.3)
)  

mymonitoredparamslist <- c( "mean","resdev","sd","rho_mu","rho_sigma","mupred","pr","xhatt","xhatf","dt","df","devt","devf","lambda")#"resdevm")#

##Run the model 
set.seed(213)
model3.sim <- bugs(data.names, inits1, model.file = "model_boxcox_full.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 1500000, n.burnin=800000, n.thin=15,  bugs.directory = winbugs.dir, debug=TRUE)

print(model3.sim,3)

#save matrix with posterior samples
bugsres<-as.data.frame(model3.sim$sims.matrix)

# calculate number of effective parameters
pD_calc(bugsres=bugsres,tp=tp,fp=fp,Tc=Tc,ns=ns,N=N)

################################################################################
########  code for threshold plot (Figure 7) ###################################
################################################################################

nsims = length(model2.sim$sims.list$resdev)# total number of simulations
# save required posterior samples
m_mu <- m_sigma <- matrix(NA, nsims, 2)
m_mu_pred <- m_sigma_pred <- matrix(NA, nsims, 2)
m_mu[,1] <- model2.sim$sims.list$mean[,1]
m_mu[,2] <- model2.sim$sims.list$mean[,2]
m_sigma[,1] <- model2.sim$sims.list$mean[,3]
m_sigma[,2] <- model2.sim$sims.list$mean[,4]
lambda<-model2.sim$sims.list$lambda #comment this line out if log-link used

m_mu_pred[,1] <- model2.sim$sims.list$mupred[,1]
m_mu_pred[,2] <- model2.sim$sims.list$mupred[,2]
m_sigma_pred[,1] <- model2.sim$sims.list$mupred[,3]
m_sigma_pred[,2] <- model2.sim$sims.list$mupred[,4]

#create vector with threshold values for which we wish to calculate TPF and FPF for
threshold<-na.omit(as.vector(C))
thres<-sort(threshold)

thres<-seq(min(thres),max(thres),length.out=800)
m<-length(thres)

#for each mcmc iteration and each threshold calculate TPF and FPF along with their predictive values
g <- matrix(NA, nsims, m) # comment this out if log-link used
logit_pr <- pr <- array(NA, c(nsims, 2, m )) 
logit_pr_pred <- pr_pred <- array(NA, c(nsims, 2, m )) 
for(k in 1:nsims){
  for(t in 1:m){
    if(lambda[k] == 1){g[k,t] <- log(thres[t])} # comment this line out if log-link used
    if(lambda[k] != 1){g[k,t] <- ((thres[t])^lambda[k] - 1)/lambda[k] } # comment this line out if log-link used
    for(j in 1:2){
      logit_pr[k,j,t]<-(m_mu[k,j]-g[k,t])/exp(m_sigma[k,j])
      logit_pr_pred[k,j,t]<-(m_mu_pred[k,j]-g[k,t])/exp(m_sigma_pred[k,j])
      # if log-link used comment out the 2 previous lines and use the following instead
      #logit_pr2[k,j,t]<-(m_mu[k,j]-log(thres[t]))/exp(m_sigma[k,j])
      #logit_pr_pred2[k,j,t]<-(m_mu_pred[k,j]-log(thres[t]))/exp(m_sigma_pred[k,j])
      
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

#### summary sense spec for each thres 
sumtprfpr<-round(data.frame(thres, summ_tpr, summ_fpr), 2)
write.table(sumtprfpr, "Summary_sens_spec_boxcox_ind.txt", sep=",")

sumtprfpr_pred<-round(data.frame(thres, summ_tpr_pred, summ_fpr_pred), 2)
write.table(sumtprfpr_pred, "Summary_sens_spec_boxcox_pred_ind.txt", sep=",")

#create the plot
#the plot will be saved in the working directory unless otherwise specified within the svg command
svg("thresplot.svg", width = 12, height = 11, pointsize = 14) # open file to output to
#create base of plot
plot(log(C[,1:maxthr]),sens[,1:maxthr],col="white",type="l",xlab=expression("Threshold ("*mu *"g/g)"),ylab="probability of a positive test result",main = "Box-Cox version",ylim=c(0,1),xlim=c(min(log(thres)),max(log(thres))),frame.plot =FALSE,xaxt="n",cex.lab=1.2, cex.axis=1.2,cex.main=1.2)
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
wnos<-wno^(1/2)#sqrt(wno)
symbols(x=xno, y=yno, circles=wnos, inches=1/7,
        ann=F, bg=mycol, fg=NULL,add=TRUE)
wnos2<-wno2^(1/2)#sqrt(wno2)
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
           lwd = c(NA,NA,2,2),bty = "n",cex=1.1)    
    
    dev.off() # turn off device   

