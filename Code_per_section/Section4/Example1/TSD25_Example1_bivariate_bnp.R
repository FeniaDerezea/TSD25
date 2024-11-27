################################################################################
#### TSD25 - Example 1 ###########################################################
################################################################################
# load required libraries
library(R2WinBUGS)
require(mcmcplots)
library(ellipse)
#read data
data<-read.table("BNP_data.txt", header = TRUE, sep = ",")
# number of diseased
data$nA<-data$TP+data$FN
# number of healthy
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
parameter.names <- c( 'theta','rho','sd','sumtpf','sumfpf','spec','Theta','delta','beta','predtpr','predfpf','Lambda','vartheta','varalpha','sdtheta','sdalpha')


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

#run model through WinBUGS and safe output
#model1.txt contains the WinBUGS code for the bivariate model and should be saved in your current working directory
model1.sim <- bugs( data.names, inits1, model.file = "model_bivariate.txt", parameters = parameter.names,
                    n.chains = n.chain, n.iter = n.iters, n.burnin=n.burn, n.thin=n.th,  bugs.directory = winbugs.dir, debug=FALSE)

# print model summary rounded to 3 decimals
print(model1.sim,3)
#trace plot for saved parameters
mcmcplot(model1.sim)

################################################################################
### code for coupled forest plots ##############################################
################################################################################
#This code shows how to obtain a forest plot using the DTAplots library
#The code for producing the version of the forest plot shown in the main document is provided in the seperate R file ForestDTA
library(DTAplots) #call library DTAplots

datf<-data[,3:6] #dataframe with only the 2x2 counts
studies<-paste(data$Author,data$Year) #creates a vector with study IDs

#import pooled results from R2WinBUGS output
pooledresults <- matrix(NA, 1, 6)

sensres<-quantile(model1.sim$sims.list$sumtpf, c(0.5,0.025,0.975))
specres<-quantile(model1.sim$sims.list$spec, c(0.5,0.025,0.975))

pooledresults[1,] <- c( sensres[1] ,  sensres[2]  , sensres[3],
                        
                        specres[1] ,  specres[2]  , specres[3])

#an svg copy of the plot will be saved at the working directory
svg("forest.svg", width = 14, height = 12, pointsize = 14) 
Forest(datf, study =  studies,
       
       se.axis = c(0, 1),
       
       sp.axis = c(0, 1),
       
       summary = pooledresults,
       
       digits = 2,
       
       summary_label = "Summary estimates (Bivariate)"
       
) 

dev.off()



#### ellipses plot using DTAplots ##############################################
# This code shows how to obtain ellipses plots using the DTAplots library
# the code for producing the version in the TSD document is provided below
mu<-cbind(model1.sim$sims.list$theta[,1],-model1.sim$sims.list$theta[,2])
tau<-cbind(model1.sim$sims.list$sd[,1],model1.sim$sims.list$sd[,2])
rho<-matrix(-model1.sim$sims.list$rho,ncol=1)
Summary_Se<-matrix(model1.sim$sims.list$sumtpf,ncol=1)
Summary_Sp<-matrix(model1.sim$sims.list$spec,ncol=1)
se<-model1.sim$sims.list$delta[,,1]
sp<--model1.sim$sims.list$delta[,,2]
Xsroc<-cbind(mu,tau,rho,Summary_Se,Summary_Sp,se,sp)

colse<-colsp<-rep(NA,ns)
for (i in 1:ns){
  colse[i]<-paste("se","[",i,"]",sep="")
  colsp[i]<-paste("sp","[",i,"]",sep="")
}

colnames(Xsroc)<-c("mu[1]","mu[2]","tau[1]","tau[2]","rho","Summary_Se","Summary_Sp",colse,colsp)

svg("sroc_bnp.svg", width = 8, height = 8, pointsize = 12) 
SROC_rjags(X=Xsroc,n=ns,ref_std = TRUE,cred_region = TRUE,dataset=data[,3:6], predict_region = TRUE)
dev.off()

#### ellipses TSD version ###############################################
#calculate sens and spec for each study

senss<-data$TP/data$nA
fpfs<-data$FP/data$nB


#calculate the posterior correlation in the MCMC chains between logit sens and fpf
mcmc1<-model1.sim$sims.list$theta[,1]
mcmc2<-model1.sim$sims.list$theta[,2]
r<-cor(mcmc1,mcmc2)

#save summary sensitivity and FPF
sens<-median(model1.sim$sims.list$sumtpf)
fpf<-median(model1.sim$sims.list$sumfpf)

#collect reqired values for the plot in the same object
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
#the plot will be saved in the working directory unless otherwise specified within the svg command
svg("roc_logit2.svg", width = 14, height = 7, pointsize = 14) 
par(mfrow=c(1,2))
cex <- 1.2
plot(fpfs, senss,col="white",  pch = 20, xlim = c(0,1), ylim = c(0,1),
     cex = cex, cex.lab = cex, xlab = "1-Specificity (FPF)", ylab = "Sensitivity",
     axes = T,frame.plot =FALSE)
#add observed data with vrrying size according to sample size
mycol <- rgb(0, 0, 0, max = 255, alpha = 55, names = "myblack")
wno<-data$nA+data$nB
wnos<-wno^(1/2)#adjust scale of points
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
mcmc1_pred<-model1.sim$sims.list$predtpr
mcmc2_pred<-model1.sim$sims.list$predfpf
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

legend(x = "bottomright", cex = 1.2,          # Position
       legend = c("Observed (TPF, FPF)","Summary (TPF,FPF)", "95% credible region","95% prediction region"),  # Legend texts
       lty = c(NA,NA,2,2),           # Line types
       col = c(mycol, "red","red","blue"), # Line colors
       pch = c(20,17, NA,NA),
       lwd = c(NA,NA,2,2),bty = "n")   

## same plot logit scale 
#calculate logit-sens and logit-spec for each study

logitsenss<-log(senss/(1-senss))
logitfpfs<-log(fpfs/(1-fpfs))

plot(logitfpfs, logitsenss,col="white",  pch = 20,xlim=c(-5,5),ylim=c(-5,5),
     cex = cex, cex.lab = cex, xlab = "Logit(FPF)", ylab = "Logit(Sensitivity)",
     axes = T,frame.plot =FALSE)


symbols(x=logitfpfs, y=logitsenss, circles=wnos, inches=1/7,
        ann=F, bg=mycol, fg=NULL,add=TRUE)

points(results_biv$mean1[1], results_biv$mean0[1], col ="red", pch = 17,cex=cex)
lines(logitellipse1[,2],logitellipse1[,1], lty = 2, lwd = 2, col = "red") 
lines(logitellipse2[,2],logitellipse2[,1], lty = 2, lwd = 2, col = "blue") 
lines(-5:5,-5:5,lty=3,type="l")
legend(x = "bottomright", cex = 1.2,          # Position
       legend = c("Observed (logit-TPF, logit-FPF)","Summary (m1,m2)", "95% credible ellipse","95% prediction ellipse"),  # Legend texts
       lty = c(NA,NA,2,2),           # Line types
       col = c(mycol, "red","red","blue"), # Line colors
       pch = c(20,17, NA,NA),
       lwd = c(NA,NA,2,2),bty = "n")  
dev.off()
###############################################################################
#### ppv npv plot ###########################################################

nsims<-length(model1.sim$sims.list$spec) #mcmc chain length
#posterior chains of TPF and spec
sumtpr<-model1.sim$sims.list$sumtpf
spec<-model1.sim$sims.list$spec
#vector with possible prevalence values
prev <- seq(0,1,by=0.01)  
#calculate ppv and npv for each mcmc iteration and each prevalence
ppv<-npv<-matrix(NA,ncol=length(prev),nrow = nsims)
for (i in 1:nsims){
  for (j in 1:length(prev)){
    ppv[i,j] <- (prev[j]*sumtpr[i])/((prev[j]*sumtpr[i])+((1-prev[j])*(1-spec[i])))                         # ppv
    npv[i,j] <- ((1-prev[j])*spec[i])/(((1-prev[j])*spec[i])+(prev[j]*(1-sumtpr[i])))  
  }
}
#median  ppv and npv values and 95% CrIs
sumppv<-sumnpv<-matrix(NA,ncol=3,nrow=length(prev))
for(j in 1:length(prev)){
  sumppv[j,]<-quantile(ppv[,j],c(0.5,0.025,0.975))
  sumnpv[j,]<-1 - quantile(npv[,j],c(0.5,0.025,0.975))
}

#create the plot
svg("ppv.svg", width = 5, height = 5, pointsize = 10) 
plot(prev,sumppv[,1],type = "l",col="red",lwd=2,xlab = "Pre-test probability",ylab="Post-test probability",frame.plot =FALSE,cex.lab=1.2, cex.axis=1)
lines(prev,sumnpv[,1],type = "l",col="blue",lwd=2)

m<-length(prev)
coltp<-"red"
  colfp<-"blue"  
    polygon( c(prev[1:m], rev(prev[1:m])),
             c(sumppv[,2], rev(sumppv[,3])),
             col=adjustcolor(coltp, alpha.f = 0.3),
             border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
    
    polygon( c(prev[1:m], rev(prev[1:m])),
             c(sumnpv[,2], rev(sumnpv[,3])),
             col=adjustcolor(colfp, alpha.f = 0.2),
             border = NA, xlab = "", ylab = "", main = "",ylim = c(0,1))
    
    legend(0,1,          # Position
           legend = c("Prevalence or pre-test probability", "Positive test", "Negative test"),  # Legend texts
           lty = c(3, 1, 1),           # Line types
           col = c("black", "red", "blue"), # Line colors
           pch = c(NA, NA, NA),
           lwd = c(3, 2,2),bty = "n")   
    
    lines(c(0,1), c(0,1), lty = 3, lwd = 3)
    
    
    dev.off()
    
    ## program ends