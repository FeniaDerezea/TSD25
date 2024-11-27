#############################################################################################
### p-value code for detecting conflict between the data at threshold 10 and all the rest ###
#############################################################################################
library(R2WinBUGS);library(MASS);library(MCMCpack)
require(mcmcplots)
#############################################################################################

#read data at threshold 10 only
Study<-1:10
FP<-c(103,417,1159,180,7176,489,699,100,151,742)
TP<-c(10,16,67,34,461,51,75,12,25,59)  
N1<-c(155,2875,5267,694,33180,3408,3506,346,722,4470)
N2<-c(11,17,74,38,514,54,90,12,28,73)

#read data at thr 100
#Study<-1:5
#FP<-c(33,78,263,1628,208)
#TP<-c(6,10,53,329,58)  
#N1<-c(155,2875,5267,33180,3506)
#N2<-c(11,17,74,514,90)

rc<-cbind(TP,FP)
# matrix with No of healthy in 1st column and No of diseased in 2nd
Nb<-cbind(N2,N1)
#calculate No of rows (which is equal to No of studies)
nsb<-nrow(rc)

#read data at all the other thresholds
data <- read.csv("OC_Sensor_excl10.csv",na="NA",header = TRUE)
#data <- read.csv("OC_Sensor_excl100.csv",na="NA",header = TRUE)

## bring data in form needed for the model
ns <- length(data$ID)

N <- as.matrix(data[, grepl("N", names(data))])
N <- cbind(N[,2],N[,1])
C <- as.matrix(data[, grepl("C", names(data))])
tp <- as.data.frame(data[, grepl("tp", names(data))])
fp <- as.data.frame(data[, grepl("fp", names(data))])

# table with number of thresholds per study
Tc <- data$T
# remove any consecutive zero FP counts
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

x <- array(NA, c(ns, 2, max(Tc[,1])))
for(i in 1:ns){
  for(t in 1:max(Tc[,1])){
    x[i, 2, t] <- fp[i,t]
    x[i, 1, t] <- tp[i,t]
  }
}

dataList = list(ns = ns, Tc = Tc, N = N, C = C, x = x, rc=rc, Nb=Nb, nsb=nsb)
data.names<-c(names(dataList) )




inits1<-list(  
  list(mean = c(5, 3,log(1), log(1)), 
       sd = c(1, 1, 1, 1),
        lambda=-0.8, theta =c(0,-1) , sdb=c(1,1),rhob=0.1) ,
  list(mean = c(6, 0,log(0.5), log(0.5)), 
       sd = c(0.1, 0.1, 0.1, 0.1),
        lambda=-0.5, theta =c(1,1) , sdb=c(0.6,0.6),rhob=0.2) ,
  list(mean = c(4, 2,log(1.2), log(1.2)), 
       sd = c(0.1, 0.2, 0.1, 0.2),
       lambda=-0.3, theta =c(0.5,0.5) , sdb=c(1.5,1.5),rhob=0.3)
) 

mymonitoredparamslist <- c( "pr_sum","sumtprb","sumfprb","pvalsens","pvalfpf","theta")#

# usual directory
winbugs.dir <- "C:/MyDirectory/WinBUGS14/"
##Run the model 
set.seed(213)
modelp.sim <- bugs(data.names, inits1, model.file = "models_pval.txt", parameters = mymonitoredparamslist ,
                   n.chains = 3, n.iter = 600000, n.burnin=300000, n.thin=10,  bugs.directory = winbugs.dir, debug=FALSE)

print(modelp.sim,3)

#p-values
prse<-modelp.sim$mean$pvalsens
pvalue_sens<-2*min(prse,1-prse)

prfp<-modelp.sim$mean$pvalfpf
pvalue_fpf<-2*min(prfp,1-prfp)
