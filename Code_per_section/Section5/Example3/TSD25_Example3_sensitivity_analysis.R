################################################################################
### TSD25- Example 3: sensitivity analysis #####################################
################################################################################
#load required libraries
library(R2WinBUGS);library(MASS);library(MCMCpack)
require(mcmcplots)

################################################################################
#Fit the bivariate model at data only from threshold 10 ########################
#########################################################################################

# the relevant data
Study<-1:10
FP<-c(103,417,1159,180,7176,489,699,100,151,742)
TP<-c(10,16,67,34,461,51,75,12,25,59)  
N1<-c(155,2875,5267,694,33180,3408,3506,346,722,4470)
N2<-c(11,17,74,38,514,54,90,12,28,73)

r<-cbind(TP,FP)
# matrix with No of healthy in 1st column and No of diseased in 2nd
N<-cbind(N2,N1)
#calculate No of rows (which is equal to No of studies)
ns<-nrow(r)

dat<-list(r=r,N=N)
data.names<-c(names(dat) ,'ns')
#create vector with names of parameters to monitor
parameter.names <- c( 'theta','rho','sd','sumtpr','sumfpr','spec','Theta','delta','beta','ppv','npv','predtpr','predfpf','logitsp','deltasp','Lambda','vartheta','varalpha')


# defining the directory of WinBUGS, this should be replaced with the full path to the directory WinBUGS is located
winbugs.dir <- "C:/My_winbugs_directory/WinBUGS14/"

# select No of chains, iterations and burnin
n.chain=3
n.iters=100000
n.burn=20000
n.th=3
#create list containing initial value for 3 chains
inits1<-list(  
  list(theta =c(0,-1) , sd=c(1,1) ),
  list(theta =c(1,1) , sd=c(0.6,0.6) ),
  list(theta =c(0.5,0.5) , sd=c(1.5,1.5) )
)

#run model through WinBUGS and safe output
model0.sim <- bugs( data.names, inits1, model.file = "model1.txt", parameters = parameter.names,
                    n.chains = n.chain, n.iter = n.iters, n.burnin=n.burn, n.thin=n.th,  bugs.directory = winbugs.dir, debug=FALSE)

# print model summary rounded to 3 decimals
print(model0.sim,3)
#trace plot for saved parameters
mcmcplot(model0.sim)


################################################################################
#Fit the bivariate model at data only from threshold 100 ########################
################################################################################

Study<-1:5
FP<-c(33,78,263,1628,208)
TP<-c(6,10,53,329,58)  
N1<-c(155,2875,5267,33180,3506)
N2<-c(11,17,74,514,90)

r<-cbind(TP,FP)
# matrix with No of healthy in 1st column and No of diseased in 2nd
N<-cbind(N2,N1)
#calculate No of rows (which is equal to No of studies)
ns<-nrow(r)

dat<-list(r=r,N=N)
data.names<-c(names(dat) ,'ns')
#create vector with names of parameters to monitor
parameter.names <- c( 'theta','rho','sd','sumtpr','sumfpr','spec','Theta','beta','predtpr','predfpf','logitsp','Lambda','vartheta','varalpha')


# defining the directory of WinBUGS, this should be replaced with the full path to the directory WinBUGS is located
winbugs.dir <- "C:/My_winbugs_directory/WinBUGS14/"

# select No of chains, iterations and burnin
n.chain=3
n.iters=100000
n.burn=20000
n.th=3
#create list containing initial value for 3 chains
inits1<-list(  
  list(theta =c(0,-1) , sd=c(1,1) ),
  list(theta =c(1,1) , sd=c(0.6,0.6) ),
  list(theta =c(0.5,0.5) , sd=c(1.5,1.5) )
)

#run model through WinBUGS and safe output
model0.sim <- bugs( data.names, inits1, model.file = "model1.txt", parameters = parameter.names,
                    n.chains = n.chain, n.iter = n.iters, n.burnin=n.burn, n.thin=n.th,  bugs.directory = winbugs.dir, debug=FALSE)

# print model summary rounded to 3 decimals
print(model0.sim,3)
#trace plot for saved parameters
mcmcplot(model0.sim)


