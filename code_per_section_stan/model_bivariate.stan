data {
  int<lower=1> ns; // number of studies
  
  array[ns,2] int<lower=0> N; // number of healthy and diseased individuals
  array[ns,2] int<lower=0> r; // TP,FP,FN,TN data
}

parameters {
  vector[2] theta; // mean logit(TPF) and logit(FPF)
  vector<lower=0>[2] sd; // between-studies SD in logit(TPF) and logit(FPF)
  real<lower=-1,upper=1> rho; // between-studies correlation
  
  matrix[ns,2] delta; // random effects of sensitivity and FPF
  
  real predtpr; // logit sensitivity predictive value
  real predfpf; // logit fpf predictive value
}

transformed parameters {
  matrix[ns,2] p; // sensitivity and FPF for the model
  for (i in 1:ns){
    p[i,1] = inv_logit(delta[i,1]);
    p[i,2] = inv_logit(delta[i,2]);
  }
  
  vector[2] prec; // logit(TPF) and logit(FPF) precision
  prec[1] = sd[1]^(-2);
  prec[2] = sd[2]^(-2);
  
  vector[ns] condmean; // conditional mean
  for (i in 1:ns){
    condmean[i] = theta[2] + rho * (sd[2]/sd[1]) * (delta[i,1]-theta[1]);
  }
  
  real condprec; // conditional precision
  condprec = 1/((1-rho^2)*sd[2]^2);
  
  real<lower=0,upper=1> sumtpf; // summary TPF
  sumtpf = exp(theta[1])/(1+exp(theta[1]));
  
  real<lower=0,upper=1> sumfpf; // summary FPF
  sumfpf = exp(theta[2])/(1+exp(theta[2]));
  
  real<lower=0,upper=1> spec; // summary specificity
  spec = 1-sumfpf;
  
  real beta; // beta parameter for the HSROC model
  beta = log(sd[2]/sd[1]);
  
  real Theta; // Theta threshold parameter for the HSROC model
  Theta = ((sqrt(sd[2]/sd[1]))*theta[1]+sqrt((sd[1]/sd[2]))*theta[2])*0.5;
  
  real Lambda; // Lambda accuracy parameter for the HSROC model
  Lambda = (sqrt(sd[2]/sd[1]))*theta[1] - (sqrt(sd[1]/sd[2]))*theta[2];
  
  real vartheta; // variance of theta HSROC parameter
  vartheta = 0.5*(sd[1]*sd[2] + sd[1]*sd[2]*rho);
  
  real varalpha; // variance of alpha HSROC parameter
  varalpha = 2*(sd[1]*sd[2] - sd[1]*sd[2]*rho);
  
  real sdtheta; // standard deviation of theta HSROC parameter
  sdtheta = sqrt(vartheta);
  
  real sdalpha; // standard deviation of alpha HSROC parameter
  sdalpha = sqrt(varalpha);
  
  real condmeanpred; // predicitive conditional mean
  condmeanpred = theta[2] + rho*(sd[2]/sd[1])*(predtpr-theta[1]);
  
  real<lower=0,upper=1> pred_tpr; // predictive sensitivity
  pred_tpr = inv_logit(predtpr);
  
  real<lower=0,upper=1> pred_fpf; // predictive fpf
  pred_fpf = inv_logit(predfpf);
  
  real<lower=0,upper=1> predspec; // predictive specificity
  predspec = 1-pred_fpf;
  
  real predlambda; // predictive Lambda
  predlambda = (sqrt(sd[2]/sd[1]))*predtpr - (sqrt(sd[1]/sd[2]))*predfpf;
}

model {
  // priors
  theta[1] ~ logistic(0,1); // prior for the mean logit(TPF)
  theta[2] ~ logistic(0,1); // prior for the mean logit(FPF)
  
  sd[1] ~ uniform(0,5); // Prior for between-studies SD in logit(TPF)
  sd[2] ~ uniform(0,5); // Prior for between-studies SD in logit(FPF)
  
  rho ~ uniform(-1,1); // Prior for between-studies correlation
  
  predtpr ~ normal(theta[1], 1/sqrt(prec[1])); //logit sensitivity predicitive value
  predfpf ~ normal(condmeanpred, 1/sqrt(condprec)); // logit fpf predicitive value
  
  // likelihoods
  for (i in 1:ns){
    r[i,1] ~ binomial(N[i,1],p[i,1]); // binomial likelihood for TP counts
    r[i,2] ~ binomial(N[i,2],p[i,2]); // binomial likelihood for FP counts
  }
  
  // random effects distributions
  for (i in 1:ns){
    delta[i,1] ~ normal(theta[1], 1/sqrt(prec[1])); // random effects on sensitivity
    delta[i,2] ~ normal(condmean[i], 1/sqrt(condprec)); // random effects on FPF
  }
}