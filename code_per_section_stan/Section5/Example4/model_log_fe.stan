data {
  int<lower=1> ns;
  int<lower=1> max_T;
  array[ns,2] int Tc;
  array[ns,2] int N;
  matrix[ns,max_T] C;
  array[ns,2,max_T] int x;
}

transformed data {
  array[ns,2,max_T] int n;
  for (i in 1:ns){
    for (j in 1:2){
      n[i,j,1] = N[i,j];
      for (t in 2:Tc[i,j]){  // Note (1) (see below)
        n[i,j,t] = x[i,j,t-1];
      }
    }
  }
}

parameters {
  vector[4] mean;
}

transformed parameters {
  array[ns,2,max_T] real d;
  for (i in 1:ns){
    for (j in 1:2){
		  for (t in 1:Tc[i,j]){
		    if (C[i,t] == 0){
		      d[i,j,t] = (mean[j]) / exp(mean[2+j]);
		    }
		    else {
			    d[i,j,t] = (mean[j] - log(C[i,t]) ) / exp(mean[2+j]);
		    }
		  }
    }
  }
  
  array[ns,2,max_T] real pr; // Note (3)
  for (i in 1:ns){
    for (j in 1:2){
      for (t in 1:Tc[i,j]){
        pr[i,j,t] = inv_logit(fmin(10.0, fmax(-10.0, d[i,j,t])));
      }
    }
  }
  
  array[ns,2,max_T] real p;
  for (i in 1:ns){
	  for (j in 1:2){
		  p[i,j,1] = pr[i,j,1]; 
		  for (t in 2:Tc[i,j]){
		  	p[i,j,t] =  pr[i,j,t] / pr[i,j,t-1];	 
			}
	  }
  }
  
  array[ns,2,max_T] real xhat; // Fitted values
  for (i in 1:ns){
    for (j in 1:2){
		  for (t in 1:Tc[i,j]){
			  xhat[i,j,t] = p[i,j,t]*n[i,j,t];
		  }
    }
  }
  
  array[ns,2,max_T] real dev; // Residual deviance contribution
  for (i in 1:ns){
    for (j in 1:2){
      for (t in 1:Tc[i,j]){
        if (x[i,j,t] == n[i,j,t]){
          dev[i,j,t] = 2*(x[i,j,t]*(log(x[i,j,t]) - log(xhat[i,j,t])));
        }
        else {
          dev[i,j,t] = 2*(x[i,j,t]*(log(x[i,j,t]) - log(xhat[i,j,t])) + (n[i,j,t] - x[i,j,t])*(log(n[i,j,t] - x[i,j,t]) - log(n[i,j,t] - xhat[i,j,t])));  
        }
      }
    }
  }
  
  vector[ns] rd; // Residual deviance study i
  for (i in 1:ns){
    rd[i] = sum(dev[i,1,1:Tc[i,1]]) + sum(dev[i,2,1:Tc[i,2]]);	
  }
  
  real resdev; // Total residual deviance
  resdev = sum(rd[]);    
}

model {
  // Priors
  for (r in 1:4){
    mean[r] ~ normal(0,1/sqrt(0.001));
  }
  
  // Likelihood
  for (i in 1:ns){
    for (j in 1:2){
		  for (t in 1:Tc[i,j]){
        x[i,j,t] ~ binomial(n[i,j,t],p[i,j,t]);
		  }
    }
  }
}