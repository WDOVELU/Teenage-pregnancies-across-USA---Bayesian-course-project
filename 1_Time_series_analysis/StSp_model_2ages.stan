//#######################################################
//### State-space one state 2 age groups model ##########
//#######################################################

///////////////////////// DATA //////////////////////////
data {
	int<lower = 0> N;       // number of data
	int<lower = 0> A;       // number of age groups
	matrix[N, A] Y;         // response vector
	real sigma2phi1;
	real sigma2phi2;
	real a_sigma2;
	real b_sigma2;
	real sigma2m0;
}

//////////////////// PARAMETERS /////////////////////////
parameters {
	vector[A] phi;
	vector[A] sigma2;
	real m0;
}


////////////////// MODEL ////////////////////////////////
model {

	// Likelihood     
	for (t in 2:N) {
	  for (a in 1:A){
	  	Y[t,a] ~ normal(m0 + phi[a] * Y[t - 1,a], pow(sigma2[a], 0.5));
	  } 
	}
	m0 ~ normal(0, sigma2m0);
	phi[1] ~ normal(0, sigma2phi1);
	phi[2] ~ normal(0, sigma2phi2);	
	sigma2[1] ~ inv_gamma(a_sigma2, b_sigma2);
	sigma2[2] ~ inv_gamma(a_sigma2, b_sigma2);
}

generated quantities {
  	matrix[N,A] log_lik;
  	for (t in 1:N) {
  	  for (a in 1:A) {
  	    if(t == 1){
  	      log_lik[t,a] = normal_lpdf(Y[t,a] | Y[t,a], 
  	                                          pow(sigma2[a], 0.5)); 
  	    } else{
  	      log_lik[t,a] = normal_lpdf(Y[t,a] | m0 + phi[a] * Y[t - 1,a], 
    		                                    pow(sigma2[a], 0.5));
  	    }
  	  }
  	}
}
