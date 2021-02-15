//#######################################################
//### State-space one state 2 age groups model ##########
//#######################################################

///////////////////////// DATA //////////////////////////
data {
	int<lower = 0> N;       // number of data
	int<lower = 0> A;       // number of age groups
	matrix[N, A] Y;         // response vector
	real sigma2m0_1;
	real sigma2m0_2;
	real sigma2phi;
	real a_sigma2;
	real b_sigma2;
}

//////////////////// PARAMETERS /////////////////////////
parameters {
	vector[A] m0;
	real phi;
	vector[A] sigma2;
}


////////////////// MODEL ////////////////////////////////
model {

	// Likelihood     
	for (t in 2:N){
	  for (a in 1:A){
	  	Y[t,a] ~ normal(m0[a] + phi * Y[t - 1,a], pow(sigma2[a], 0.5));
	  } 
	}
	m0[1] ~ normal(0, sigma2m0_1);
	m0[2] ~ normal(0, sigma2m0_2);
	phi ~ normal(0, sigma2phi);
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
  	      log_lik[t,a] = normal_lpdf(Y[t,a] | m0[a] + phi * Y[t - 1,a], 
    		                                    pow(sigma2[a], 0.5));
  	    }
  	  }
  	}
}
