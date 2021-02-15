//#######################################################
//### State-space one state 2 age groups model ##########
//#######################################################

// d_tj = alpha + phi1_j * d_t-1,j + phi2_j * d_t-2,j 
//              + beta1_j * eps_t-1,j +beta2_j * eps_t-2,j+ beta1_j * eps_t-3,j
//              + eps_t




///////////////////////// DATA //////////////////////////
data {
	int<lower = 0> N;       // number of data
	int<lower = 0> A;       // number of age groups
	matrix[N, A] Y;         // response vector
	real a_sigma2;
	real b_sigma2;
	real sigma2m0;
}

//////////////////// PARAMETERS /////////////////////////
parameters {
	vector[A] phi1;
    vector[A] phi2;
    vector[A] beta1;
    vector[A] beta2;
    vector[A] beta3;
	vector[A] sigma2;
	real m0;
}

transformed parameters {
	matrix[N, A] epsilon;        
	matrix[N, A] nu;   
    matrix[N,A] log_lik;

	for (a in 1:A){
	   nu[1,a] = m0+ phi1[a] * m0 + phi2[a] * m0;
       epsilon[1,a]=Y[1,a]-nu[1,a];
       nu[2,a] = m0+ phi1[a] * Y[1,a] + phi2[a] * m0 + beta1[a] * epsilon[1,a];
       epsilon[2,a]=Y[2,a]-nu[2,a];
       nu[3,a] = m0+ phi1[a] * Y[2,a] + phi2[a] * Y[1,a] + beta1[a] * epsilon[2,a]+ beta2[a] * epsilon[1,a];
       epsilon[3,a]=Y[3,a]-nu[3,a];
	}
    for (t in 4:N) {
	  for (a in 1:A){
	  	nu[t,a] = m0+ phi1[a] * Y[t-1,a] + phi2[a] * Y[t-2,a]  + beta1[a] * epsilon[t-1,a] + beta2[a] * epsilon[t-2,a]+beta3[a] * epsilon[t-3,a];
        epsilon[t,a]=Y[t,a]-nu[t,a];
	  } 
	}


  	for (t in 4:N) {
  	  for (a in 1:A) {
  	    if(t == 1){
  	      log_lik[t,a] = normal_lpdf(Y[t,a] | Y[t,a], 
  	                                          pow(sigma2[a], 0.5)); 
  	    } else{
  	      log_lik[t,a] = normal_lpdf(Y[t,a] | m0 + phi1[a] * Y[t-1,a] + phi2[a] * Y[t-2,a]  + beta1[a] * epsilon[t-1,a] + beta2[a] * epsilon[t-2,a]+beta3[a] * epsilon[t-3,a], pow(sigma2[a], 0.5));
  	    }
  	  }
  	}
}


////////////////// MODEL ////////////////////////////////
model {

	// Likelihood     
	for (t in 4:N) {
	  for (a in 1:A){
	  	Y[t,a] ~ normal(m0 + phi1[a] * Y[t-1,a] + phi2[a] * Y[t-2,a]  + beta1[a] * epsilon[t-1,a] + beta2[a] * epsilon[t-2,a]+beta3[a] * epsilon[t-3,a], pow(sigma2[a], 0.5));
	  } 
	}


	m0 ~ normal(0, sigma2m0);
	phi1[1] ~ normal(0, 2);
	phi1[2] ~ normal(0, 2);	
    phi2[1] ~ normal(0, 2);
	phi2[2] ~ normal(0, 2);
    beta1[1] ~ normal(0,2);
	beta1[2] ~ normal(0,2);
    beta2[1] ~ normal(0,2);
	beta2[2] ~ normal(0,2);
    beta3[1] ~ normal(0,2);
	beta3[2] ~ normal(0,2);
//    epsilon ~ normal(0,e_sigma);
//    e_sigma ~ 
	sigma2[1] ~ inv_gamma(a_sigma2, b_sigma2);
	sigma2[2] ~ inv_gamma(a_sigma2, b_sigma2);
}




// generated quantities {
//   	matrix[N,A] log_lik;
//   	for (t in 4:N) {
//   	  for (a in 1:A) {
//   	    if(t == 1){
//   	      log_lik[t,a] = normal_lpdf(Y[t,a] | Y[t,a], 
//   	                                          pow(sigma2[a], 0.5)); 
//   	    } else{
//   	      log_lik[t,a] = normal_lpdf(Y[t,a] | m0 + phi1[a] * Y[t-1,a] + phi2[a] * Y[t-2,a]  + beta1[a] * epsilon[t-1,a] + beta2[a] * epsilon[t-2,a]+beta3[a] * epsilon[t-3,a], pow(sigma2[a], 0.5));
//   	    }
//   	  }
//   	}
// }
