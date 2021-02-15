//#######################################################
//### State-space one state 2 age groups model ##########
//#######################################################

///////////////////////// DATA //////////////////////////
  data {
    int<lower = 0> N;       // number of data
    int<lower = 0> A;       // number of age groups
    matrix[N, A] Y;         // response vector
    
    // hyperparameters
    real alpha_m0;
    real alpha_phi;

    real alpha_m0_var;
    real alpha_phi_var;

    real a_sigma2;
    real b_sigma2;
  }

//////////////////// PARAMETERS /////////////////////////
  parameters {
    vector[A] m0;
    vector[A] phi;
    vector[A] sigma2;
    real mu_m0;
    real mu_phi;
    real sigma2m0;
    real sigma2phi;
  }


////////////////// MODEL ////////////////////////////////
  model {
    
    // Likelihood     
    for (t in 2:N){
      for (a in 1:A){
        Y[t,a] ~ normal(m0[a] + phi[a] * Y[t - 1,a], pow(sigma2[a], 0.5));
      } 
    }
    for (a in 1:A){
      m0[a] ~ normal(mu_m0, sigma2m0);
      phi[a] ~ normal(mu_phi, sigma2phi);
      sigma2[a] ~ inv_gamma(a_sigma2, b_sigma2);
    }
    
    mu_m0 ~ normal(alpha_m0, alpha_m0_var);
    sigma2m0 ~ inv_gamma(a_sigma2, b_sigma2);
    
    mu_phi ~ normal(alpha_phi, alpha_phi_var);
    sigma2phi ~ inv_gamma(a_sigma2, b_sigma2);
  }

generated quantities {
  matrix[N,A] log_lik;
  for (t in 1:N) {
    for (a in 1:A) {
      if(t == 1){
        log_lik[t,a] = normal_lpdf(Y[t,a] | Y[t,a], 
                                   pow(sigma2[a], 0.5)); 
      } else{
        log_lik[t,a] = normal_lpdf(Y[t,a] | m0[a] + phi[a] * Y[t - 1,a], 
                                   pow(sigma2[a], 0.5));
      }
    }
  }
}
