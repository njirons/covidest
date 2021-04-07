data{
  int<lower=1> n_samples_IN;
  int<lower=1> n_samples_OH;
  
  // posterior phi samples
  real<lower=0> phi_IN[n_samples_IN];
  real<lower=0> phi_OH[n_samples_OH];
  real<lower=0> phi_IN_sd; //posterior std. dev. of phi
  real<lower=0> phi_OH_sd; //posterior std. dev. of phi
  real<lower=0> sigma_phi_IN[n_samples_IN];
  real<lower=0> sigma_phi_OH[n_samples_OH];
  real<lower=0> sigma_phi_IN_sd; //posterior std. dev. of sigma_phi
  real<lower=0> sigma_phi_OH_sd; //posterior std. dev. of sigma_phi
  
  // posterior ifr samples
  real<lower=0> ifr_min;
  real<lower=0> ifr_max;  
  real<lower=ifr_min,upper=ifr_max> ifr_IN[n_samples_IN]; //posterior IFR samples
  real<lower=ifr_min,upper=ifr_max> ifr_OH[n_samples_OH]; //posterior IFR samples
  real<lower=0> sigma_ifr_IN; //posterior std. dev. of IFR
  real<lower=0> sigma_ifr_OH; //posterior std. dev. of IFR
}

parameters{
  // test positivity parameters
  ////real<lower=0,upper=1> mu_phi_IN;
  ////real<lower=0,upper=1> mu_phi_OH;
  ////real<lower=0,upper=1> mu_phi;
  ////real<lower=0> tau_phi;
  
  // ifr parameters
  real<lower=ifr_min,upper=ifr_max> mu_ifr_IN; // IN IFR
  real<lower=ifr_min,upper=ifr_max> mu_ifr_OH; // OH IFR
  real<lower=ifr_min,upper=ifr_max> mu_ifr; //hyperprior mean of IFR
  real<lower=0> tau_ifr; //hyperprior precision of IFR
}

transformed parameters{
  real<lower=0> sigma_ifr = 1/sqrt(tau_ifr);
}

model{
  //// preferential testing BHM
  //likelihood
  ////phi_IN ~ normal(mu_phi_IN, phi_IN_sd);
  ////phi_OH ~ normal(mu_phi_OH, phi_OH_sd);
  ////sigma_phi_IN ~ normal(mu_sigma_phi_IN, sigma_phi_IN_sd);
  ////sigma_phi_OH ~ normal(mu_sigma_phi_OH, sigma_phi_OH_sd);
  
  //priors
  ////mu_phi_IN ~ normal(mu_phi,tau_phi) T[0,1];
  ////mu_phi_OH ~ normal(mu_phi,tau_phi) T[0,1];
  
  ////mu_sigma_phi_IN ~ normal(mu_sigma_phi, )
  
  //hyperpriors
  ////mu_phi ~ uniform(0,1);
  ////tau_phi ~ uniform(0,1);
  
  //// ifr BHM
  //likelihood
  ifr_IN ~ normal(mu_ifr_IN, sigma_ifr_IN);
  ifr_OH ~ normal(mu_ifr_OH, sigma_ifr_OH);
  
  //priors
  mu_ifr_IN ~ normal(mu_ifr,sigma_ifr) T[ifr_min,ifr_max];
  mu_ifr_OH ~ normal(mu_ifr,sigma_ifr) T[ifr_min,ifr_max];
  
  //hyperpriors
  mu_ifr ~ normal(0.00675,0.000725) T[ifr_min,ifr_max];
  tau_ifr ~ exponential(0.0000008206947);
}
