data {
  int<lower=0> cutoff;
  int<lower=1> n_days;
  int<lower=1> days[n_days];
  int<lower=1> n_d_days;
  int<lower=1> d_days[n_d_days];
  real<lower=0> N;      // population
  int<lower=0> d[n_d_days];              // deaths
  vector<lower=0>[cutoff+1] tau_rev; // time to death distribution
  real<lower=0> ifr_min;
  real<lower=0> ifr_max;
  
  int<lower=1> n_v; // number of viral tests administered
  int<lower=1> n_s; // number of serology tests administered
  real<lower=0, upper=1> vp_hat;  // viral prevalence
  real<lower=0, upper=1> sp_hat; //seroprevalence
  int T1;  // start of phase 1 testing
  int T2;  // end of phase 1 testing
  
  int<lower=0> n_per;
  int<lower=0> n_test_days;
  int<lower=0> test_days[n_test_days];
  vector<lower=0>[n_per] cc;
  vector<lower=0>[n_per+1] tests;
  int<lower=0> per_ends[n_per+1];
  
  real<lower=0> gamma_inv_IN_mean;
  real<lower=0> gamma_inv_IN_sd;
  
  real<lower=0> phi_IN_mean;
  real<lower=0> phi_IN_sd;
}
parameters {
  real<lower=5.5, upper = 11.5> gamma_inv;
  real<lower=0> beta[n_days];
  ////real<lower=0> beta_tilde[n_days];
  real sigma_tilde;
  real<lower=ifr_min,upper=ifr_max> ifr;
  real y0_tilde[2];
  real<lower = 0,upper=2> phi;
  real<lower = 0> sigma_phi;
}
transformed parameters{
  ////real<lower=0> beta[n_days];
  matrix<lower=0>[n_days+1, 3] y;
  vector<lower=0>[n_days] nu;
  real<lower=0,upper=1> sigma = Phi(sigma_tilde);
  real<lower=0, upper=1> vp;
  real<lower=0,upper=1> sp;
  real<lower=0> vp_se;
  real<lower=0> sp_se;
  real<lower=1.0/11.5,upper=1.0/5.5> gamma = 1/gamma_inv;
  real<lower=0> death_par[n_days];
  vector<lower=0>[n_per] cc_mean;
  vector<lower=0>[n_per] cc_sd;
  {
    y[1,1] = N*0.99 + 0.01*N*Phi(y0_tilde[1]);
    y[1,2] = 0.01*N*(1-Phi(y0_tilde[1]))*Phi(y0_tilde[2]);
    y[1,3] = N - y[1,1] - y[1,2];
    ////beta[1] = Phi(beta_tilde[1]);
    for(t in 1:n_days){
      ////if(t >= 2){
      ////  beta[t] = beta[t-1] + sigma * beta_tilde[t];
      ////}
      y[t+1,1] = y[t,1] * (1 - beta[t] * y[t,2] / N);
      y[t+1,2] = y[t,2] * (1 + beta[t] * y[t,1] / N - gamma);
      y[t+1,3] = y[t,3] + gamma * y[t,2];
      nu[t] = y[t,1] - y[t+1,1];
      if(t <= cutoff){
      death_par[t]=ifr*sum(tau_rev[(cutoff-t+2):(cutoff+1)] .*nu[1:t])+ifr*sum(tau_rev[1:(cutoff-t+1)])*(y[1,2]+y[1,3]);
      }else{
      death_par[t]=ifr*sum(tau_rev .*nu[(t-cutoff):t]);
      }
    }
  }
  
  for(t in 1:n_per){
    cc_mean[t] = phi*(sqrt(sum(tests[1:(t+1)])/N) * (N-y[per_ends[t+1]+1,1]) - sqrt(sum(tests[1:t])/N) * (N-y[per_ends[t]+1,1]));
    cc_sd[t] = sigma_phi * sqrt(tests[t+1]/N);
  }  

  vp = mean(y[(T1+1):(T2+1),2])/N;
  vp_se = sqrt(vp*(1-vp)/n_v);
  
  sp = mean(y[(T1+1):(T2+1),3])/N;
  sp_se = sqrt(sp*(1-sp)/n_s);
}
model {
  //priors
  beta[1] ~ uniform(0,1);
  beta[2:n_days] ~ normal(beta[1:(n_days-1)], sigma);
  ////beta_tilde ~ normal(0,1);
  gamma_inv ~ normal(gamma_inv_IN_mean, gamma_inv_IN_sd) T[5.5,11.5];
  sigma_tilde ~ normal(0,1);
  ifr ~ uniform(ifr_min,ifr_max);
  y0_tilde ~ normal(0,1);
  phi ~ normal(phi_IN_mean,phi_IN_sd) T[0,2];
  sigma_phi ~ uniform(0,3);
  
  //likelihoods
  d ~ poisson(death_par[1:n_d_days]);
  vp_hat ~ normal(vp,vp_se) T[0,1];
  sp_hat ~ normal(sp,sp_se) T[0,1];
  100*cc/N ~ normal(100/N * cc_mean, cc_sd);
}
