### load data
abbrev <- read_csv("data/state_abbrev.csv")
states <- c(abbrev$abbreviation,"DC")
load("data/gamma_inv.RData")
load("data/phi.RData")
load("data/tau_ifr.RData")
load("data/mu_ifr.RData")

# come back to GA (10)? used len_per=3 for that one. lots of divergent transitions.
# also used len_per=3 for OK (36) and OR (37)

for(s in (1:26)[-c(14,35)]){ #Run for all states, other than IN,OH
  State <- states[s]
  
  ### load state data
  load(paste0("data/",State,"_ctp.RData"))

  ####################### Run Stan
  state_data <- list(n_days = ctp_data$n_days, 
                     n_d_days = ctp_data$n_d_days, 
                     d_days = ctp_data$d_days,
                     days = ctp_data$days, 
                     N = ctp_data$N, d = ctp_data$d, 
                     tau_rev = rev(ctp_data$ttd),
                     ifr_min = ctp_data$ifr_min, ifr_max = ctp_data$ifr_max,
                     n_test_days = ctp_data$n_test_days,
                     test_days = ctp_data$shifted_test_days,
                     cc = ctp_data$cum_cc_per,
                     tests = ctp_data$cum_tests_per,
                     cutoff = ctp_data$cutoff, 
                     n_per = ctp_data$n_per,
                     per_ends = ctp_data$per_ends,
                     gamma_inv_mean = gamma_inv$mean,
                     gamma_inv_sd = gamma_inv$sd,
                     phi_mean = phi$mean,
                     phi_sd = phi$sd,
                     tau_ifr_shape = tau_ifr$shape,
                     tau_ifr_rate = tau_ifr$rate,
                     mu_ifr_mean = mu_ifr$mean,
                     mu_ifr_sd = mu_ifr$sd)
  init_fn <- function() {
    list(gamma_inv = 11, beta = rep(0.1,ctp_data$n_days),
         ifr=0.008, phi=0.5,sigma_phi=0.5)
  }
  fit <- stan(file = "code/OtherStates_final.stan", init = init_fn,
              data = state_data, chains = 4, iter = 10000,
              control = list(max_treedepth = 15,adapt_delta = 0.995))           
              #data = state_data, chains = 4, iter = 20000,
              #control = list(max_treedepth = 15,adapt_delta = 0.995))
  #print(warnings())
  
  saveRDS(fit, paste0("fits/",State,"_fit.rds"))
  
  sims <- rstan::extract(fit)
  save(sims, file = paste0("fits/",State,"_sims.RData"))
  
  print(paste0(State," complete, ",round(100*s/51),"%"))
  
  rm(fit,sims,ctp_data,state_data,init_fn,State)
}

