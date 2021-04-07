State <- "IN"
load(paste0("data/",State,"_ctp.RData"))

# specifying parameters
days <- ctp_data$days
test_days = ctp_data$shifted_test_days

T1 <- which(days == yday("2020-04-25"))
T2 <- which(days == yday("2020-04-29"))

# phase I testing results
vp_hat <- 1.74/100; vp_se <- (2.5-1.1)/(200*1.96); n_v <- 3605;
sp_hat <- 1.09/100; sp_se <- (1.5-0.8)/(200*1.96); n_s <- 3518;

####################### Run Stan
IN_dat <- list(n_days = ctp_data$n_days, 
               n_d_days = ctp_data$n_d_days, 
               d_days = ctp_data$d_days,
               days = ctp_data$days, 
               N = ctp_data$N, d = ctp_data$d, 
               tau_rev = rev(ctp_data$ttd),
               ifr_min = ctp_data$ifr_min, ifr_max = ctp_data$ifr_max,
               n_v = n_v, n_s = n_s,
               vp_hat = vp_hat, sp_hat = sp_hat,
               T1= T1, T2=T2,
               n_test_days = ctp_data$n_test_days,
               test_days = test_days,
               cc = ctp_data$cum_cc_per,
               tests = ctp_data$cum_tests_per,
               cutoff = ctp_data$cutoff, 
               n_per = ctp_data$n_per,
               per_ends = ctp_data$per_ends)
init_fn <- function() {
  list(gamma_inv = 11, beta = rep(0.1,ctp_data$n_days),
       ifr=0.008, phi=0.5,sigma_phi=0.5)
}
fit <- stan(file = 'code/Indiana_final.stan', init = init_fn, 
            data = IN_dat, chains = 4, iter = 10000,
            control = list(max_treedepth = 15,adapt_delta = 0.995))
            # data = IN_dat, chains = 4, iter = 20000,
            # control = list(max_treedepth = 15,adapt_delta = 0.99))
saveRDS(fit, "fits/IN_fit.rds")

sims <- rstan::extract(fit)
save(sims,file ="fits/IN_sims.RData")

####### fit truncated normal dist to posterior gamma_inv
tn_mll <- function(mu,sigma){
  mll <- -mean(dnorm(sims$gamma_inv,mean=mu,sd=sigma,log=TRUE)) +
    log(pnorm(11.5,mean=mu,sd=sigma)-pnorm(5.5, mean=mu,sd=sigma))
  return(mll)  
}

tn_fit <- mle(tn_mll,start = c(mean(sims$gamma_inv),
                               sd(sims$gamma_inv)))

gamma_inv_IN <- list(mean = coef(tn_fit)[1], 
                     sd = coef(tn_fit)[2])
save(gamma_inv_IN,file = "data/gamma_inv_IN.RData")

####### fit normal dist to posterior phi
phi_IN <- list(mean = mean(sims$phi), 
                     sd = sd(sims$phi))
save(phi_IN,file = "data/phi_IN.RData")
