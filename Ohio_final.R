State <- "OH"
load(paste0("data/",State,"_ctp.RData"))

load("data/gamma_inv_IN.RData")
load("data/phi_IN.RData")

days <- ctp_data$days
test_days = ctp_data$shifted_test_days

T1 <- which(days == yday("2020-07-09"))
T2 <- which(days == yday("2020-07-28"))

vp_hat <- 0.9/100; vp_se <- (2.0-0.1)/(200*1.96); n_v <- 727
#sp_hat <- 1.5/100; sp_se <- (2.9-0.3)/(200*1.96);
sp_hat <- 1.3/100; sp_se <- (2.7-0.2)/(200*1.96); n_s <- 667

####################### Run Stan
OH_dat <- list(n_days = ctp_data$n_days, 
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
               per_ends = ctp_data$per_ends,
               gamma_inv_IN_mean = gamma_inv_IN$mean,
               gamma_inv_IN_sd = gamma_inv_IN$sd,
               phi_IN_mean = phi_IN$mean,
               phi_IN_sd = phi_IN$sd)
init_fn <- function() {
  list(gamma_inv = 11, beta = rep(0.1,ctp_data$n_days),
       ifr=0.008, phi=0.5,sigma_phi=0.5)
}
fit <- stan(file = 'code/Ohio_final.stan', init = init_fn, 
            data = OH_dat, chains = 4, iter = 10000,
            control = list(max_treedepth = 15,adapt_delta = 0.995))
            # data = OH_dat, chains = 4, iter = 20000,
            # control = list(max_treedepth = 15,adapt_delta = 0.99))
saveRDS(fit, "fits/OH_fit.rds")

sims <- rstan::extract(fit)
save(sims,file = "fits/OH_sims.RData")

####### fit truncated normal dist to posterior gamma_inv
tn_mll <- function(mu,sigma){
  mll <- -mean(dnorm(sims$gamma_inv,mean=mu,sd=sigma,log=TRUE)) +
    log(pnorm(11.5,mean=mu,sd=sigma)-pnorm(5.5, mean=mu,sd=sigma))
  return(mll)  
}

tn_fit <- mle(tn_mll,start = c(mean(sims$gamma_inv),
                               sd(sims$gamma_inv)))

gamma_inv <- list(mean = coef(tn_fit)[1], 
                     sd = coef(tn_fit)[2])
save(gamma_inv,file = "data/gamma_inv.RData")

####### fit normal dist to posterior phi
phi <- list(mean = mean(sims$phi), 
               sd = sd(sims$phi))
save(phi,file = "data/phi.RData")
