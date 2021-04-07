OH_fit <- readRDS("fits/OH_fit.rds")
IN_fit <- readRDS("fits/IN_fit.rds")
OH_sims <- rstan::extract(OH_fit)
IN_sims <- rstan::extract(IN_fit)

load("data/IN_ctp.RData")

### BHM fit
BHM_data <- list(n_samples_IN = length(IN_sims$ifr),
            n_samples_OH = length(OH_sims$ifr),
            phi_IN = IN_sims$phi,
            phi_OH = OH_sims$phi,
            phi_IN_sd = sd(IN_sims$phi),
            phi_OH_sd = sd(OH_sims$phi),             
            sigma_phi_IN = IN_sims$sigma_phi,
            sigma_phi_OH = OH_sims$sigma_phi,  
            sigma_phi_IN_sd = sd(IN_sims$sigma_phi),
            sigma_phi_OH_sd = sd(OH_sims$sigma_phi),  
            ifr_min = ctp_data$ifr_min, ifr_max = ctp_data$ifr_max,
            ifr_IN = IN_sims$ifr,
            ifr_OH = OH_sims$ifr,
            sigma_ifr_IN = sd(IN_sims$ifr),
            sigma_ifr_OH = sd(OH_sims$ifr)) 
init_fn <- function() {
  list(#mu_phi = 25, mu_phi_IN = 20, mu_phi_OH=30, sigma_phi=5,
       mu_ifr = 0.007, mu_ifr_IN = 0.007, mu_ifr_OH=0.007,
       sigma_ifr=0.002)
}
fit <- stan(file = 'code/BHM_final.stan', init = init_fn,
            data = BHM_data, chains = 4, iter = 10000,
            control = list(max_treedepth = 15,adapt_delta = 0.995))
            #data = BHM_data, chains = 4, iter = 20000,
            #control = list(max_treedepth = 15,adapt_delta = 0.99))
saveRDS(fit, "fits/BHM_fit.rds")

sims <- rstan::extract(fit)
save(sims, file = "fits/BHM_sims.RData")

# check posteriors
# stan_dens(fit,"mu_phi")
# stan_dens(fit,"tau_phi")
# stan_dens(fit,"mu_ifr")
# stan_dens(fit,"mu_ifr_IN")
# stan_dens(fit,"mu_ifr_OH")
# stan_dens(fit,"sigma_ifr")
# plot(density(sims$sigma_ifr), xlim=c(0,0.005))
# plot(density(1/sqrt(rexp(100000,rate = 0.0000008206947))),xlim=c(0,0.005))
# stan_dens(fit,"tau_ifr")

####### fit normal dist to posterior mu_ifr
# mll <- function(mu,sigma){
#   minusloglik <- -mean(dnorm(sims$mu_ifr,mean=mu,
#                              sd=sigma,log=TRUE))
#   return(minusloglik)  
# }
# 
# mll_fit <- mle(mll,start = c(mean(sims$mu_ifr),sd(sims$mu_ifr)),
#                control=list(trace=1))
# summary(mll_fit)

# plot(density(sims$mu_ifr))
# x <- seq(0,0.02,by=0.0001)
# #p <- dnorm(x,mean=coef(mll_fit)[1],sd=coef(mll_fit)[2])
# p_med <- dnorm(x,median(sims$mu_ifr),sd(sims$mu_ifr))
# p_mean <- dnorm(x,mean(sims$mu_ifr),sd(sims$mu_ifr))
# lines(x,p_mean,col = "red")
# lines(x,p_med,col = "blue")

#mu_ifr <- coef(mll_fit)
mu_ifr <- list(mean = mean(sims$mu_ifr),
               sd = sd(sims$mu_ifr))
save(mu_ifr,file = "data/mu_ifr.RData")

####### fit gamma dist to posterior tau_ifr
# mll <- function(alpha,beta){
#   minusloglik <- -mean(dgamma(sims$tau_ifr,shape=alpha,
#                              rate=beta,log=TRUE))
#   return(minusloglik)  
# }
# 
# mll_fit <- mle(mll,start = c((mean(sims$tau_ifr)/sd(sims$tau_ifr))^2,
#                    mean(sims$tau_ifr)/(sd(sims$tau_ifr)^2)), 
#                    lower = 0, control=list(trace=1))
# summary(mll_fit)

# plot(density(sims$tau_ifr))
# x <- seq(0,4*10^6,by=1000)
# #p <- dnorm(x,mean=coef(mll_fit)[1],sd=coef(mll_fit)[2])
# p <- dgamma(x, shape = (mean(sims$tau_ifr)/sd(sims$tau_ifr))^2,
#             rate = mean(sims$tau_ifr)/(sd(sims$tau_ifr)^2))
# #p <- dexp(x, rate = 1/mean(sims$tau_ifr))
# lines(x,p,col = "red")

tau_ifr <- list(shape = 
                  (mean(sims$tau_ifr)/sd(sims$tau_ifr))^2,
              rate = 
                mean(sims$tau_ifr)/(sd(sims$tau_ifr)^2))
save(tau_ifr,file = "data/tau_ifr.RData")

# ####### fit lognormal dist to posterior sigma_phi
# mll <- function(mu,sigma){
#   minusloglik <- -mean(dlnorm(sims$sigma_phi,meanlog=mu,
#                               sdlog=sigma,log=TRUE))
#   return(minusloglik)  
# }
# 
# mll_fit <- mle(mll,start = c(mean(log(sims$sigma_phi)),
#                              sd(log(sims$sigma_phi))))#, 
#                #lower = 0, control=list(trace=1))
# summary(mll_fit)
# 
# plot(density(sims$sigma_phi))
# x <- seq(0,40,by=0.01)
# p <- dlnorm(x,meanlog=coef(mll_fit)[1],sdlog=coef(mll_fit)[2])
# p <- dlnorm(x,mean(log(sims$sigma_phi)),sd(log(sims$sigma_phi)))
# lines(x,p,col = "red")
# 
# sigma_phi <- list(meanlog = 
#                     mean(log(sims$sigma_phi)),
#                   sdlog = 
#                     sd(log(sims$sigma_phi)))
# save(sigma_phi,file = "data/sigma_phi.RData")
# 
# ####### fit normal dist to posterior mu_phi
# mll <- function(mu,sigma){
#   minusloglik <- -mean(dnorm(sims$mu_phi,mean=mu,
#                              sd=sigma,log=TRUE))
#   return(minusloglik)  
# }
# mll_fit <- mle(mll,start = c(mean(sims$mu_phi),sd(sims$mu_phi)))
# summary(mll_fit)
# 
# plot(density(sims$mu_phi))
# x <- seq(0,50,by=0.01)
# #p <- dnorm(x,mean=coef(mll_fit)[1],sd=coef(mll_fit)[2])
# p_med <- dnorm(x,median(sims$mu_phi),sd(sims$mu_phi))
# p_mean <- dnorm(x,mean(sims$mu_phi),sd(sims$mu_phi))
# lines(x,p_mean,col = "red")
# lines(x,p_med,col = "blue")
# 
# #mu_ifr <- coef(mll_fit)
# mu_phi <- list(mean = mean(sims$mu_phi),
#                sd = sd(sims$mu_phi))
# save(mu_phi,file = "data/mu_phi.RData")
