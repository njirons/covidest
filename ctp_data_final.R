#####################
##                 ##
##  Get CTP Data   ##
##                 ##
#####################

# moving average function
ma <- function(x,n){#n = number to average over on each tail
  L <- length(x)
  res = x
  for(k in 1:L){
    i1 <- max(c(1,k-n)); i2 <- min(c(L,k+n));
    res[k] = mean(x[i1:i2])
  }
  return(res)
}

DF <- read_csv("https://covidtracking.com/data/download/all-states-history.csv")
DF <- DF %>% mutate(day = yday(ymd(date))+366*(year(date)==2021))

abbrev <- read_csv("data/state_abbrev.csv")
#states <- c(abbrev$abbreviation,"DC","PR")
states <- c(abbrev$abbreviation,"DC")
#names <- c(abbrev$state, "District of Columbia","Puerto Rico")
names <- c(abbrev$state, "District of Columbia")
pop <- read_csv("data/pop_data.csv")

start_date <- "2020-01-01"
#start_date <- "2020-03-01"
end_date <- Sys.Date()
#end_date <- "2020-09-15"
#end_date <- "2020-06-06"

ifr_min <- 0
ifr_max <- 0.03

#p <- 0.01  # cutoff fraction of population tested
p <- 0

len_per_temp <- 6 # length of periods for cases/tests
#len_per_temp <- 2

cutoff <- 40 # longest time to death

#add_days <- max(c(0,cutoff - (first_death - first_day)))
add_days <- 0
cut <- 0  # for predicting deaths

#time to death distribution (truncated negative binomial)
alpha <- 21 
beta <- 1.1 
denom <- pnbinom(cutoff, size = alpha, prob=1/(beta+1))
ttd <- dnbinom(0:cutoff,size=alpha, prob = 1/(beta+1))/denom

for(s in 1:length(states)){
  State <- states[s]
  name <- names[s]
  N = pop$Pop[pop$State == name] #population of state
  
  df <- DF %>% filter(state == State) %>%
    select(state,date,day,positive,positiveIncrease,
           totalTestResults, totalTestResultsIncrease,
           death,deathIncrease)
  df[is.na(df)] <- 0
  
  if(State == "CO"){
    df <- DF %>% filter(state == State, date >= "2020-03-11") %>%
      select(state,date,day,positive,positiveIncrease,
             totalTestResults, totalTestResultsIncrease,
             death,deathIncrease)
    df[is.na(df)] <- 0
  }
  if(State == "GA"){
    df <- DF %>% filter(state == State, date >= "2020-03-09") %>%
      select(state,date,day,positive,positiveIncrease,
             totalTestResults, totalTestResultsIncrease,
             death,deathIncrease)
    df[is.na(df)] <- 0
  }
  if(State == "KY"){
    df <- DF %>% filter(state == State) %>%
      select(state,date,day,positive,positiveIncrease,
             totalTestResults, totalTestResultsIncrease,
             death,deathIncrease)
    df[is.na(df)] <- 0
    
    problem_day <- which(df$date == "2020-11-06")
    replacement <- sum(c(df$totalTestResultsIncrease[problem_day],
                         df$totalTestResultsIncrease[problem_day-1]))
    df$totalTestResultsIncrease[c(problem_day-1,problem_day)] =
      c(0,replacement)
  }
  if(State == "MD"){
    df <- DF %>% filter(state == State, date >= "2020-03-18") %>%
      select(state,date,day,positive,positiveIncrease,
             totalTestResults, totalTestResultsIncrease,
             death,deathIncrease)
    df[is.na(df)] <- 0
  }
  if(State == "MA"){
    df <- DF %>% filter(state == State) %>%
      select(state,date,day,positive,positiveIncrease,
             totalTestResults, totalTestResultsIncrease,
             death,deathIncrease)
    df[is.na(df)] <- 0
    
    prob_case_day_1 <- which(df$date == "2020-06-01")
    prob_case_day_2 <- which(df$date == "2020-09-02")
    reapportionment <- sum(c(df$positiveIncrease[prob_case_day_1],
                             df$positiveIncrease[prob_case_day_2]))
    df$positiveIncrease[c(prob_case_day_1,prob_case_day_2)] = c(0,0)
    df$positiveIncrease[prob_case_day_2:nrow(df)] =
      df$positiveIncrease[prob_case_day_2:nrow(df)]*
      (1+reapportionment/sum(df$positiveIncrease[prob_case_day_2:nrow(df)]))
  }
  if(State == "NC"){
    df <- DF %>% filter(state == State) %>%
      select(state,date,day,positive,positiveIncrease,
             totalTestResults, totalTestResultsIncrease,
             death,deathIncrease)
    df[is.na(df)] <- 0
    
    problem_day <- which(df$date == "2020-08-12")
    reapportionment <- df$totalTestResultsIncrease[problem_day]
    df$totalTestResultsIncrease[problem_day] = 0
    df$totalTestResultsIncrease[problem_day:nrow(df)] =
      df$totalTestResultsIncrease[problem_day:nrow(df)]*
      (1+reapportionment/sum(df$totalTestResultsIncrease[problem_day:nrow(df)]))
  }
  if(State == "WA"){
    df <- DF %>% filter(state == State, date >= "2020-02-26") %>%
      select(state,date,day,positive,positiveIncrease,
             totalTestResults, totalTestResultsIncrease,
             death,deathIncrease)
    df[is.na(df)] <- 0
  }
  
  ### clean data
  
  cum_deaths <- c(df$death,0)
  #cum_cases  <- c(df$positive,0)
  #cum_tests  <- c(df$totalTestResults,0)
  for (i in 2:length(cum_deaths)) {
    cum_deaths[i] <- min(cum_deaths[1:i], na.rm = T)
    #cum_cases[i]  <- min(cum_cases[1:i], na.rm = T)
    #cum_tests[i] <- min(cum_tests[1:i], na.rm = T)
  }
  df$deathIncrease <- rev(diff(rev(cum_deaths)))
  #df$positiveIncrease <- rev(diff(rev(cum_cases)))
  #df$totalTestResultsIncrease <- rev(diff(rev(cum_tests)))
  
  if(State =="NY"){
    problem_day <- which(df$date == "2020-05-07")
    extra_deaths <- df$deathIncrease[problem_day]
    df$deathIncrease[problem_day] <- mean(df$deathIncrease[c(problem_day-1,problem_day+1)])
    extra_deaths <- extra_deaths - df$deathIncrease[problem_day]
    df$deathIncrease[problem_day:nrow(df)] <- df$deathIncrease[problem_day:nrow(df)]*
      (1+extra_deaths/sum(df$deathIncrease[problem_day:nrow(df)]))
    df$deathIncrease <- round(df$deathIncrease)
  }
  
  for(i in nrow(df):1){
    if(df$positiveIncrease[i] < 0){
      cc_temp <- df$positiveIncrease[i]
      df$positiveIncrease[i] <- 0
      df$positiveIncrease[i:nrow(df)] = df$positiveIncrease[i:nrow(df)]*
        (1+cc_temp/sum(df$positiveIncrease[i:nrow(df)]))
    }
    if(df$totalTestResultsIncrease[i] < 0){
      test_temp <- df$totalTestResultsIncrease[i]
      df$totalTestResultsIncrease[i] <- 0
      df$totalTestResultsIncrease[i:nrow(df)] =
        df$totalTestResultsIncrease[i:nrow(df)]*
        (1+test_temp/sum(df$totalTestResultsIncrease[i:nrow(df)]))
    }
  }
  
  df$positiveIncrease <- ma(df$positiveIncrease,0)
  df$totalTestResultsIncrease <- ma(df$totalTestResultsIncrease,0)
  
  df$positive <- rev(cumsum(rev(df$positiveIncrease)))
  df$totalTestResults <- rev(cumsum(rev(df$totalTestResultsIncrease)))
  
  df <- df %>% filter(date >= start_date, date <= end_date)
  
  ###
  first_death <- min(df$day[df$deathIncrease > 0])
  first_day <- min(df$day)
  last_day <- max(df$day)
  new_first_day <- first_day - add_days
  
  if(add_days > 0){d_days <- c(new_first_day:(first_day-1),rev(df$day))}
  if(add_days==0){d_days <- rev(df$day)}
  
  ### for predicting deaths
  new_death_day <- last_day - cut # yday(end_date)-cut
  d_days <- d_days[d_days <= new_death_day]
  n_d_days <- length(d_days)
  d <- c(rep(0,add_days),rev(df$deathIncrease[df$day <= new_death_day]))
  
  ### cases/tests
  
  # test_days <- rev(df$day)
  # cum_tests <- rev(df$totalTestResults)
  # tests <- rev(df$totalTestResultsIncrease)
  # cum_cc <- rev(df$positive)
  # cc <- rev(df$positiveIncrease)
  
  tests_pre <- rev(df$totalTestResultsIncrease[df$totalTestResults/N <= p])
  cc_pre <- rev(df$positiveIncrease[df$totalTestResults/N <= p])
  cc_pre_total <- sum(cc_pre)
  tests_pre_total <- sum(tests_pre)
  
  cc_post <- rev(df$positiveIncrease[df$totalTestResults/N > p])
  tests_post <- rev(df$totalTestResultsIncrease[df$totalTestResults/N > p])
  test_days <- rev(df$day[df$totalTestResults/N > p])
  
  cum_tests <- rev(df$totalTestResults)
  cum_cc <- rev(df$positive)
  tests <- rev(df$totalTestResultsIncrease)
  cc <- rev(df$positiveIncrease)
  
  n_test_days <- length(test_days)
  
  new_first_day <- min(c(min(d_days),min(test_days)))
  shifted_d_days <- d_days - new_first_day + 1
  shifted_test_days <- test_days - new_first_day + 1
  
  n_days <- max(c(max(shifted_d_days),max(shifted_test_days)))
  days <- new_first_day:last_day
  
  dates = as.Date(days,origin="2019-12-31")
  
  # segment case/test data into periods
  len_per <- len_per_temp  # length of periods in days
  cum_tests_per <- c(0,0)
  
  #while(is.element(0, cum_tests_per[-1])){
  while(prod(cum_tests_per[-1] > 0)==0){
    len_per <- len_per+1
    cc_per <- split(cc_post, ceiling(seq_along(cc_post)/len_per)) # cases in each period
    tests_per <- split(tests_post, ceiling(seq_along(tests_post)/len_per)) # tests in each period
    if(p>0){
      cum_tests_per <- c(0,tests_pre_total, 
                        unlist(lapply(tests_per,sum),use.names = FALSE))
      cum_cc_per <- c(cc_pre_total,unlist(lapply(cc_per,sum),use.names = FALSE))
      per_ends <- c(0,shifted_test_days[1]-1,
                    shifted_test_days[cumsum(unlist(lapply(cc_per,length)))])  
    }else{
      cum_cc_per <- unlist(lapply(cc_per,sum),use.names = FALSE) # sum of cases in each period
      cum_tests_per <- c(tests_pre_total, unlist(lapply(tests_per,sum),use.names = FALSE)) # sum of tests in each period
      per_ends <- c(shifted_test_days[1]-1,
                    shifted_test_days[cumsum(unlist(lapply(cc_per,length)))]) # (shifted) end date of each period
    }
    n_per <- length(per_ends)-1
  }
  
  ctp_data <- list(state = State, name=name, N=N,
                   cutoff = cutoff, cut = cut, 
                   cc = cc, cum_cc = cum_cc,
                   tests = tests, cum_tests = cum_tests,
                   d = d, days = days,
                   n_days = n_days, d_days = d_days,
                   n_d_days = n_d_days, shifted_d_days = shifted_d_days,
                   test_days = test_days, n_test_days=n_test_days,
                   shifted_test_days=shifted_test_days, ttd = ttd,
                   start_date = start_date,end_date = end_date, 
                   dates = dates, len_per = len_per, n_per = n_per, 
                   cc_per = cc_per,cum_cc_per = cum_cc_per, 
                   tests_per = tests_per,cum_tests_per = cum_tests_per, 
                   per_ends = per_ends, tests_pre = tests_pre,
                   cc_pre = cc_pre, cc_pre_total = cc_pre_total,
                   tests_pre_total = tests_pre_total, p=p,
                   cc_post = cc_post, tests_post = tests_post,
                   ifr_min=ifr_min, ifr_max=ifr_max)
  save(ctp_data, file = paste0("data/",State,"_ctp.RData"))
}


