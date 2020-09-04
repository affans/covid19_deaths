rm(list = ls())
library(nimble)
library(coda) # For manipulation of MCMC results.
library(mgcv)
library(ggplot2)
#library(ggpubr)
library(data.table)
library(reshape)
library(dplyr)
library(boot)
library(ggthemes)
library(gridExtra)
library(GetoptLong)
library(qgam)  # for fitting

# write data at the end
WRITE_DATA = F

# the order of the columns from Seyed's files 
validstates = c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")

# to get the population of each state, read the city_state_metadata file. 
city_metadata = fread("city_state_metadata.csv", colClasses = 'character')
city_metadata$statepop = as.numeric(gsub(",", "", city_metadata$statepop)) # turn to numeric

## lets load the data/constants from all the states 
all_death_data = fread("Death.csv")
all_cfr_data = fread("meanCFR.csv")
all_test_data = fread("Test.csv")
# flip the data as its in descending order 
# and remove the date column from test
all_test_data = all_test_data %>% mutate(sortk = seq(nrow(all_test_data), 1)) %>% arrange(sortk)
all_test_data[, sortk := NULL]
all_test_data[, date := NULL]
names(all_death_data) <- validstates
names(all_cfr_data) <- validstates
names(all_test_data) <- validstates

# problem 52 columns.. 
all_cases_data = fread("/data/incidence_positive_cases_lancetid.csv")
all_cases_data = all_cases_data[1:190, ]
#names(all_cases_data) <- validstates # extra state?

# estimated case data from Seyed
ur_case_data = fread("./meanINC.csv")
names(ur_case_data) <- validstates
ur_case_data = ur_case_data[1:190, ]
# turn columns into cumulative
#ur_case_data = as.data.table(lapply(ur_case_data, cumsum))

## define functions
skewnormal <- nimbleFunction(
  run = function(x=double(0), xi=double(0, default=1), omega=double(0, default=1), 
                 alpha=double(0, default=0), tau=double(0, default=0), log=double(0, default=0)){
    returnType(double(0))
    z <- (x-xi)/omega
    logN <- (-log(sqrt(2*pi)) -log(omega) - z^2/2)
    logS <- pnorm((alpha*z), log.p=TRUE)
    logPDF <- as.numeric(logN + logS - pnorm(tau, log.p=TRUE))
    #logPDF <- replace(logPDF, abs(x) == Inf, -Inf)
    #logPDF <- replace(logPDF, omega <= 0, NaN)
    if(log) return(logPDF) else return(exp(logPDF))
  }
)

sk_vec <- nimbleFunction(
  run = function(n=double(0), xi=double(0, default=1), omega=double(0, default=1), 
                  alpha=double(0, default=0), tau=double(0, default=0), log=double(0, default=0)){
    returnType(double(0))
    sknm <- rep(0, n)
    for (i in 1:n){
      sknm[i] <- skewnormal(i, xi, omega, alpha)
    }
    return(max(sknm))
  }
)

fitlmn <- function(dat, qu_val){
  # qgam paper: Statistical distribution fitting to the number of COVID-19 deaths in South Africa
  # follows the code samples from that paper -- can be optimized and cleaned up
  adma3 = data.table(t = 1:length(dat), De_New = dat)
  
  # take 100 rows for training, the rest as test. 
  De_data_test <- 100:nrow(adma3) 
  data_train <- adma3[-De_data_test, ] 
  data_test <- adma3[De_data_test, ]
  set.seed(2356)
  
  ## bs = "ps"
  tun <- tuneLearnFast(form=De_New~s(t,bs="ad"), err = 0.1, qu = qu_val, data = data_train) 
  fit1 <- qgam(De_New~s(t,bs="ad"), err = 0.1, qu = qu_val, lsig = tun$lsig, data = adma3) 
  
  #summary(fit1, se="ker")
  plot(adma3$De_New, xlab="Day, t", ylab="Reported deaths in SA", col="blue", type="b")
  lines(fit1$fit, col="red")
  return(fit1)
}

## function to get state specific vectors
get_cfr_death_data <- function(st, ma=F){
  ## get the data data, and remove the zeros
  cfd = all_cfr_data[, get(st)]
  ddd = all_death_data[, get(st)]
  urc = ur_case_data[, get(st)]
  
  firstnonzero = min( which ( ddd != 0 )) 
  lastelement = length(ddd)
  
  if(ma){
    death_rm = rollmean(ddd[(firstnonzero-2): lastelement], 3)
    cfr_rm = rollmean(cfd[(firstnonzero-2): lastelement], 3)
    urc_rm = rollmean(urc[(firstnonzero-2) : lastelement], 3)
  } else {
    death_rm = ddd[firstnonzero : lastelement]
    cfr_rm = cfd[firstnonzero : lastelement] 
    urc_rm = urc[firstnonzero : lastelement]
  }
  return(list(cfr=cfr_rm, death=death_rm, urcases=urc_rm))
}

minmax <- function(dat, a, b){
  a + (dat - min(dat))*(b-a) / (max(dat) - min(dat))
}


# model run ---------------------------------------------------------------

# define the model
lmod_code = nimbleCode({
  for(i in 1:n){
    ## data generating process
    lambda[i] <- a1*deaths[i]
    y[i] ~ dpois(lambda[i]) 
    
    ## out of true counts, binomial thinning
    #pi[i] <- ilogit(b0 + b1*(skewnormal(x=i, peak, s0, s1) * prop[i]) + log(K))
   
    pi[i] <- 1 -  (skewnormal(x=i, peak, s0, s1)/nax) * prop[i]
    z[i] ~ dbin(prob = pi[i], size = y[i])
  }
  
  a1 ~ T(dnorm(10, sd=10), 1,) 
  
  s0 ~ T(dnorm(log(K), sd = 0.01), 0, )
  s1 ~ dnorm(0, sd = 0.01)
  nax <- sk_vec(n, peak, s0, s1)
  #nax <- max(skewnormal(x=1:n, peak, s0, s1))

  #b0 ~ dnorm(-log(K), sd = 1)
  #b1 ~ dnorm(-log(K), sd = 1)
})

# data processing
#args = commandArgs(trailingOnly=TRUE)
#arg_st = validstates[as.numeric(args[1])]
arg_st = "NJ"
popsize = (city_metadata %>% filter(abbr == arg_st))$statepop


cfr_data = get_cfr_death_data(arg_st, ma=F)$cfr
deathdata_ma = get_cfr_death_data(arg_st, ma=F)$death
deathdata = get_cfr_death_data(arg_st, ma=F)$death
urcases = get_cfr_death_data(arg_st, ma=F)$urcases

nobs = length(deathdata)

#fit the data (fit to rolling avg to smooth out, seyed's idea)
poly_tim = fitlmn(deathdata_ma, 0.5)$fitted.values
poly_tim[poly_tim < 0] = 1e-5 # for model stability

length(deathdata) == length(deathdata_ma) 
length(cfr_data) == length(urcases)
length(deathdata) == length(cfr_data)
length(poly_tim) == length(deathdata)

peak= which(poly_tim == max(poly_tim[1:(nobs-45)]))
est_true = cfr_data[nobs] * urcases
est_prop = deathdata_ma / (est_true)
# to do, check fit for negative numbers
est_prop_fit = fitlmn(est_prop, 0.5)$fitted.values
est_prop_fit = minmax(est_prop_fit, 0, 0.75)
est_prop_fit[est_prop_fit < 0] = 1e-5

lmod_constants = list(n = nobs, K=popsize, deaths=poly_tim, peak=peak, mprop=max(est_prop), prop=est_prop_fit)
lmod_data = list(z = deathdata)

# Set initial values.
inits1=list(a1=10, s0 = log(popsize) - 1, s1 = 1, y=deathdata)
inits2=list(a1=10, s0 = log(popsize), s1 = 1,  y=deathdata)
inits3=list(a1=10, s0 = log(popsize), s1 = 1,  y=deathdata)
inits=list(chain1=inits1, chain2=inits2, chain3=inits3)

# Build the model.
model <- nimbleModel(lmod_code, lmod_constants, lmod_data, inits)
compiled_model <- compileNimble(model, resetFunctions = TRUE)
mcmc_conf <- configureMCMC(compiled_model, monitors=c("lambda", "pi", "a1", "s0", "s1"), print=T)
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = model)

mcmc_samples = runMCMC(compiled_mcmc, inits=inits, nchains=3, nburnin=10000, niter=50000, thin=10, 
                       samplesAsCodaMCMC = TRUE, summary=T, setSeed=c(1, 2, 3))

# Run the model (a few hours).
df = data.table(do.call('rbind', mcmc_samples$samples))

lstr = unlist(lapply(seq(1, nobs), function(x){paste0("lambda[", x, "]")}))
posterior_lambda=df[, ..lstr]
posterior_lambda
posterior_y=t(apply(posterior_lambda,1,function(x)rpois(nobs,x)))
mns1 = apply(posterior_y, 2, mean)

lstr = unlist(lapply(seq(1, nobs), function(x){paste0("pi[", x, "]")}))
posterior_pi=df[, ..lstr]

lst <-  list(posterior_pi, posterior_y)
ar1 <- array(unlist(lst), dim = c(dim(posterior_pi), length(lst)))
posterior_z <-  apply(aperm(ar1, c(3, 1, 2)), c(2,3), FUN = function(x) rbinom(n = 1, prob = x[1], size = x[2]))
mns2 = apply(posterior_z, 2, mean)

# get cumulative numbers to inspect
c_data = sum(deathdata)
c_mns1 = round(sum(mns1), 2)
c_mns2 = round(sum(mns2), 2)
c_pinc = round((c_mns1 - c_data)/c_data, 2)
c_finc = round((c_mns1 - sum(poly_tim))/sum(poly_tim), 2)

c_binc = round((c_mns1 - c_mns2)/c_mns2, 2)
xvals = length(deathdata)
c_Str = qq(" st: @{arg_st}, data: @{c_data}, true: @{c_mns1}, obs: @{c_mns2}, %inc (data): @{c_pinc}, %inc (blue): @{c_binc}, %inc (green): @{c_finc}")
print(c_Str)
gg = ggplot()
gg = gg + geom_point(aes(x=1:xvals, y=deathdata), color="#636363")
gg = gg + geom_line(aes(x=1:xvals, y=deathdata), color="#636363")
gg = gg + geom_line(aes(x=1:xvals, y=mns1), color="#e41a1c", size=1.2)
#gg = gg + ylim(0, 500)
gg = gg + geom_line(aes(x=1:xvals, y=mns2), color="#377eb8", size=1.2)
gg = gg + geom_line(aes(x=1:xvals, y=poly_tim), color="#4daf4a", size=1.2, alpha=0.8)
gg = gg + annotate("text", -Inf, Inf, label = c_Str, hjust = 0, vjust = 1)
gg = gg + xlab("time") + ylab("deaths")
gg
#ggsave(qq("/data/actualdeaths_covid19/st_@{arg_st}_final.pdf"))
# 

if (WRITE_DATA){
  posterior_y = data.table(posterior_y)
  posterior_z = data.table(posterior_z)
  fwrite(posterior_y, qq("/data/actualdeaths_covid19/st_@{arg_st}_posterior_y.dat"))
  fwrite(posterior_z, qq("/data/actualdeaths_covid19/st_@{arg_st}_posterior_z.dat"))
  print(qq("writing data for @{arg_st}\n"))
}


#b0 = mean(df$b0)
#plot(ilogit(b0 + b1(skewnormal(x=1:nobs, peak, mean(df$s0), mean(df$s1)) * est_prop) + log(popsize)))

#plot(ilogit(mean(df$b0) + mean(df$b1)*(skewnormal(x=1:nobs, peak, mean(df$s0), mean(df$s1)) * est_prop) + log(popsize)))


skn = skewnormal(x=1:nobs, peak, mean(df$s0), mean(df$s1))
plot(skn)
skn = (skn - min(skn) )/(max(skn) - min(skn))
sknp =  skn * est_prop_fit
plot(1 - sknp)

# 
# plot(est_prop)
# 
# fitlmn(sknp)
# 
# adma3 = data.table(t = 1:length(est_prop), De_New = est_prop)
# 
# # take 100 rows for training, the rest as test. 
# De_data_test <- 100:nrow(adma3) 
# data_train <- adma3[-De_data_test, ] 
# data_test <- adma3[De_data_test, ]
# set.seed(2356)
# 
# ## bs = "ps"
# tun <- tuneLearnFast(form=De_New~s(t,bs="ad"), err = 0.02, qu = 0.5, data = data_train) 
# fit1 <- qgam(De_New~s(t,bs="ad"), err = 0.02, qu = 0.5, lsig = tun$lsig, data = adma3) 
# 
# #summary(fit1, se="ker")
# plot(adma3$De_New, xlab="Day, t", ylab="Reported deaths in SA", col="blue", type="b")
# lines(fit1$fit, col="red")
# 
# plot(fit1$fit)
# 
# 
# plot(1 - skn*fit1$fit)
# plot uncertainty
# pym = melt(posterior_pi, id.vars = NULL, measure.vars = NULL )
# ggplot(pym) + geom_boxplot(aes(x=variable, y=value))
