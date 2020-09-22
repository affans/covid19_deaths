## MCMC analysis for estimating the true COVID-19 deaths 
## Affan Shoukat 2020
## 
## Raw data files are downloaded using download_data.R
##   incidence*, cumulative* => sent to Seyed to process underreporeted cases and testing data
##
##   returns three files to be used as input in this model
##   1) meanINC.csv = Incidence (including cases that may be underreporting)
##   2) Death.csv = contains the reported deaths
##   3) Test.csv =  contains the total number of tests administered
##   4) hosp_data = hospital prevalence (Seyed does not need this file)
##   in all these files, the columns represent each state. 
##
## Outputs:
##   - Creates posterior_y, posterior_z data files for y, z (true, observed) for each state in the /data folder. 
##.  - Creates dataplot_{state}.pdf files to visualize the input raw data
##   - Use the file create_plots.R to plot these output files in a panel.
##.  - Use the julia script "sample_national_level.jl" to get process data at the national level

## How to run: 
##   - Set the arg_st variable to a valid state acronym (ie from `validstates`)
##   - automation if cluster with slurm
##   --- parallel_bsh_script will submit jobs to cluster ~ 5 - 8 minutes for all states
##   --- remember to sbatch the parallel script. running it as a bash script will error out.
## 

rm(list = ls())
library(nimble)
library(coda) # For manipulation of MCMC results.
library(mgcv)
library(ggplot2)
library(ggpubr)
#library(tidyverse) # not on the cluster
library(ggmcmc)
#library(ggpubr)
library(data.table)
library(reshape)
library(dplyr)
library(boot)
library(ggthemes)
library(gridExtra)
library(GetoptLong)
library(qgam)  # for fitting
library(zoo) 
library(stringr)

# load data, constants, and functions ---------------------------------------------------------------
# get state metadata
city_metadata = fread("city_state_metadata.csv", colClasses = 'character')
city_metadata$statepop = as.numeric(gsub(",", "", city_metadata$statepop)) # turn to numeric

# state order in meanCFR,Death, and other files from Seyed 
validstates = c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
          
# load all the data
all_death_data = fread("Death.csv")
all_test_data = fread("Test.csv")
# for test data: flip as its in descending order and remove the date column
# all_test_data = all_test_data %>% mutate(sortk = seq(nrow(all_test_data), 1)) %>% arrange(sortk)
# all_test_data[, sortk := NULL]
ur_case_data = fread("./meanINC.csv")
# ur_case_data = as.data.table(lapply(ur_case_data, cumsum)) # turn columns into cumulative

# get the hospital data... since this is raw data from covid tracking (and not processed by Seyed)
# we have to arrange by date column and only get the number of rows as same as Seyeds files.
# and also remove the date column 
all_hosp_data = fread("hosp_data.csv")
all_hosp_data = all_hosp_data %>%  mutate(date = as.Date(as.character(date), "%Y%m%d")) %>%
  arrange(date) %>% head(nrow(ur_case_data)) %>% select(-!!c("date"))

# check if the order of states match, although this is is important since we select by column name and not index
names(all_hosp_data) == validstates

# attach column headers to all data frames
names(all_death_data) <- validstates
names(all_test_data) <- validstates
names(ur_case_data) <- validstates

## define global functions
skewnormal <- nimbleFunction(
  run = function(x=double(0), xi=double(0, default=1), omega=double(0, default=1), 
                 alpha=double(0, default=0), tau=double(0, default=0), log=double(0, default=0)){
    returnType(double(0))
    z <- (x-xi)/omega
    logN <- (-log(sqrt(2*pi)) -log(omega) - z^2/2)
    logS <- pnorm((alpha*z), log.p=TRUE)
    logPDF <- as.numeric(logN + logS - pnorm(tau, log.p=TRUE))
    if(log) return(logPDF) else return(exp(logPDF))
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
get_state_data_vectors <- function(st, ma=F){
  ## get the data data, and remove the zeros
  ddd = all_death_data[, get(st)]
  urc = ur_case_data[, get(st)]
  tst = all_test_data[, get(st)]
  hos = all_hosp_data[, get(st)]
  
  firstnonzero = min( which ( ddd != 0 )) 
  lastelement = length(ddd)
  
  if(ma){
    death_rm = rollmean(ddd[(firstnonzero-2): lastelement], 3)
    urc_rm = rollmean(urc[(firstnonzero-2) : lastelement], 3)
    tst_rm = rollmean(tst[(firstnonzero-2) : lastelement], 3)
    hos_rm = rollmean(hos[(firstnonzero-2) : lastelement], 3)
  } else {
    death_rm = ddd[firstnonzero : lastelement]
    urc_rm = urc[firstnonzero : lastelement]
    tst_rm = tst[firstnonzero : lastelement]
    hos_rm = hos[firstnonzero : lastelement]
  }
  # TX and GA have two NAs at start of data, replace with their 3rd index
  if(st == "TX" || st == "GA")
    tst_rm[is.na(tst_rm)] = tst_rm[3]
  return(list(death=death_rm, urcases=urc_rm, tstdata=tst_rm, hospdata=hos_rm))
}

fill_hospital_data <- function(hospdata, deathdata){
  # get the first non NA + 7 days of data, also get the death data
  na_indx = which(is.na(hospdata))
  if(length(na_indx) == 0){
    return(hospdata)
  }
  firstnonna = max(na_indx) + 1 
  f_hospdata =  hospdata[firstnonna:(firstnonna+7)]
  f_deathdata = deathdata[firstnonna:(firstnonna+7)]
  f_urcases = urcases[firstnonna:(firstnonna+7)]
  # take the ratio of death
  #f_ratio = f_deathdata / f_hospdata
  f_ratio =  f_hospdata / f_urcases

  # get trhe mean of the ratios
  f_mean = mean(f_ratio)
  if (f_mean < 0.0005) 
    f_mean = 1 
  
  # use the death data to estimate the hospitalized for the NA part 
  #estimated_hosp = deathdata[na_indx] / f_mean 
  estimated_hosp = urcases[na_indx] * f_mean
  
  
  # the slight problem here is if the deaths = 0, this makes hospitalized = 0 
  # so solution to this would be to take a moving average (which still may not solve the problem)
  # or do regression on the data 
  # either way, append to the hospital data first 
  hospdata[na_indx] <- estimated_hosp
  hospdata <- rollapply(hospdata, width = length(na_indx), FUN = mean, align = "center", partial = TRUE)
  
  #simpler solution, for any 0.. find the first non-zero on the left of the data 
  # since the deathdata is always >0 at day 1, this should always work 
  # zero_indx = which(hospdata[na_indx] == 0)
  # for (i in zero_indx){
  #   hospdata[i] = hospdata[i - 1]
  # }
  
  plot(x=(firstnonna+1):length(hospdata), y=hospdata[(firstnonna+1):length(hospdata)], 
       xlab="time", ylab="prevalence hosp", xlim=c(1, length(hospdata)), ylim=c(0, max(hospdata)))
  points(x=1:firstnonna, y=hospdata[1:firstnonna], pch=16, cex=1.2, col="red")
  return(hospdata)
}

# helper function for transforming data between a and b
minmax <- function(dat, a, b){
  a + (dat - min(dat))*(b-a) / (max(dat) - min(dat))
}

# function to read excess deaths from CDC data 
read_cdc_excess_death_data <- function(){
  # cdc file "excess_deaths.csv" downloaded on September 9th, 2020
  excess <- fread("excess_deaths.csv", select = c("Week Ending Date", "Year", "State", "Type", "Outcome", "Exceeds Threshold", "Observed Number", "Average Expected Count", "Percent Excess Lower Estimate", "Percent Excess Higher Estimate"))
  names(excess) <- c("week", "year", "state", "type", "outcome", "thresh", "observed", "expected", "percent_lo", "percent_hi")
  excess$week = as.Date(excess$week)
  excess = excess %>% arrange(week)
  ef = excess %>% filter(year == 2020 & outcome == "All causes, excluding COVID-19" & 
                           thresh==T & week < "2020-07-30")
  
  pmeans = ef %>% group_by(state) %>% 
    summarise(percent_val = (sum(observed) - sum(expected))/sum(expected))
  pmeans = pmeans %>% left_join(city_metadata, by="state")
  
}
fmeans = read_cdc_excess_death_data()

get_normalized_deaths <- function(st){
  cumdeaths = sum(get_state_data_vectors(st, ma=F)$death)
  popsize = (city_metadata %>% filter(abbr == st))$statepop
  rt = cumdeaths / popsize * 100000 # normalize to population 
  rt = log(rt) # take the log 
  return(rt)
}

# old values for b0, informed using cdc but not really used in the model 
bzvals = list("AK" = 10.5, "AL" = 2.0, "AR" = 3.2, "AZ" = 1.5, "CA" = 2.4, "CO" = 1.6, "CT" = 1.5, "DC" = 1.6, "DE" = 1.3, 
              "FL" = 1.6, "HI" = 10.5, "IA" = 1.8, "ID" = 2.1, "IL" = 1.6, "IN" = 2.2, "KS" = 1.8, "KY" = 1.8, "LA" = 1.4,
              "MA" = 1.2, "MD" = 1.0, "ME" = 10.5,  "MI" = 1.0, "MN" = 1.5,  "MO" = 2.5, "MS" = 1.5, "MT" = 4.0, "NC" = 2.3, "ND" = 2.2, 
              "NE" = 2.0, "NH" = 1.3, "SD" = 10.5, "NJ" = 0.85, "NM" = 1.8, "NV" = 1.6, "NY" = 1.1, "OH" = 2.0, "OK" = 2.0, "OR" = 2.7, "PA" = 3.5,
              "RI" = 1.0, "SC" = 1.5, "TN" = 2.5, "UT" = 3.5, "VA" = 2.2, "VT" = 3.8, "WA" = 2.2, "WI" = 3.5,  "WV" = 2.3, "WY" = 4.4, "TX" = 2.0, "GA" = 1.9)
# convert the bzvals into a table for joining purposes 
l_tib <- data.table(bzvals %>% unlist(recursive = FALSE) %>% tibble::enframe())
names(l_tib) = c("abbr", "bzval")

# get normalized deaths for each state 
st_norm_deaths <- unlist(purrr::map(validstates, get_normalized_deaths))
st_norm_deaths <- data.table(states=as.factor(validstates), logdeaths=st_norm_deaths)
# join the old cdc values for comparison purposes 
st_norm_deaths <- st_norm_deaths %>% left_join(l_tib, by = c("states" = "abbr"))
#st_norm_deaths <- st_norm_deaths %>% mutate(meanval=mean(logdeaths))
st_norm_deaths <- st_norm_deaths %>% mutate(redblue_diff = abs(bzval - logdeaths)/2)
#seyed_dt <- seyed_dt %>% mutate(diffdis = dat - linept) %>% mutate(newpt = linept - diffdis)

# we need to adjust the logdeaths values to calibrate the model. 
# the model should be calibrated so that NY gives the reported underreporting ~30%
correction_factor = (1.88 + 1.1) # where 1.88 is the diff between NY bzval and CDC excess death.
st_norm_deaths <- st_norm_deaths %>% mutate(newbz = 2*correction_factor - logdeaths)

# have to manually adjust PA as the logdeaths produces a very nonrealistic number
# we pick 3.5 as that gets us to the maximum of excess deaths reported by CDC
st_norm_deaths <- st_norm_deaths %>% mutate(newbz=replace(newbz, states=="PA", 3.5))

#smm <- melt(st_norm_deaths, id.vars = "states", measure.vars = c("logdeaths", "newbz"))
gg <- ggplot(st_norm_deaths)
colors <- c("logdeaths" = "red", "cdc_estim" = "blue", "corrected" = "orange")
gg <- gg + geom_point(aes(x=states, y=logdeaths, color="logdeaths"), alpha=0.5)
gg <- gg + geom_line(aes(x=states, y=logdeaths, group=1, color="logdeaths"), alpha=0.2)
gg <- gg + geom_point(aes(x=states, y=bzval, color="cdc_estim"))
gg <- gg + geom_line(aes(x=states, y=bzval, group=1, color="cdc_estim"), linetype = "dashed", alpha = 0.8)
gg <- gg + geom_point(aes(x=states, y=newbz,  color="corrected"))
gg <- gg + geom_line(aes(x=states, y=newbz, group=1,  color="corrected"))
gg <- gg + scale_color_manual(values=colors)
gg <- gg + ylim(c(0, 6.0))
gg <- gg + geom_hline(aes(yintercept=correction_factor))
#ggsave(filename = "bzero_informed_values.pdf", plot = gg, height = 7, width=16)
gg

# model def and implementation ---------------------------------------------------------------

# define the model in NIMBLE
lmod_code = nimbleCode({
  for(i in 1:n){
    ## data generating process
    lambda[i] <- a0 + log(K) + a1*deaths[i]
    y[i] ~ dpois(lambda[i]) 
    ## out of true counts, binomial thinning
    pi[i] <- ilogit(b0 + b1*x1[i] + b2*x2[i] + b3*x3[i])
    z[i] ~ dbin(prob = pi[i], size = y[i])
  }
  a0 ~ dnorm(0, sd=log(K)) 
  a1 ~ T(dnorm(10, sd=1), 1,) 
  b0 ~ dnorm(bzval, sd = 0.1) 
  b1 ~ dnorm(0, sd = 10)
  b2 ~ dnorm(0, sd = 10)
  b3 ~ dnorm(0, sd = 10) 
})

# informative mean values for b0 in the MCMC model. 
# these values are estimated using CDC excess deaths.
# see function read_cdc_excess_death_data()

# get command line argument for state
args = commandArgs(trailingOnly=TRUE)
arg_st = validstates[as.numeric(args[1])]
#arg_st = "IA" 
WRITE_DATA = T # write data files at the end
print(qq("Working with state: @{arg_st}"))

bzvalue = st_norm_deaths %>% filter(states == arg_st) %>% select(newbz) %>% .$newbz
print(qq("bzval: @{bzvalue}"))
popsize = (city_metadata %>% filter(abbr == arg_st))$statepop
  
deathdata = get_state_data_vectors(arg_st, ma=F)$death
nobs = length(deathdata)
# ma = t for the next three
deathdata_ma = get_state_data_vectors(arg_st, ma=T)$death
urcases = get_state_data_vectors(arg_st, ma=T)$urcases
tstdata = get_state_data_vectors(arg_st, ma=T)$tstdata
hospdata = get_state_data_vectors(arg_st, ma=T)$hospdata
hospdata = fill_hospital_data(hospdata, deathdata) 

#fit the data (fit to rolling avg to smooth out, seyed's idea)
poly_tim = fitlmn(deathdata_ma, 0.5)$fitted.values
poly_tim[poly_tim < 0] = 1e-5 # for model stability

## plot the data to review later 
dat_df = data.table(time=1:nobs, death=as.numeric(deathdata), cases=urcases, tests=tstdata, hosp=hospdata)
dat_df_m = melt(dat_df, id.vars = "time", measure.vars = c("death", "cases", "tests", "hosp"))
gg = ggplot(dat_df_m) 
gg = gg + geom_line(aes(x=time, y=value))
gg = gg + facet_wrap(facets=vars(variable), nrow = 2, ncol = 2, scales="free") 
# labeller = as_labeller(c(A = "Currents (A)", V = "Voltage (V)") )
gg = gg + ylab(NULL) + xlab("Time")
gg
if (WRITE_DATA){
  ggsave(qq("/data/actualdeaths_covid19/dataplot_@{arg_st}.pdf"), plot = gg, width=6.54, height=3.96)  
}

# x1: proportion of people dead out of total infected (basically cfr, but with timing issue)
prop_death = deathdata_ma / urcases

# x2: interpret: estimated cases that were tested
risk_death = urcases / tstdata 

# x3: death divided by hosp - proportion of deaths that were hospitalized
# the ratio of death to hospitalized cases. informing to some degree the proportion of under reporting
risk_hosp = deathdata_ma / hospdata

# Set initial values.
inits1=list(a0=0, a1=1, b0 = 0, b1 = 7.5, b2=8.5, b3=1, y=deathdata)
inits2=list(a0=0, a1=1, b0 = 0, b1 = 8, b2=8, b3=1, y=deathdata)
inits3=list(a0=0, a1=1, b0 = 0, b1 = 8.5, b2=7.5, b3=1, y=deathdata)
inits=list(chain1=inits1, chain2=inits2, chain3=inits3)

# load the constants and data
lmod_constants = list(n = nobs, K=popsize, deaths=poly_tim, x1=risk_death, x2=prop_death, x3=risk_hosp, bzval=bzvalue)
lmod_data = list(z = deathdata)

#Build the model.
model <- nimbleModel(lmod_code, lmod_constants, lmod_data, inits)
compiled_model <- compileNimble(model, resetFunctions = TRUE)
mcmc_conf <- configureMCMC(compiled_model, monitors=c("lambda", "pi", "a0", "a1", "b0", "b1", "b2", "b3"), print=T)
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = model)

mcmc_samples = runMCMC(compiled_mcmc, inits=inits, nchains=3, nburnin=10000, niter=50000, thin=10,
                       samplesAsCodaMCMC = TRUE, summary=T, setSeed=c(1, 2, 3))

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

lstr = c("a0", "a1", "b0", "b1", "b2")
trace_plots = df[, ..lstr]
trace_plots = trace_plots %>% mutate(itr=rep(seq(1, 4000), 3)) %>% mutate(chain=rep(c(1, 2, 3), 1, each=4000))
# tpm = melt(trace_plots, measure.vars = lstr)
# ggs(tpm)
#tpm = ggs(mcmc_samples$samples)
#ggmcmc(tpm)

if (WRITE_DATA){
  #bzfs = str_replace_all(BZVAL,"[.]","-")
  bzfs = "00"
  posterior_y = data.table(posterior_y)
  posterior_z = data.table(posterior_z)
  fwrite(posterior_y, qq("/data/actualdeaths_covid19/st_@{arg_st}_@{bzfs}_posterior_y.dat"))
  fwrite(posterior_z, qq("/data/actualdeaths_covid19/st_@{arg_st}_@{bzfs}_posterior_z.dat"))
  fwrite(trace_plots, qq("/data/actualdeaths_covid19/st_@{arg_st}_@{bzfs}_traceplots.dat"))
  print(qq("writing data for @{arg_st}\n"))
}

# get cumulative numbers to inspect/plot
c_data = sum(deathdata)
c_mns1 = round(sum(mns1), 2)
c_mns2 = round(sum(mns2), 2)
c_pinc = round((c_mns1 - c_data)/c_data, 2)
c_finc = round((c_mns1 - sum(poly_tim))/sum(poly_tim), 2)

c_binc = round((c_mns1 - c_mns2)/c_mns2, 2)
cdc_pexcess = round(as.numeric(fmeans %>% filter(abbr == arg_st) %>% select("percent_val") ), 2)
xvals = length(deathdata)
c_Str = qq(" st: @{arg_st}, data: @{c_data}, true: @{c_mns1}, obs: @{c_mns2}, %inc (data): @{c_pinc}, %inc (blue): @{c_binc}, cdc: @{cdc_pexcess}")
print(c_Str)

plot(df$a0)
plot(df$a1)
plot(df$b0)
plot(df$b1)
plot(df$b2)
plot(df$b3)
plot(apply(posterior_pi, 2, mean))

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



# loopvals = c("0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.1", "1.2", "1.3", "1.4", "1.5")
# loopvals = c("1.6", "1.7", "1.8", "1.9", "2.0", "2.1", "2.2", "2.3", "2.4", "2.5")
# loopvals = c("0.1", "0.2", "0.3", "0.4", "2.6", "2.7", "2.8", "2.9", "3.0", "3.1", "3.2", "3.3", "3.4", "3.5")
# for(i in loopvals){
#  print(qq("running loop value: @{i}"))
#  mcmc_samples = process_mcmc(i)
# }

