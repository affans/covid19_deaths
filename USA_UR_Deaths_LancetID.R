## MCMC analysis for estimating the true COVID-19 deaths 
## 
## Raw data files
##   Downloaded using download_data.R
##   incidence*, cumulative* => raw data files thats sent to Seyed to make CFR
##
## Inputs:
##   1) meanCFR.csv = a daily CFR file Seyed creates using daily # tests, incidence, and reported deaths
##   2) Death.csv = a file that contains the reported deaths 
##   in both these files, the columns represent each state. 
##
## Outputs:
##   - Creates posterior_y, posterior_z data files for y, z (true, observed) for each state in the /data folder. 
##   - Use the file USA_UR_Deaths_Plots.R to plot these output files in a panel.
##
## How to run: 
##   - Set the arg_st variable to a valid state acronym (ie from `validstates`)
##   - caveat: some states require different prior values. So make sure the right lines are commented and not
##   - automation 
##   --- serial_bsh_script will run all the states serially ~1 hour. 
##   --- parallel_bsh_script will submit jobs to cluster ~ 5 - 8 minutes for all states
##   --- remember to sbatch the parallel script. running it as a bash script will error out.
##   --- note the caveat. After all states have finished running, manually run the 7 states that require a different prior and overwrite the results
## 
## What happens after?
##   - the posterior files are sent to Seyed who uses a matlab script to make the time "uniform"
##   - that is.. pad the posterior data tables with 0 columns until there are 1 to 141 columns in all files
##   - I have code for this too (see at the end)... So might be worth it to just uncomment the code and use the parallel/script script without sending to Seyed
##   - Once Seyed (or my script) pads the data tables, 
##   - use my Julia script to calculate national level. 

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

# this file modified on june 21 to work on state level data. 
# uses command line args to get the state variable
# if running script locally, just comment this out and pass in the state manually
args = commandArgs(trailingOnly=TRUE)

## lets load the data/constants from all the states 
all_death_data = fread("Death.csv")
all_cfr_data = fread("meanCFR.csv")

## valid states and rename the columns
validstates = c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
names(all_death_data) <- validstates
names(all_cfr_data) <- validstates

arg_st = validstates[as.numeric(args[1])]
arg_st = "FL"

# ID ## ask mean = 0.2
# ME
# ND
# SD
# VT
# WV
# OR

writedata = F
print(arg_st)
popsizes = c(4887871, 737438, 7171646, 3013825, 39557045, 5695564, 3572665, 967171, 702455, 21299325, 10519475, 1420491, 1754208, 12741080, 6691878, 3156145, 2911505, 4468402, 4659978, 1338404, 6042718, 6902149, 9995915, 5611179, 2986530, 6126452, 1062305, 1929268, 3034392, 1356458, 8908520, 2095428, 19542209, 10383620, 760077, 11689442, 3943079, 4190713, 12807060, 1057315, 5084127, 882235, 6770010, 28701845, 3161105, 626299, 8517685, 7535591, 1805832, 5813568, 577737)

## ORIGINAL LANCET ID SUBMISSION DATA FOR ALL OF USA
#usa_death_data = c(1, 0, 5, 3, 2, 1, 3, 4, 3, 4, 4, 8, 3, 7, 10, 15, 22, 26, 50, 68, 70, 65, 135, 180, 268, 303, 354, 496, 644, 497, 697, 1124, 1191, 1175, 1256, 1537, 1401, 1496, 2219, 2156, 2101, 2226, 2013, 1715, 1714, 2553, 2618, 2176, 2528, 1867, 1561, 1939, 3022, 2358, 2340, 1957, 2065, 1156, 1383, 2470, 2390, 2201, 1897, 1691, 1153, 1324, 2350, 2528, 2129, 1687, 1422, 750, 1060, 1871, 1822, 1753, 1602)
#usa_death_data = c(1, 0, 5, 3, 2, 1, 3, 4, 3, 4, 4, 8, 3, 7, 10, 15, 22, 26, 50, 68, 70, 65, 135, 180, 268, 303, 354, 496, 644, 497, 697, 1124, 1191, 1175, 1256, 1537, 1401, 1496, 2219, 2156, 2101, 2226, 2013, 1715, 1714, 2553, 2618, 2176, 2528, 1867, 1561, 1939, 3022, 2358, 2340, 1957, 2065, 1156, 1383, 2470, 2390, 2201, 1897, 1691, 1153, 1324, 2350, 2528, 2129, 1687, 1422, 750, 1060, 1871, 1822, 1753, 1602, 1218, 865, 1003, 1552, 1403, 1418, 1293)
#chi_death_data = c(0, 3, 11, 0, 9, 15, 15, 25, 25, 26, 38, 43, 46, 45, 57, 65, 66, 72, 73, 86, 89, 97, 108, 97, 97, 97, 143, 142, 105, 98, 139, 112, 118, 109, 98, 150, 70, 52, 29, 44, 47, 35, 42, 32, 37, 31, 30, 28, 27, 23, 17, 22, 11, 7, 15, 9, 13, 9, 16, 8, 4, 6, 6, 9, 7, 4, 6, 5, 3, 5, 2, 3, 1, 6, 10, 4)
#sk_death_data = c(1, 0, 1, 3, 2, 3, 2, 1, 0, 4, 1, 4, 6, 4, 3, 7, 2, 6, 1, 3, 6, 6, 1, 5, 3, 0, 6, 3, 7, 3, 8, 2, 7, 9, 6, 5, 8, 5, 8, 6, 4, 3, 4, 5, 3, 6, 3, 6, 8, 4, 4, 3, 3, 3, 5, 3, 4, 1, 2, 2, 2, 1, 1, 2, 0, 0, 2, 1, 1, 2, 1, 1, 2, 0, 2, 2, 1, 1, 0, 0, 0, 0, 2, 1, 1, 0)
#cfr_means = c(0.007189492, 0.006474445, 0.029590088, 0.03627556, 0.035262704, 0.028662143, 0.025290245, 0.023706567, 0.022075661, 0.019712234, 0.01630829, 0.015726791, 0.013317298, 0.011742014, 0.010796488, 0.010263982, 0.010392547, 0.00950655, 0.00921385, 0.008603148, 0.007678606, 0.007128255, 0.006772992, 0.006770771, 0.007167756, 0.007298304, 0.007128662, 0.007517678, 0.007825363, 0.007567504, 0.007653341, 0.007987095, 0.00810327, 0.008047587, 0.007904116, 0.007732384, 0.007689936, 0.007534368, 0.007573694, 0.007630293, 0.007561125, 0.007515714, 0.00750592, 0.007359943, 0.007316432, 0.007439684, 0.007584599, 0.0076522, 0.007755632, 0.007790195, 0.007751095, 0.007855259, 0.008150941, 0.008333219, 0.00853811, 0.008644491, 0.008721969, 0.00874491, 0.008849346, 0.009125988, 0.00934131, 0.009511871, 0.009695766, 0.009805003, 0.00981613, 0.009893181, 0.010142782, 0.010391629, 0.010571942, 0.010689025, 0.010774242, 0.010884333, 0.01087724, 0.011044553, 0.011270986, 0.01133691, 0.011461102, 0.011347715, 0.011334453, 0.011456548, 0.01153147, 0.01168622, 0.011724216, 0.011782563)
#cfr_data = cfr_means
#usethisdata = usa_death_data

rowstoadd <<- 0 ## how many "0"s to add at the end to the posterior for uniformity 
## otherwise every state has a different start date

## define basic functions
get_cfr_death_data <- function(st){
  ## get the data data, and remove the zeros
  cfd = all_cfr_data[, get(st)]
  ddd = all_death_data[, get(st)]
  firstnonzero = min( which ( ddd != 0 ))
  rowstoadd <<- firstnonzero - 1
  lastelement = length(ddd)
  death_rm = ddd[firstnonzero : lastelement]
  cfr_rm = cfd[firstnonzero : lastelement] 
  return(list(cfr=cfr_rm, death=death_rm))
}

cfr_data = get_cfr_death_data(arg_st)$cfr
usethisdata = get_cfr_death_data(arg_st)$death
 
# if (arg_st == "WA"){
#   usethisdata[which(usethisdata == 73)] = 9
# }

UR_N = length(usethisdata)
qplot(x=1:UR_N, y=usethisdata, geom="line")
qplot(x=1:UR_N, y=cfr_data, geom="line")

skewnormal <- nimbleFunction(
  run = function(x=double(0), xi=double(0, default=1), omega=double(0, default=1), 
                 alpha=double(0, default=0), tau=double(0, default=0), log=double(0, default=0)){
    returnType(double(0))
    #x = x %% 91

    z <- (x-xi)/omega
    logN <- (-log(sqrt(2*pi)) -log(omega) - z^2/2)
    logS <- pnorm((alpha*z), log.p=TRUE)
    logPDF <- as.numeric(logN + logS - pnorm(tau, log.p=TRUE))
    #logPDF <- replace(logPDF, abs(x) == Inf, -Inf)
    #logPDF <- replace(logPDF, omega <= 0, NaN)
    if(log) return(logPDF) else return(exp(logPDF))
  }
)


# define the model
UR_code=nimbleCode({
  for(i in 1:n){
    #pi[i] <- ilogit(b0)
    #pi[i] <- ilogit(skewnormal(x=i, locA, scaleA, shapeA) * u[i])
    #pi[i] <- ilogit(skewnormal(x=i, locA, scaleA, shapeA) * 0.0178) ## okay results
    pi[i] <- ilogit(b0 + b1*u[i])
    
    # mod for two waves
    #lambda[i] <- exp(p + log(skewnormal(x=i, loc, scale, shape)) + log(K))
    lambda[i] <- exp(p + log(skewnormal(x=i, loc, scale, shape)) + log(K))
    
    z[i] ~ dpois(lambda[i]*pi[i])
    #z[i] ~ dbin(prob = 1, size = y[i])
  }
  
  # values for all states except a few selected
  asd ~ dnorm(mean=0.052348, sd=(0.62348)/10) # hyperparameter
  p ~ dnorm(mean=-12.19369, sd=abs(asd))
  locsd ~ dnorm(mean=0.6688287, sd=0.6688287/10)
  loc ~ dnorm(36.74267, sd=locsd)

  #scalesd ~ dnorm(mean=0.1527184, sd=0.1527184/10)
  scalesd ~ dnorm(mean=0.1627184, sd=0.1527184/10)
  scale ~ T(dnorm(9.238223, scalesd),0.001,)

  #shapesd ~ dnorm(mean=0.04587738, sd=0.04587738/10)
  shapesd ~ dnorm(mean=0.8587738, sd=0.04587738/10)
  shape ~ dnorm(0.006196017, sd=shapesd) #1] -0.1743705

  # ID ## ask mean = 0.2
  # ME
  # ND
  # SD
  # VT
  # WV
  # OR
  # asd ~ dnorm(mean=0.2, sd=(0.62348)/10) # hyperparameter
  # p ~ dnorm(mean=-12.19369, sd=abs(asd))
  # 
  # locsd ~ dnorm(mean=0.6688287, sd=0.6688287/10)
  # loc ~ dnorm(36.74267, sd=locsd)
  # 
  # scalesd ~ dnorm(mean=0.1627184, sd=0.1527184/10)
  # scale ~ T(dnorm(9.238223, scalesd),0.001,)
  # 
  # shapesd ~ dnorm(mean=0.8587738, sd=0.04587738/10)
  # shape ~ dnorm(0.006196017, sd=shapesd) #1] -0.1743705
  
  # # 
  # WA specifically.
  # asd ~ dnorm(mean=0.20, sd=(0.62348)/10) # hyperparameter
  # p ~ dnorm(mean=-12.19369, sd=abs(asd))
  # 
  # locsd ~ dnorm(mean=0.6688287, sd=0.6688287/10)
  # loc ~ dnorm(36.74267, sd=locsd)
  # 
  # scalesd ~ dnorm(mean=0.167184, sd=0.1527184/10)
  # scale ~ T(dnorm(9.238223, scalesd),0.001,)
  # 
  # shapesd ~ dnorm(mean=0.8587738, sd=0.04587738/10)
  # shape ~ dnorm(0.006196017, sd=shapesd) #1] -0.1743705

  # 
  b0 ~ dnorm(0,sd=0.6)  
  b1 ~ dnorm(0,sd=0.6)  
  #locA ~ dnorm(0, sd=10)
  #scaleA ~ T(dnorm(0, sd=10), 0.001,)
  #shapeA ~ dnorm(0, sd=10)
})


# Set up data for NIMBLE.
popsize_state = popsizes[which(validstates==arg_st)]
UR_constants=list(n=UR_N, K=popsize_state, u=cfr_data)
UR_data=list(z=usethisdata)

# Set initial values.
UR_inits1=list(asd=0, locsd=0.6647207, shapesd=0.04561767, scalesd=0.1527172, b0 = 0, b1 = 0, p = 2, loc=0, shape=1, scale=0.5 )#gval=skewnormal(x=1:UR_N))
UR_inits2=list(asd=0.2, locsd=0.6647207, shapesd=0.04561767, scalesd=0.1527172, b0 = 0, b1 = 0, p = 2.5, loc=0, shape=1, scale=1) #gval=skewnormal(x=1:UR_N))
UR_inits3=list(asd=0.5, locsd=0.6647207, shapesd=0.04561767, scalesd=0.1527172, b0 = 0, b1 = 0, p = 2, loc=-1, shape=1, scale=0.1) #gval=skewnormal(x=1:UR_N))
UR_inits4=list(asd=0, locsd=0.6647207, shapesd=0.04561767, scalesd=0.1527172, b0 = 0, b1 = 0, p = 1.4, loc=1, shape=1, scale=1.5) #gval=skewnormal(x=1:UR_N))
UR_inits5=list(asd=0, locsd=0.6647207, shapesd=0.04561767, scalesd=0.1527172, b0 = 0, b1 = 0, p = 1.9, loc=1, shape=1, scale=2) #gval=skewnormal(x=1:UR_N))

UR_inits=list(chain1=UR_inits1, chain2=UR_inits2, chain3=UR_inits3, chain4=UR_inits4, chain5=UR_inits5)

# Build the model.
UR_model <- nimbleModel(UR_code, UR_constants, UR_data, UR_inits)
UR_compiled_model <- compileNimble(UR_model, resetFunctions = TRUE) #dont need to do this.

# Set up MCMC samplers.
#UR_mcmc_conf <- configureMCMC(UR_model, monitors=c('a0', 'a1', 'b0', 'b1', 'epsilon', 'pi','lambda'),useConjugacy = TRUE)
UR_mcmc_conf <- configureMCMC(UR_model, monitors = c("lambda","pi",  "p", "loc", "scale", "shape"), useConjugacy = TRUE, print=TRUE)

UR_mcmc <- buildMCMC(UR_mcmc_conf)
UR_compiled_mcmc <- compileNimble(UR_mcmc, project = UR_model, resetFunctions = TRUE)

# Run the model (a few hours).
UR_samples=runMCMC(UR_compiled_mcmc,inits=UR_inits,
                   nchains = 5, nburnin=10000,niter = 50000, samplesAsCodaMCMC = TRUE,thin=10,
                   summary = FALSE, WAIC = FALSE,setSeed=c(978, 977, 944, 951, 172)) 

#  
df = data.table(do.call('rbind', UR_samples))

# Compute posterior quantities.
lstr = unlist(lapply(seq(1, UR_N), function(x){paste0("lambda[", x, "]")}))
posterior_lambda=df[, ..lstr]

lstr = unlist(lapply(seq(1, UR_N), function(x){paste0("pi[", x, "]")}))
posterior_pi=df[, ..lstr]

# Simulate z.
posterior_z=t(apply(posterior_pi*posterior_lambda,1,function(x)rpois(UR_N,x)))
mns1 = apply(posterior_z, 2, mean)
#posterior_lmse=apply(posterior_z,1,function(x) log(mean((x-usa_death_data)^2)))

# simulate y (true count)
# this formula comes from the TB under reporting paper
posterior_y=t(apply(posterior_lambda*(1-posterior_pi),1,function(x)rpois(UR_N,x)+mns1))
mns2 = apply(posterior_y, 2, mean)

# get cumulative numbers to inspect
c_data = sum(usethisdata)
c_mns1 = sum(mns1)
c_mns2 = sum(mns2)
c_Str = qq("data: @{c_data}, red: @{c_mns1}, blue: @{c_mns2}")
print(c_Str)
gg = ggplot()
gg = gg + geom_line(aes(x=1:UR_N, y=usethisdata))
gg = gg + geom_line(aes(x=1:UR_N, y=mns1), color="red")
gg = gg + geom_line(aes(x=1:UR_N, y=mns2), color="blue")
gg = gg + annotate("text", -Inf, Inf, label = c_Str, hjust = 0, vjust = 1)
gg


posterior_y = data.table(posterior_y)
posterior_z = data.table(posterior_z)

means_per_day = apply(posterior_y, 2, mean)
sums_per_sims = apply(posterior_y, 1, sum)

quantile(sums_per_sims, probs = c(0.025, 0.975))

## added this code on 24th june. 
## before this the posterior files have to be sent to Seyed
## he would pad the files so all files have 141 columns and time is uniform
## better to do it on my end so I can just run the Julia script after.

#mycols=lapply(1:rowstoadd, function(x) rep(0, 20000)) %>% as.data.frame
#names(mycols) = paste0("W", 1:rowstoadd)
#nys = cbind(mycols, nys)

#print("writing files \n")
if (writedata){
  fwrite(posterior_y, qq("/data/lancetid_actualdeaths_covid19/st_@{arg_st}_posterior_y.dat"))
  fwrite(posterior_z, qq("/data/lancetid_actualdeaths_covid19/st_@{arg_st}_posterior_z.dat"))
  print("writing data\n")
}




