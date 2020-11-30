## MCMC analysis for estimating the true COVID-19 deaths 
## Downloaded files here are used by Seyed/Ellie to generate CFR and UR incidence. 
## The file hosp_data.csv is not sent to Seyed, but still used in the Bayesian model

rm(list=ls())
library(data.table)
library(tidyverse)
library(GetoptLong)

ANALYSIS_DATE = "2020-11-20"
# science submission: data up to october 9

make_incidence <- function(cum_data){
  len = length(cum_data)
  incidence = numeric(len)
  for (i in 1:(len-1)) {
    incidence[i] = cum_data[i] - cum_data[i+1]
  }
  incidence[len] = 0
  return(incidence)
}
raw_data_dl <- fread('https://covidtracking.com/api/v1/states/daily.csv')

#mydat <- raw_data_dl[, c(1:5, 8, 19, 26, 27)] 
mydat <-raw_data_dl[, c("date", "state", "positive", "negative", "hospitalizedCurrently", "death", "deathConfirmed", "deathProbable")]
mydat <- mydat %>% mutate(total = positive + negative)

validstates = c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
mydat <- mydat %>% filter(state %in% validstates)
length(unique(mydat$state))

## lets convert the date column and sort the data
mydat = mydat %>% mutate(date = as.Date(as.character(date), "%Y%m%d")) 
mydat = mydat %>% arrange(state, desc(date)) 
mydat = mydat %>% filter(date <= ANALYSIS_DATE)

# compare covid-tracking with nytimes data
nyt_data = fread('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv')
nyt_data = nyt_data %>% group_by(state) %>% arrange(fips, desc(date)) %>% ungroup()
nyt_data = nyt_data %>% filter(date <= ANALYSIS_DATE)

# get nytimes data on the last date of downloaded 
bb = nyt_data %>% filter(date == max(date))
names(state.abb) <- state.name   
state.abb["District of Columbia"] = "DC"
bb = bb %>% mutate(abbr = state.abb[bb$state])
bb <- bb %>% filter(abbr %in% validstates)

aa = mydat %>% filter(date == max(date)) %>% select(date, state, death)

# check if the downloaded data ends on the same day.
unique(bb$date) == unique(aa$date)

# join the tables together
ab = aa %>% left_join(bb, by = c("state" = "abbr")) %>% 
  select(date.x, state, cvd=death, nyt=deaths) %>% 
  mutate(diff = abs(cvd - nyt))

# visual inspection shows covidtracking is same as nytimes
# except for new york. So we replace the nyc data in cvd with nytimes 
View(ab)

# fix NY data using NYtimes and CovidTracker
# nytimes NY data starts march 1, while covid tracking starts march 4.
# we match the covid tracking date since we need the rest of data also. 
mindate = mydat %>% filter(state == "NY") %>% select(date) %>% summarize(min(date)) %>% first
nyc_data_nyt = nyt_data %>% filter(state == "New York" & date >= mindate) %>% arrange(desc(date)) %>% pull(deaths)

# check size of vectors
length(nyc_data_nyt) == nrow(mydat %>% filter(state == "NY"))

# keep the old covid tracking data 
old_deathcol_ny = mydat[state == "NY"]$death

# replace the data
mydat[state == "NY"]$death = nyc_data_nyt
mydat[state == "NY"]$deathConfirmed = old_deathcol_ny
mydat[state == "NY"]$deathProbable = mydat[state == "NY"]$death - mydat[state == "NY"]$deathConfirmed

# replace NA with 0
#mydat = mydat %>% replace_na(list(death = 0, total = 0, positive=0))
# go through all states. If the deathConfirmed and deathProbably is ALL NA, then 
# death == deathConfirmed. 
for(st in validstates){
  deathT = mydat %>% filter(state == st) %>% pull(death)
  deathC = mydat %>% filter(state == st) %>% pull(deathConfirmed)
  deathP = mydat %>% filter(state == st) %>% pull(deathProbable)
  Ccheck = length(which(is.na(deathC))) == length(deathC)
  Pcheck = length(which(is.na(deathP))) == length(deathP)
  if (Ccheck && Pcheck){
    # replace the data
    mydat[state == st]$deathConfirmed = mydat[state == st]$death
    mydat[state == st]$deathProbable = mydat[state == st]$death - mydat[state == st]$deathConfirmed  
  }
}

# create an analysis worksheet 
# manual work, the death_analysis file is modified manually
# filling in all the missing information 
# either by cross referencing other sources or using regression methods
# this part is very time consuming and tedious. 
analysis_sheet = mydat %>% pivot_wider(id_cols = "date", names_from = state, values_from = c(death, deathConfirmed, deathProbable))
newcolorder = unlist(lapply(validstates, function(x) return(list(qq("death_@{x}"), qq("deathConfirmed_@{x}"), qq("deathProbable_@{x}")))))
analysis_sheet = analysis_sheet %>% select(all_of(newcolorder))
# lets add a date column arranged by date 
analysis_sheet$date = rev(seq( from = as.Date(ANALYSIS_DATE) - nrow(analysis_sheet) + 1, to = as.Date(ANALYSIS_DATE), by=1))
#fwrite(analysis_sheet, "death_analysis.csv") # save the data for interm analysis

# append the data Seyed/Affan fixed. This was done uptil October 9th. 
# once the new analysis is done, we'll have a new updated file to use for the next analysis. 
death_analysis_filled = fread("deaths_oct9_fixed.csv")
# the october 9th file dosn't have a date column so add it. 
#   -- I added it and re-saved the file for future - Affan Nov 26
#death_analysis_filled$date = rev(seq( from = as.Date("2020-10-09") - nrow(death_analysis_filled) + 1, to = as.Date("2020-10-09"), by=1))
#fwrite(death_analysis_filled, "deaths_oct9_fixed.csv")
nextdate = max(as.Date(death_analysis_filled$date)) + 1
as = analysis_sheet %>% filter(date >= nextdate)
# before appending, check if the columns are aligned
names(as) == names(death_analysis_filled)
# convert the date column to date type 
death_analysis_filled$date = as.Date(death_analysis_filled$date)
as_new = rbind(as, death_analysis_filled)
fwrite(as_new, "deaths_with_missing_data.csv") 


### BEFORE CONTINUING !!! 
# impute data deaths_with_missing_data.csv in excel, 
# then save as "deaths_DATE_imputed.csv" to persist
# THEN CONTINUE !!! 
# If reproducing, there should already be files checked in to continue analysis

# once the death analysis file is done, 
# read that in, and create a dataframe similar to deathdf 
# i.e. don't use mydat anymore. 
# death_df has colummns [date, states...]
# create two death_df (one for confirmed and one for total)
death_analysis_filled = fread("deaths_nov20_imputed.csv")
# if the top row is October 9, get the dates 
#death_analysis_filled$date = 
 # rev(seq( from = as.Date(ANALYSIS_DATE) - nrow(death_analysis_filled) + 1, to = as.Date("2020-10-09"), by=1))
# for november 20 dataset, the date column is already added, so just convert type
death_analysis_filled$date = as.Date(death_analysis_filled$date, format = "%m/%d/%Y") 

death_df <- data.table(death_analysis_filled)
#death_df = mydat %>% pivot_wider(id_cols = "date", names_from = state, values_from = death) %>% data.table
death_df_incidence = death_df %>% 
  select(-!!c("date")) %>%  
  map(make_incidence) %>% 
  as.data.frame
# add the date column back
death_df_incidence = cbind(date=death_df$date, death_df_incidence)

death_df[is.na(death_df)] <- 0
death_df_incidence[is.na(death_df_incidence)] <- 0

death_df[death_df < 0] <- 0
death_df_incidence[death_df_incidence < 0] <- 0

# check for negative numbers, should return zero. 
death_df_incidence %>% filter_all(any_vars(. < 0))

death_df_incidence <- death_df_incidence %>% relocate(date)
death_df_incidence <- death_df_incidence %>% arrange(date)

# might be worth checking the order of states here as well 
#names(death_df_incidence)[2:52] == validstates


total_df = mydat %>% pivot_wider(id_cols = "date", names_from = state, values_from = total) %>% data.table
total_df_incidence = total_df %>% 
  select(-!!c("date")) %>%  
  map(make_incidence) %>% 
  as.data.frame
# add the date column back
total_df_incidence = cbind(date=total_df$date, total_df_incidence)

total_df[is.na(total_df)] <- 0
total_df_incidence[is.na(total_df_incidence)] <- 0


positive_df = mydat %>% pivot_wider(id_cols = "date", names_from = state, values_from = positive) %>% data.table
positive_df_incidence = positive_df %>% 
  select(-!!c("date")) %>%  
  map(make_incidence) %>% 
  as.data.frame
# add the date column back
positive_df_incidence = cbind(date=positive_df$date, positive_df_incidence)

positive_df[is.na(positive_df)] <- 0
positive_df_incidence[is.na(positive_df_incidence)] <- 0


hosp_df = mydat %>% pivot_wider(id_cols = "date", names_from = state, values_from = hospitalizedCurrently) %>% data.table

# write the files. 
fwrite(total_df, "/data/actualdeaths_covid19/downloaded_data/cumulative_total_tests_lancetid.csv")
fwrite(positive_df, "/data/actualdeaths_covid19/downloaded_data/cumulative_positive_cases_lancetid.csv")
fwrite(death_df, "/data/actualdeaths_covid19/downloaded_data/cumulative_deaths_lancetid.csv")
fwrite(hosp_df, "/data/actualdeaths_covid19/downloaded_data/hosp_data.csv")
#fwrite(hosp_df, "~/actual_covid19_deaths/hosp_data.csv")

# write the files. 
fwrite(total_df_incidence, "/data/actualdeaths_covid19/downloaded_data/incidence_total_tests_lancetid.csv")
fwrite(positive_df_incidence, "/data/actualdeaths_covid19/downloaded_data/incidence_positive_cases_lancetid.csv")
fwrite(death_df_incidence, "/data/actualdeaths_covid19/downloaded_data/incidence_deaths_lancetid.csv")

## split the death files into confirmed and probably 
## this is useful for Seyed's plotting
all_colmns = unlist(map(validstates, function(x) qq("death_@{x}")))
cnf_colmns = unlist(map(validstates, function(x) qq("deathConfirmed_@{x}")))
all_dd = death_df_incidence %>% select(all_colmns)
cnf_dd = death_df_incidence %>% select(cnf_colmns)
fwrite(all_dd, "/data/actualdeaths_covid19/downloaded_data/total_deaths.csv")
fwrite(cnf_dd, "/data/actualdeaths_covid19/downloaded_data/confirmed_deaths.csv")



