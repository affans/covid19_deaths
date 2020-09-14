## MCMC analysis for estimating the true COVID-19 deaths 
## Downloaded files here are used by Seyed/Ellie to generate CFR and UR incidence. 
## The file hosp_data.csv is not sent to Seyed, but still used in the Bayesian model

rm(list=ls())
library(data.table)
library(tidyverse)

make_incidence <- function(cum_data){
  len = length(cum_data)
  incidence = numeric(len)
  for (i in 1:(len-1)) {
    incidence[i] = cum_data[i] - cum_data[i+1]
  }
  incidence[len] = 0
  return(incidence)
}
mydat <- fread('https://covidtracking.com/api/v1/states/daily.csv')
mydat <- mydat[, c(1:4, 6, 17)] # col 17 is death
mydat <- mydat %>% mutate(total = positive + negative)

validstates = c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
mydat <- mydat %>% filter(state %in% validstates)
length(unique(mydat$state))

## how many times points?
length(unique(mydat$date))

# none of these long -> wide work
# reshape probably works but row numbers are useless
# reshape(mydat, idvar="date", timevar="state", v.names="total", direction="wide", sep="_")
# spread is the wrong function for this. 
# mydat %>% spread(state, total)

# pivot_wider and dcast work
total_df = mydat %>% pivot_wider(id_cols = "date", names_from = state, values_from = total) %>% data.table
total_df_incidence = total_df %>% 
  select(-!!c("date")) %>%  
  map(make_incidence) %>% 
  as.data.frame
# add the date column back
total_df_incidence = cbind(date=total_df$date, total_df_incidence)

positive_df = mydat %>% pivot_wider(id_cols = "date", names_from = state, values_from = positive) %>% data.table
positive_df_incidence = positive_df %>% 
  select(-!!c("date")) %>%  
  map(make_incidence) %>% 
  as.data.frame
# add the date column back
positive_df_incidence = cbind(date=positive_df$date, positive_df_incidence)

death_df = mydat %>% pivot_wider(id_cols = "date", names_from = state, values_from = death) %>% data.table
death_df_incidence = death_df %>% 
  select(-!!c("date")) %>%  
  map(make_incidence) %>% 
  as.data.frame
# add the date column back
death_df_incidence = cbind(date=death_df$date, death_df_incidence)

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



