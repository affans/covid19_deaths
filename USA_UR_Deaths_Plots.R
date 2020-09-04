## this is a helper file to USA_UR_Deaths_LancetID to create plots of the data
## MCMC analysis for estimating the true COVID-19 deaths 
## See file USA_UR_Deaths_LancetID.R for details.

## Affan Shoukat
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

validstates = c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
#validstates = c("AL", "AK", "AZ" ,"AR" ,"CA" ,"CO" ,"CT" ,"DE")

# read death data
all_death_data = fread("Death.csv")
names(all_death_data) <- validstates

get_cfr_death_data <- function(st){
  ## get the data data, and remove the zeros
  ddd = all_death_data[, get(st)]
  firstnonzero = min( which ( ddd != 0 ))
  lastelement = length(ddd)
  death_rm = ddd[firstnonzero : lastelement]
  return(death_rm)
}

totaldeaths <<- 0 
totaldeaths_red <<- 0 

make_a_plot <- function(st){
  dd =  get_cfr_death_data(st)
  
  posterior_z = fread(qq("/data/actualdeaths_covid19/st_@{st}_posterior_z.dat"))
  mns1 = apply(posterior_z, 2, mean)
  
  
  posterior_y = fread(qq("/data/actualdeaths_covid19/st_@{st}_posterior_y.dat"))
  mns2 = apply(posterior_y, 2, mean)
  UR_N = length(mns1)
  
  c_dd = sum(dd)
  c_mns1 = sum(mns1)
  c_mns2 = sum(mns2)
  totaldeaths <<- totaldeaths + c_mns2
  totaldeaths_red <<- totaldeaths_red + c_dd
  c_Str = qq("data: @{c_dd}, fit: @{c_mns1}, blue: @{c_mns2}")
  gg = ggplot()
  gg = gg + annotate("text", -Inf, Inf, label = c_Str, hjust = 0, vjust = 1)
  gg = gg + geom_line(aes(x=1:UR_N, y=dd), color="black")
  gg = gg + geom_line(aes(x=1:UR_N, y=mns1), color="red")
  gg = gg + geom_line(aes(x=1:UR_N, y=mns2), color="blue")
  gg = gg + ggtitle(st)
  gg
}

print(qq("total: @{totaldeaths_red}"))

mplots = map(validstates, make_a_plot)
#a_large_plot = ggpubr::ggarrange(plotlist = mplots, ncol=4)
a_large_plot = cowplot::plot_grid(plotlist = mplots, ncol = 4)
ggsave(filename = qq("all_states_lancetid.pdf"), plot=a_large_plot, width=12.5, height=20.5, units="in")


get_chad_request <- function(st, find_date){
  posterior_z = fread(qq("/data/lancetid_actualdeaths_covid19/st_@{st}_posterior_z.dat"))
  mns1 = apply(posterior_z, 2, mean)
  posterior_y = fread(qq("/data/lancetid_actualdeaths_covid19/st_@{st}_posterior_y.dat"))
  mns2 = apply(posterior_y, 2, mean)
  
  # start date is january 22. 
  # end date was july 3
  # thats 164 elements, but the first xx elements are zero so they are removed.
  ddd = all_death_data[, get(st)]
  firstnonzero = min( which ( ddd != 0 ))
  
  # so january 22 + firstnonzero gives the start date for the state. 
  start_date = as.Date("2020-01-22")
  data_date = start_date + firstnonzero
  #find_date = as.Date("2020-04-26")
  diff_date = as.Date(find_date) - data_date
  smns1 = sum(mns1[1:diff_date])
  smns2 = sum(mns2[1:diff_date])
  print(qq("st: @{st}, date: @{find_date}, estimated: @{smns1}"))
}






