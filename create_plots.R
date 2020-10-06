## MCMC analysis for estimating the true COVID-19 deaths 
## Affan Shoukat 2020

## to do: eventually move everything to Julia
## creates a ggplot for all states after the model is run. 
## -- uses functions make_a_plot() and create_all_plots() 
## creates a table file (that can be further processed in Excel) of the results
## -- mean, lo, and high CI 
## -- uses functions create_table() and create_table_states()

## NOTE: source main_model_run.R to bring in functions and data. 
## In particular, the function `get_state_data_vectors` and vector `validstates` are required


library(tidyverse) # not on the cluster
library(bayestestR)

extract_data <- function(st){
  deathdata = get_state_data_vectors(st, ma=F)$death
  
  posterior_y = fread(qq("/data/actualdeaths_covid19/st_@{st}_00_posterior_y.dat"))
  posterior_z = fread(qq("/data/actualdeaths_covid19/st_@{st}_00_posterior_z.dat"))
  mns1 = apply(posterior_y, 2, mean)
  mns2 = apply(posterior_z, 2, mean)
  
 
  c_data = sum(deathdata)
  c_mns1 = round(sum(mns1), 2) # estimated true
  c_mns2 = round(sum(mns2), 2) # estimated fit 
  c_pinc = round((c_mns1 - c_mns2)/c_mns2, 2) # % inc from fit 
  c_dinc = round((c_mns1 - c_data)/c_data, 2) # % inc from data 

  # idea is to use the %increase from fit to true 
  # and use that % to translate to increase in the reported number
  d_inc = c_data * (1 + c_pinc)
  return(list(deathdata=deathdata,mns1=mns1, mns2=mns2, c_data=c_data,c_mns1=c_mns1,c_mns2=c_mns2, c_pinc=c_pinc, c_dinc=c_dinc, d_inc=d_inc))
}

make_a_plot <- function(st){
  gdat = extract_data(st)
  deathdata = gdat$deathdata
  mns1=gdat$mns1; mns2=gdat$mns2
  c_data = gdat$c_data
  c_mns1 = gdat$c_mns1
  c_mns2 = gdat$c_mns2
  c_pinc = gdat$c_pinc
  c_dinc = gdat$c_dinc
  d_inc  = gdat$d_inc
  
  
  diff_dinc = round(c_mns1 - d_inc, 0)
  cdc_pexcess = round(as.numeric(fmeans %>% filter(abbr == st) %>% select("percent_val") ), 2)
  if(is.na(cdc_pexcess))
    cdc_pexcess = "NA"
  #c_finc = round((c_mns1 - sum(poly_tim))/sum(poly_tim), 2)
  c_Str = qq(" st: @{st}, data: @{c_data}, fit: @{c_mns2}, red: @{c_mns1} \n %inc (data): @{c_dinc} %inc (fit): @{c_pinc} \n data inc: @{d_inc}, diff: @{diff_dinc} \n cdc: @{cdc_pexcess}")
  print(c_Str)
  xvals = length(mns1)
  gg = ggplot()
  gg = gg + geom_point(aes(x=1:xvals, y=deathdata), color="#636363")
  gg = gg + geom_line(aes(x=1:xvals, y=deathdata), color="#636363")
  gg = gg + geom_line(aes(x=1:xvals, y=mns1), color="#e41a1c", size=1.2)
  #gg = gg + ylim(0, 500)
  gg = gg + geom_line(aes(x=1:xvals, y=mns2), color="#377eb8", size=1.2)
  #gg = gg + geom_line(aes(x=1:xvals, y=poly_tim), color="#4daf4a", size=1.2, alpha=0.8)
  gg = gg + annotate("text", -Inf, Inf, label = c_Str, hjust = 0, vjust = 1)
  gg = gg + xlab("time") + ylab("deaths")
  gg
}

create_all_plots <- function(){
  mplots = map(validstates, make_a_plot)
  #a_large_plot = ggpubr::ggarrange(plotlist = mplots, ncol=4)
  a_large_plot = cowplot::plot_grid(plotlist = mplots, ncol = 4)
  ggsave(filename = qq("all_states_lancetid.pdf"), plot=a_large_plot, width=20.5, height=40.5, units="in")
}

## next two functions looks at the posterior_y/posterior_z files and 
## creates a table of results (mean values + HDI intervals)
create_table <- function(st){
  # read all posterior_y files and create a data table
  gdat = extract_data(st)
  
  deathdata = gdat$deathdata
  postdata  = round(gdat$d_inc)
  posty = round(gdat$c_mns1)
  
  # cant use gdat to create the daily intervals 
  posterior_y = fread(qq("/data/actualdeaths_covid19/st_@{st}_00_posterior_y.dat"))
  samp_sums = apply(posterior_y, 1, sum)
  cred_int = ci(samp_sums, method = "HDI",ci = 0.95)
  lo = cred_int[[2]]
  hi = cred_int[[3]]
  #hist(samp_sums)
  c_Str = qq(" st: @{st}, mean: @{posty} lo: @{lo}, hi: @{hi}")
  print(c_Str)
  
  return(list(st=st, death=sum(deathdata), posty = posty, postdata=postdata,
              lo=lo, hi=hi))
}

create_table_states <- function(){
  aa = data.table(rbindlist(map(validstates, create_table)))
  aa = aa %>% left_join(select(city_metadata, c(state, abbr)), by = c("st" = "abbr"))
  aa = aa %>% relocate(state, .before = death)
  fwrite(aa, "table_results.csv")
}

