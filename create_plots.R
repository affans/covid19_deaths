## MCMC analysis for estimating the true COVID-19 deaths 
## Affan Shoukat 2020

## creates a ggplot for the posterior y/z results.
## NOTE: THIS FILE REQUIRES main_model_run to be sourced. 
## In particular, the function `get_state_data_vectors` and vector `validstates` are required

library(tidyverse) # not on the cluster
library(bayestestR)

validstates = c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
subfldr <- ""

make_a_plot <- function(st){
  deathdata = get_state_data_vectors(st, ma=F)$death
  
  posterior_y = fread(qq("/data/actualdeaths_covid19/@{subfldr}/st_@{st}_00_posterior_y.dat"))
  posterior_z = fread(qq("/data/actualdeaths_covid19/@{subfldr}/st_@{st}_00_posterior_z.dat"))
  mns1 = apply(posterior_y, 2, mean)
  mns2 = apply(posterior_z, 2, mean)
  
  xvals = length(mns1)
  c_data = sum(deathdata)
  c_mns1 = round(sum(mns1), 2)
  c_mns2 = round(sum(mns2), 2)
  c_pinc = round((c_mns1 - c_data)/c_data, 2)
  cdc_pexcess = round(as.numeric(fmeans %>% filter(abbr == st) %>% select("percent_val") ), 2)
  if(is.na(cdc_pexcess))
    cdc_pexcess = "NA"
  #c_finc = round((c_mns1 - sum(poly_tim))/sum(poly_tim), 2)
  c_Str = qq(" st: @{st}, data: @{c_data}, true: @{c_mns1} \n obs: @{c_mns2}, %inc (data): @{c_pinc} \n cdc: @{cdc_pexcess}")
  print(c_Str)
  
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
mplots = map(validstates, make_a_plot)
#a_large_plot = ggpubr::ggarrange(plotlist = mplots, ncol=4)
a_large_plot = cowplot::plot_grid(plotlist = mplots, ncol = 4)
ggsave(filename = qq("all_states_lancetid.pdf"), plot=a_large_plot, width=15.5, height=30.5, units="in")



create_table <- function(st){
  # read all posterior_y files and create a data table
  deathdata = get_state_data_vectors(st, ma=F)$death
  csum = sum(deathdata)
  
  posterior_y = fread(qq("/data/actualdeaths_covid19/@{subfldr}/st_@{st}_00_posterior_y.dat"))
  samp_sums = apply(posterior_y, 1, sum)
  posty = round(mean(samp_sums))
  cred_int = ci(samp_sums, method = "HDI",ci = 0.95)
  lo = cred_int[[2]]
  hi = cred_int[[3]]
  hist(samp_sums)
  print(qq("@{st}"))
  return(list(st=st, death=sum(deathdata), posty = posty, 
              lo=lo, hi=hi))
}

aa = data.table(rbindlist(map(validstates, create_table)))
aa = aa %>% left_join(select(city_metadata, c(state, abbr)), by = c("st" = "abbr"))
aa = aa %>% relocate(state, .before = death)
fwrite(aa, "table_results.csv")


