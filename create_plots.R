## this is a helper file to USA_UR_Deaths_LancetID to create plots of the data
## MCMC analysis for estimating the true COVID-19 deaths 
## See file USA_UR_Deaths_LancetID.R for details.

## Affan Shoukat
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

totaldeaths <<- 0 
totaldeaths_red <<- 0 

make_a_plot <- function(st){
  deathdata = get_cfr_death_data(st, ma=F)$death
  
  posterior_y = fread(qq("/data/actualdeaths_covid19/st_@{st}_00_posterior_y.dat"))
  posterior_z = fread(qq("/data/actualdeaths_covid19/st_@{st}_00_posterior_z.dat"))
  mns1 = apply(posterior_y, 2, mean)
  mns2 = apply(posterior_z, 2, mean)
  
  xvals = length(mns1)
  c_data = sum(deathdata)
  c_mns1 = round(sum(mns1), 2)
  c_mns2 = round(sum(mns2), 2)
  c_pinc = round((c_mns1 - c_data)/c_data, 2)
  #c_finc = round((c_mns1 - sum(poly_tim))/sum(poly_tim), 2)
  c_Str = qq(" st: @{st}, data: @{c_data}, true: @{c_mns1} \n obs: @{c_mns2}, %inc (data): @{c_pinc}")
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






