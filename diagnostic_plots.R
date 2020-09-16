## MCMC analysis for estimating the true COVID-19 deaths 
## Affan Shoukat 2020
## scratch file: not really used for the main analysis. 

rm(list = ls())
library(ggplot2)
library(ggmcmc)
#library(ggpubr)
library(data.table)
library(reshape)
library(dplyr)
library(ggthemes)
library(gridExtra)
library(GetoptLong)
library(stringr)


#pmeans <- pmeans %>% filter(abbr %in% levels(df$arg_st))
## delete columns not needed
#pmeans = pmeans %>%
#  mutate(abbr =  factor(abbr, levels = levels(df$arg_st))) %>%
#  arrange(abbr)   


## appendix plots 
## plotting the qgam fit for NY, California, NJ, and FL 

ny_dat = deathdata_ma
ny_fit = poly_tim

ca_dat = deathdata_ma
ca_fit = poly_tim

nj_dat = deathdata_ma
nj_fit = poly_tim

fl_dat = deathdata_ma 
fl_fit = poly_tim

appd.gg <- function(dat, fit){
  gg = ggplot() + theme_bw() + ylab("Deaths")
  #gg <- gg + theme(plot.margin=unit(c(0.1,0,0,0.1), "null"))
  #gg <- gg + theme(panel.spacing=unit(c(-5,0,0,0), "null"))
  
  xvals = 1:length(dat)
  xlabs = c(1, seq(from=0, to=length(xvals), by = 10))
  gg = gg + geom_point(aes(x=xvals, y=dat, color = "black"),  size=2, shape=5)
  gg = gg + geom_line(aes(x=xvals, y=fit,  color="blue"), size=1.2, alpha=0.8)
  gg = gg + scale_x_continuous("Time", labels = as.character(xlabs), breaks = xlabs, expand=c(0, 0))
  gg = gg +  scale_colour_manual(name=NULL, values=c("black", "blue"), labels = c("raw data", "qgam fitted"), 
                            guide = guide_legend(override.aes = list(
                               linetype = c("blank", "solid"),
                               shape = c(5, NA))))
  gg = gg + theme(legend.position = c(0.95, 0.95), legend.justification = c(1, 1))
  gg = gg +  theme(legend.key.size = unit(0.5, "cm"))
  gg
}
# remove the legends from all of them and the y labels for the columns
ny.gg <- appd.gg(ny_dat, ny_fit) +  theme(legend.position="none")
nj.gg <- appd.gg(nj_dat, nj_fit) + ylab(NULL)
ca.gg <- appd.gg(ca_dat, ca_fit) + theme(legend.position="none")
fl.gg <- appd.gg(fl_dat, fl_fit) + ylab(NULL) + theme(legend.position="none")

prow = plot_grid(ny.gg, nj.gg, ca.gg, fl.gg, labels=c("NY", "NJ", "CA", "FL"),
                 ncol = 2, nrow = 2, hjust=0)
prow
# extract the legend from one of the plots
# legend <- get_legend(
#   # create some space to the left of the legend
#   appd.gg(ny_dat, ny_fit) + 
#     theme(legend.box.margin = margin(0, 0, 0, 12)) 
#     #guides(color = guide_legend(nrow = 1)) +
#     #theme(legend.position = "bottom")
# )
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
#plot_grid(prow, legend, rel_widths = c(2, 0.1))


#### read trace plots 
arg_st = "NY"
#color_scheme_set("blue")
#mcmc_trace(posterior, pars = c("wt", "sigma"))
trace_dt = fread(qq("/data/actualdeaths_covid19/st_@{arg_st}_00_traceplots.dat"))
trace_dt_gg <- function(param){
  xvals = seq(from=0, to=12000, by=2000)
  gg = ggplot() + ylab(qq("Parameter: @{param}")) + theme_bw()
  gg <- gg + theme(plot.margin=unit(c(0,0.05,0,0), "null"))
  #gg <- gg + theme(panel.spacing=unit(c(-5,0,0,0), "null"))
  gg = gg + geom_line(aes(x=1:12000, y=trace_dt[, get(param)]), size=0.1)
  gg = gg + scale_x_continuous(name = NULL, breaks=xvals, expand=c(0, 0))
  return(gg)
}
gg_a0 = trace_dt_gg("a0")
gg_a1 = trace_dt_gg("a1")
gg_b0 = trace_dt_gg("b0")
gg_b1 = trace_dt_gg("b1")
gg_b2 = trace_dt_gg("b2")

p_top_row = plot_grid(gg_a0, gg_a1, ncol = 2, hjust=0)
p_bot_row = plot_grid(gg_b0, gg_b1, gg_b2, ncol = 3, hjust=0)
cowplot::plot_grid(p_top_row, p_bot_row, ncol=1, rel_heights=c(0.5,0.5,0.5))
