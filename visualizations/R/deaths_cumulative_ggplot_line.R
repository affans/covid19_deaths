

# death cumulative line plot ---- 

load("output_data/country_deaths.RData")
library(tidyverse);library(scales)

head(deaths_country_m)

# states_geom_and_data_m 


# cumulative_state_death <- states_geom_and_data %>% 
#   arrange((date)) %>%
#   group_by(state,reported_or_estimated) %>% 
#   mutate(cumulative_deaths = cumsum(deaths)) %>% 
#   arrange(state, date)


# 
# 
# cumulative_deaths <- cumulative_state_death %>%
#   group_by(date,reported_or_estimated) %>%
#   summarize(total_deaths = sum(cumulative_deaths)) %>%
#   ungroup() # %>%
# # pivot_wider(id_cols = date,names_from = reported_or_estimated,values_from= total_deaths)


fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}


cumulative_deaths <- deaths_country_m %>%
  mutate(`Death type` = case_when(`Death type` == "deaths_per_day_est_cum" ~ "Estimated", 
                                  `Death type` == "deaths_per_day_cum" ~ "Reported"))


cumulative_deaths$`Death type` <- factor(cumulative_deaths$`Death type`, levels = (c("Reported","Estimated")))

for(date_ind in 1:length(unique(cumulative_deaths$date))){
  
  this_date <- unique(cumulative_deaths$date)[date_ind]
  
  ggplot_line_plot_daily <- ggplot(cumulative_deaths) + 
    geom_line(aes(x=date, y=`Cumulative deaths`,color=`Death type`),size=3) +
    geom_segment(data = subset(cumulative_deaths, date == this_date),
                 aes(x=this_date,xend=this_date,y=0,yend=`Cumulative deaths`,color=`Death type`),size=2) +
    geom_point(data = subset(cumulative_deaths, date == this_date),
               aes(x=this_date,y=`Cumulative deaths`,color=`Death type`),size=5) + 
    labs(y="Total COVID-19 deaths",
         x= "Date") + 
    scale_y_continuous(labels = fancy_scientific) + 
    scale_x_date(date_breaks = "1 months",labels = date_format("%b")) +
    scale_color_manual(values = setNames(object = c("red","blue"),
                                         nm = c("Estimated","Reported"))) +
    theme(panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          text = element_text(size=18,color="white"),
          axis.text = element_text(size=18,color="white"),
          line = element_line(size=0.25,color="white"),
          legend.background = element_rect(fill = "transparent",colour = NA),
          legend.title=element_blank(),legend.position = "bottom",
          legend.key = element_rect(colour = "transparent", fill=NA)) 
  
  ggsave(filename = paste0("figures/frames/cumulative_deaths_ggplot_line",date_ind,".png"),
         plot = ggplot_line_plot_daily,bg = "transparent",height = 4,width = 8)
  print(date_ind)
}
