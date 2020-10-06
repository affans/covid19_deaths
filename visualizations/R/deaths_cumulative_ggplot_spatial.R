# death cumulative line plot ---- 

load("output_data/states_geom_and_data.RData")
library(tidyverse)

per_pop <- 1e3

states_and_geom_estimate <- states_geom_and_data %>% 
  # mutate(deaths_per_person_log10 = log10(deaths_per_person)) %>%
  # mutate(deaths_per_person_log10 = case_when(deaths_per_person_log10 != -Inf ~ deaths_per_person_log10)) %>%
  mutate(deaths_per_day_cum_per_person_per_10 = (deaths_per_day_cum * per_pop)/estimate) %>%
  mutate(deaths_per_day_est_cum_per_person_per_10 = (deaths_per_day_est_cum * per_pop)/estimate) %>%
  # mutate(deaths_binned_reported = cut(deaths_per_day_cum_per_person_per_10,breaks=c(0,0.01,0.1,0.5,1,2,3,4,5),include.lowest = T)) %>%
  # mutate(deaths_binned_estimated = cut(deaths_per_day_est_cum_per_person_per_10,breaks=c(0,0.01,0.1,0.5,1,2,3,4,5),include.lowest = T))
  mutate(deaths_binned_reported = cut(deaths_per_day_cum_per_person_per_10,breaks=c(0,0.005,0.01,0.05,0.1,0.5,1,2,3),include.lowest = T)) %>%
  mutate(deaths_binned_estimated = cut(deaths_per_day_est_cum_per_person_per_10,breaks=c(0,0.005,0.01,0.05,0.1,0.5,1,2,3),include.lowest = T))


cols_in_plot_estimated <- setNames(
  object = scales::seq_gradient_pal(
    low="#EECB8B",high = "darkred")(seq(0,1,length.out = length(unique(states_and_geom_estimate$deaths_binned_estimated)))),
  nm = levels(states_and_geom_estimate$deaths_binned_estimated))


cols_in_plot_reported <- setNames(
  object = scales::seq_gradient_pal(
    low="#C3D3FB",high = "darkblue")(seq(0,1,length.out = length(unique(states_and_geom_estimate$deaths_binned_estimated)))),
  nm = levels(states_and_geom_estimate$deaths_binned_estimated))




plotter_function <- function(date_ind){
  
  data_for_spatial_ggplot <- states_and_geom_estimate %>% 
    filter(date == unique(states_and_geom_estimate$date)[date_ind])
  
  daily_ggplot_reported <- ggplot(data_for_spatial_ggplot) + 
    geom_sf(aes(fill=deaths_binned_reported)) + 
    scale_fill_manual(values = cols_in_plot_reported, 
                      name="Reported deaths\nper 1,000 people",  
                      guide = guide_legend(reverse = TRUE),drop=FALSE,
                      labels = c(
                        "[0–0.005]",
                        "(0.005–0.01]",
                        "(0.01–0.05]",
                        "(0.05–0.1]",
                        "(0.1–0.5]",
                        "(0.5–1]",
                        "(1–2]",
                        "(2–3]")
    ) + 
    theme(panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA),
          text = element_text(size=18,color="white"),
          axis.text = element_text(size=18,color="white"),
          line = element_line(size=0.25,color="white"),
          legend.text = element_text(size=18,color="white"),
          legend.title = element_text(size=18))
  ggsave(filename = paste0("figures/frames/cumulative_deaths_ggplot_spatial_reported",date_ind,".png"),
         plot = daily_ggplot_reported,bg = "transparent",width = 8, height = 4)
  
  
  
  daily_ggplot_estimated <- ggplot(data_for_spatial_ggplot) + 
    geom_sf(aes(fill=deaths_binned_estimated)) + 
    scale_fill_manual(values = cols_in_plot_estimated, 
                      name="Estimated deaths\nper 1,000 people",  
                      guide = guide_legend(reverse = TRUE),drop=FALSE,
                      labels = c(
                        "[0–0.005]",
                        "(0.005–0.01]",
                        "(0.01–0.05]",
                        "(0.05–0.1]",
                        "(0.1–0.5]",
                        "(0.5–1]",
                        "(1–2]",
                        "(2–3]")
    ) + 
    theme(panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA),
          text = element_text(size=18,color="white"),
          axis.text = element_text(size=18,color="white"),
          line = element_line(size=0.25,color="white"),
          legend.text = element_text(size=18,color="white"),
          legend.title = element_text(size=18))
  ggsave(filename = paste0("figures/frames/cumulative_deaths_ggplot_spatial_estimated",date_ind,".png"),
         plot = daily_ggplot_estimated,bg = "transparent",width = 8, height = 4)
  
  
  # print(date_ind)
}


# parallel::mclapply(X = 1:230,FUN = plotter_function,mc.cores = 6)

# pbapply::pblapply(X = 1:230,FUN = plotter_function,cl = 6)





for(date_ind in 1:length(unique(states_and_geom_estimate$date))){

  data_for_spatial_ggplot <- states_and_geom_estimate %>%
    filter(date == unique(states_and_geom_estimate$date)[date_ind])

  daily_ggplot_reported <- ggplot(data_for_spatial_ggplot) +
    geom_sf(aes(fill=deaths_binned_reported)) +
    scale_fill_manual(values = cols_in_plot_reported,
                      name="Reported deaths\nper 1,000 people",
                      guide = guide_legend(reverse = TRUE),drop=FALSE,
                      labels = c(
                        "[0–0.005]",
                        "(0.005–0.01]",
                        "(0.01–0.05]",
                        "(0.05–0.1]",
                        "(0.1–0.5]",
                        "(0.5–1]",
                        "(1–2]",
                        "(2–3]")
    ) +
    theme(panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA),
          text = element_text(size=18,color="white"),
          axis.text = element_text(size=18,color="white"),
          line = element_line(size=0.25,color="white"),
          legend.text = element_text(size=18,color="white"),
          legend.title = element_text(size=18))
  ggsave(filename = paste0("figures/frames/cumulative_deaths_ggplot_spatial_reported",date_ind,".png"),
         plot = daily_ggplot_reported,bg = "transparent",width = 8, height = 4)



  daily_ggplot_estimated <- ggplot(data_for_spatial_ggplot) +
    geom_sf(aes(fill=deaths_binned_estimated)) +
    scale_fill_manual(values = cols_in_plot_estimated,
                      name="Estimated deaths\nper 1,000 people",
                      guide = guide_legend(reverse = TRUE),drop=FALSE,
                      labels = c(
                        "[0–0.005]",
                        "(0.005–0.01]",
                        "(0.01–0.05]",
                        "(0.05–0.1]",
                        "(0.1–0.5]",
                        "(0.5–1]",
                        "(1–2]",
                        "(2–3]")
    ) +
    theme(panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA),
          text = element_text(size=18,color="white"),
          axis.text = element_text(size=18,color="white"),
          line = element_line(size=0.25,color="white"),
          legend.text = element_text(size=18,color="white"),
          legend.title = element_text(size=18))
  ggsave(filename = paste0("figures/frames/cumulative_deaths_ggplot_spatial_estimated",date_ind,".png"),
         plot = daily_ggplot_estimated,bg = "transparent",width = 8, height = 4)


  print(date_ind)
}
