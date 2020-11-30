

# paste death figures together
library(tidyverse) 
library(magick)

load("output_data/country_deaths.RData")

total_rayrender_frames <- length(dir(path = "figures/frames",pattern = "rayrender_deaths_estimated"))


for(image_index in 1:total_rayrender_frames){
  
  
  line_plot <- magick::image_read(paste0("figures/frames/cumulative_deaths_ggplot_line",image_index,".png"))
  spatial_ggplot_rep <- magick::image_read(paste0("figures/frames/cumulative_deaths_ggplot_spatial_reported",image_index,".png"))
  spatial_ggplot_est <- magick::image_read(paste0("figures/frames/cumulative_deaths_ggplot_spatial_estimated",image_index,".png"))
  spatial_rayrender_rep <- magick::image_read(paste0("figures/frames/rayrender_deaths_reported",image_index,".png")) 
  spatial_rayrender_est <- magick::image_read(paste0("figures/frames/rayrender_deaths_estimated",image_index,".png")) 
  
  
  composite_image_all <- magick::image_append(c(magick::image_resize(spatial_rayrender_rep,"x1000"),
                                                magick::image_resize(spatial_rayrender_est,"x1000"))) %>% 
    # magick::image_composite(.,magick::image_resize(line_plot,"x300"),offset = "+635") %>%
    magick::image_composite(.,magick::image_resize(line_plot,"x300"),offset = "+700") %>%
    magick::image_composite(., magick::image_resize(spatial_ggplot_rep,"x300"),offset = "+2+25") %>%
    magick::image_composite(., magick::image_resize(spatial_ggplot_est,"x300"),offset = "+1400+25") %>%
    magick::image_annotate(., paste0("Date: ",unique(deaths_country_m$date)[image_index]), 
                           size = 25, location = "+1750+950", color = "white") %>%
    magick::image_annotate(., paste0("Shoukat, A., Moghadas, S.M.,  Fitzpatrick, M.C., Wells, C.R., \nAbdollahi, E.,  Pandey, A., Cannataro, V.L.,  Galvani, A.P. (2020) \nEstimating the number of COVID-19 deaths in the United States"), 
                           size = 15, location = "+1+925", color = "white") %>%
    magick::image_annotate(., paste0("Reported deaths\nper 1,000 people in each state"), 
                           size = 25, location = "+610+775", color = "white") %>%
    magick::image_annotate(., paste0("Estimated deaths\nper 1,000 people in each state"), 
                           size = 25, location = "+1610+775", color = "white")
  
    magick::image_write(image = composite_image_all,path = paste0("figures/frames/composite_image_all_deaths_",image_index,".png"))
  
  
  
  
  # composite_image_all <- magick::image_composite(magick::image_resize(spatial_rayrender_rep,"x1000"),
  #                                                magick::image_resize(spatial_rayrender_est,"x1000"),offset = "+1000") %>% 
  #   magick::image_composite(., magick::image_resize(spatial_ggplot,"x250"),offset = "+510+25") %>%
  #   magick::image_annotate(., paste0("Date: ",unique(estimated_and_merged_m$date)[image_index]), 
  #                          size = 25, location = "+800+900", color = "white") %>%
  #   magick::image_annotate(., paste0("Galvani et al. 2020 [important information here]"), 
  #                          size = 15, location = "+1+925", color = "white") 
  # magick::image_write(image = composite_image_all,path = paste0("figures/frames/composite_image_all_deaths",image_index,".png"))
  print(image_index)
  gc()
}

av::av_encode_video(glue::glue("figures/frames/composite_image_all_deaths_{1:(total_rayrender_frames)}.png"),
                    framerate=8, output = "figures/estimated_deaths_USA.mp4",
                    vfilter = "pad=ceil(iw/2)*2:ceil(ih/2)*2")
