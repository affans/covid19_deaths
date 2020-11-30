#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse);library(lubridate)
library(sf)
library(rayrender) 

# rayrender has great tutorials, including a gist that shows how to 
# extrude polygon states along a color spectrum 
# https://gist.github.com/tylermorganwall/19b2a39ee881e27d42ba5661f0867052 

date_ind <- as.numeric(args[1])

load("output_data/states_geom_and_data.RData")

# for(date_ind in seq_along(unique(states_and_geom_estimate$date))){

this_date <- unique(states_geom_and_data$date)[date_ind]


message(this_date)
this_sf <- states_geom_and_data %>%
  dplyr::filter(date == lubridate::ymd(this_date)) %>%
  mutate(deaths_per_day_cum_per_person = deaths_per_day_cum_per_person + 1e-7,
         deaths_per_day_est_cum_per_person = deaths_per_day_est_cum_per_person + 1e-7)

# 70A1D1 sent to colleagues 
reported_deaths_scene <-
  xz_rect(y=-0.01,zwidth=500,xwidth=500,
          material = diffuse(color="grey3")) %>%
  add_object(extruded_polygon(this_sf, center=TRUE,z=-5,scale=c(1/100000,1,1/100000),
                              data_column_top = "deaths_per_day_cum_per_person",
                              scale_data = 1/max(states_geom_and_data$deaths_per_day_est_cum_per_person,na.rm=TRUE)*30,
                              material=glossy(
                                color="#4657F6",
                                # gradient_color = "#fcde9c",
                                # color="#7c1d6f",
                                gradient_color="darkblue",
                                gradient_point_start = c(0,0,0),
                                gradient_point_end = c(0,7,0),
                                reflectance=0,
                                gloss = 0.3))) %>%
  # material=diffuse(color="#754315",
  #                 gradient_color="darkred",
  #                 gradient_point_start = c(0,0,0),
  #                 gradient_point_end = c(0,8,0)))) %>%
  add_object(sphere(y=200,radius=80,material=light()))

# save(test_scene, file="output_data/test_scene_cases.RData")



reported_deaths_scene_red <-
  xz_rect(y=-0.01,zwidth=500,xwidth=500,
          material = diffuse(color="grey3")) %>%
  add_object(extruded_polygon(this_sf, center=TRUE,z=-5,scale=c(1/100000,1,1/100000),
                              data_column_top = "deaths_per_day_cum_per_person",
                              scale_data = 1/max(states_geom_and_data$deaths_per_day_est_cum_per_person,na.rm=TRUE)*30,
                              material=glossy(
                                color="#754315",
                                # gradient_color = "#fcde9c",
                                # color="#7c1d6f",
                                gradient_color="darkred",
                                gradient_point_start = c(0,0,0),
                                gradient_point_end = c(0,7,0),
                                reflectance=0,
                                gloss = 0.3))) %>%
  # material=diffuse(color="#754315",
  #                 gradient_color="darkred",
  #                 gradient_point_start = c(0,0,0),
  #                 gradient_point_end = c(0,8,0)))) %>%
  add_object(sphere(y=200,radius=80,material=light()))


estimated_deaths_scene <-
  xz_rect(y=-0.01,zwidth=500,xwidth=500,
          material = diffuse(color="grey3")) %>%
  add_object(extruded_polygon(this_sf, center=TRUE,z=-5,scale=c(1/100000,1,1/100000),
                              data_column_top = "deaths_per_day_est_cum_per_person",
                              scale_data = 1/max(states_geom_and_data$deaths_per_day_est_cum_per_person,na.rm=TRUE)*30,
                              material=glossy(
                                color="#754315",
                                # gradient_color = "#fcde9c",
                                # color="#7c1d6f",
                                gradient_color="darkred",
                                gradient_point_start = c(0,0,0),
                                gradient_point_end = c(0,7,0),
                                reflectance=0,
                                gloss = 0.3))) %>%
  # material=diffuse(color="#754315",
  #                 gradient_color="darkred",
  #                 gradient_point_start = c(0,0,0),
  #                 gradient_point_end = c(0,8,0)))) %>%
  add_object(sphere(y=200,radius=80,material=light()))






render_scene(scene = reported_deaths_scene, parallel=TRUE,lookfrom = c(0,100,-30),lookat=c(0,0,0),samples=1000,
             filename=paste0("figures/frames/rayrender_deaths_reported",date_ind,".png"),
             fov=0,clamp_value=10,ortho_dimensions = c(55,55),width=1000,height=1000)


render_scene(scene = estimated_deaths_scene, parallel=TRUE,lookfrom = c(0,100,-30),lookat=c(0,0,0),samples=1000,
             filename=paste0("figures/frames/rayrender_deaths_estimated",date_ind,".png"),
             fov=0,clamp_value=10,ortho_dimensions = c(55,55),width=1000,height=1000)

render_scene(scene = reported_deaths_scene_red, parallel=TRUE,lookfrom = c(0,100,-30),lookat=c(0,0,0),samples=1000,
             filename=paste0("figures/frames/rayrender_deaths_reported_red",date_ind,".png"),
             fov=0,clamp_value=10,ortho_dimensions = c(55,55),width=1000,height=1000)


# }





