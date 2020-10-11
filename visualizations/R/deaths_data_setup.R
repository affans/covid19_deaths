# Run this first to preprocess data

library(tidyverse);library(lubridate)

# load in Deaths data and annotate ---- 
deaths_per_state <- read.csv(file = "input_data/Death.csv",header = F)

rownames(deaths_per_state) <-  seq.Date(from = ymd("2020-01-23"),to = ymd("2020-10-03"),by = "1 day")

# dates_mat <- data.frame(date = seq.Date(from = ymd("2020-01-22"),to = ymd("2020-10-03"),by = "1 day"))

state_names <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", 
                 "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA",
                 "MD","ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE","NH",
                 "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", 
                 "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")

colnames(deaths_per_state) <- state_names

galvani_deaths_per_state <- deaths_per_state %>% 
  mutate(date = rownames(deaths_per_state)) %>%
  pivot_longer(-date,names_to="state", values_to="deaths_per_day") %>%
  mutate(date = lubridate::ymd(date))



deaths_per_state_est <- read.csv(file = "input_data/EDeath.csv",header = F)

rownames(deaths_per_state_est) <-  seq.Date(from = ymd("2020-01-23"),to = ymd("2020-10-03"),by = "1 day")

state_names <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", 
                 "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA",
                 "MD","ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE","NH",
                 "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", 
                 "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")

colnames(deaths_per_state_est) <- state_names


galvani_deaths_per_state_est <- deaths_per_state_est %>% 
  mutate(date = rownames(deaths_per_state)) %>%
  pivot_longer(-date,names_to="state", values_to="deaths_per_day_est") %>%
  mutate(date = lubridate::ymd(date))


deaths_joined <- left_join(galvani_deaths_per_state,galvani_deaths_per_state_est,
                           by = c("date", "state"))


deaths_joined <- deaths_joined %>% 
  group_by(state) %>%
  mutate(deaths_per_day_cum = cumsum(deaths_per_day),
         deaths_per_day_est_cum = cumsum(deaths_per_day_est)) %>%
  ungroup()



deaths_country <- deaths_joined %>%
  group_by(date) %>%
  summarize(deaths_per_day_cum = sum(deaths_per_day_cum),
            deaths_per_day_est_cum = sum(deaths_per_day_est_cum))




# 
# deaths_country_for_ms <- deaths_joined %>%
#   # filter(!state %in% c("HI","AK","WY","VT")) %>%
#   group_by(state) %>%
#   mutate(deaths_per_day_cum = cumsum(deaths_per_day),
#          deaths_per_day_est_cum = cumsum(deaths_per_day_est)) %>%
#   ungroup() %>%
#   group_by(date) %>%
#   summarize(deaths_per_day_cum = sum(deaths_per_day_cum),
#             deaths_per_day_est_cum = sum(deaths_per_day_est_cum))



# checking total deaths
tail(deaths_country$deaths_per_day_cum)
tail(deaths_country$deaths_per_day_est_cum)


# ran tidycensus::census_api_key() first with key sent to inbox.
# More info: https://walker-data.com/tidycensus/articles/basic-usage.html

## Ran this once and commented out so I do not have to keep bothering
## the API
# state_pop <-
#   tidycensus::get_acs(geography = "state",
#           variables = "B01003_001",
#           year = 2018,
#           geometry = TRUE)
# save(state_pop,file = "output_data/state_pop.RData")
load("output_data/state_pop.RData")


states_sf <- urbnmapr::get_urbn_map("states", sf = TRUE)


state_pop_df <- as.data.frame(state_pop)
state_pop_df <- state_pop_df %>% 
  dplyr::select(NAME,estimate)

states_sf <- left_join(states_sf,state_pop_df,by = c("state_name" = "NAME"))



states_geom_and_data <- left_join(states_sf,deaths_joined,
                                  by = c("state_abbv"="state"))

# View(states_geom_and_data)


colnames(states_geom_and_data)

states_geom_and_data <- states_geom_and_data %>% 
  mutate(deaths_per_day_per_person = deaths_per_day / estimate, 
         deaths_per_day_est_per_person = deaths_per_day_est / estimate, 
         deaths_per_day_cum_per_person = deaths_per_day_cum / estimate,
         deaths_per_day_est_cum_per_person = deaths_per_day_est_cum / estimate)


save(states_geom_and_data,file = "output_data/states_geom_and_data.RData")


deaths_country_m <- deaths_country %>%
  pivot_longer(-date, names_to = "Death type", values_to = "Cumulative deaths")


save(deaths_country_m, file = "output_data/country_deaths.RData")


# 
# deaths_country_m_for_ms <- deaths_country_for_ms %>%
#   pivot_longer(-date, names_to = "Death type", values_to = "Cumulative deaths")
# 
# 
# save(deaths_country_m_for_ms, file = "output_data/country_deaths_sep2020_for_ms.RData")



