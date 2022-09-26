library(tidyverse)
library(sf)

listings <- read_csv("stories/airbnb/data/listings.csv")

lists <- listings |>
  filter(number_of_reviews > 2, minimum_nights <= 7, id != 19249217)

# plot them?
neighbours <- read_sf("stories/airbnb/data/neighbourhoods.geojson")


ggplot(neighbours) +
  geom_sf() +
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4)

# train stations
train <- read_csv("stories/airbnb/data/StationEntrances2020_v4.csv")

ggplot(neighbours) +
  geom_sf() +
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4) +
  geom_point(data=train, aes(x=LONG, y=LAT), col='red')

# map stuff
map <- read_sf("stories/airbnb/data/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp")

nsw_boundary <- map |> filter(STE_NAME21 == "New South Wales")

ggplot(nsw_boundary) +
  geom_sf() +
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4) +
  geom_point(data=train, aes(x=LONG, y=LAT), col='red')

# Now we could compute distance to nearest train station and distance to sea (harbour)
# as some covariates for price?

