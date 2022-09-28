library(tidyverse)
library(sf)

listings <- read_csv("stories/airbnb/data/listings.csv")

lists <- listings |>
  filter(number_of_reviews > 2, minimum_nights <= 7, id != 19249217)

# plot them?
neighbours <- read_sf("stories/airbnb/data/neighbourhoods.geojson") |>
  mutate(geometry = st_make_valid(geometry))


ggplot(neighbours) +
  geom_sf() +
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4)

border <- neighbours |>
  summarise(geometry = st_union(geometry))

ggplot(border) +
  geom_sf()+
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4)


# train stations
train <- read_csv("stories/airbnb/data/StationEntrances2020_v4.csv")

train_sf <- train |>
  st_as_sf(coords = c('LONG', 'LAT'), crs=st_crs(4326))

ggplot(neighbours) +
  geom_sf() +
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4) +
  geom_point(data=train, aes(x=LONG, y=LAT), col='red')

# map stuff
map <- read_sf("stories/airbnb/data/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp")

au <- map |> mutate(geometry = st_make_valid(geometry)) |>
  summarise(geometry = st_union(geometry))

# OK, see if we can get the NSW coastline out. Idea is to intersect with
# an expanded NSW border...
buff <- border |>
  st_geometry() |>
  st_buffer(dist = 1500)
plot(buff)

au_buff <- au |> st_geometry() |>
  st_transform(crs = st_crs(buff)) |>
  st_intersection(buff)

au_line <-  au |> st_geometry() |>
  st_transform(crs = st_crs(buff)) |>
  st_cast(to = "MULTILINESTRING")

coast <- au_line |>
  st_intersection(buff)

# ok, now compute the distance from points in our grid (and our dataset)
# to points on the coast

# generate a grid
my_grid <- border |>
  st_bbox() |>
  st_as_sfc() |>
  st_make_grid(
    n = 190,
    what = "centers"
  )

plot(my_grid)
# intersect the grid with the border
syd_grid <-
  my_grid |>
  st_intersection(border)

plot(syd_grid)

# ok, now combine shit...
plot(coast)
plot(syd_grid, add=TRUE)

# OK, now compute how far each point is to the coast.
train_locs <- train_sf |> st_geometry() |> st_combine()

gridded <- tibble(dist_to_coast =
         st_distance(syd_grid, coast),
       dist_to_rail =
         st_distance(syd_grid, train_locs)) |>
  st_set_geometry(syd_grid) |>
  mutate(across(starts_with('dist'), as.numeric))

gridded |>
  ggplot() +
  geom_sf(data=coast) +
  geom_sf(aes(col=dist_to_rail, alpha=log(dist_to_coast)))

# OK, now generate our dataset with the same measures
listings <- lists |>
  st_as_sf(coords=c('longitude', 'latitude'), crs = st_crs(4326)) |>
  mutate(dist_to_coast = st_distance(geometry, coast),
         dist_to_rail  = st_distance(geometry, train_locs)) |>
  mutate(across(starts_with('dist'), as.numeric)) |>
  select(starts_with('dist'))

# OK, now do some kriging stuff with listings and gridded for price


