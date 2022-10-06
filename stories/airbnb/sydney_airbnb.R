library(tidyverse)
library(sf)

list_price <- read_csv("stories/airbnb/data/listings.csv") |>
  select(clean_price = price, id)
listings <- read_csv("stories/airbnb/data/listings.csv.gz") |>
  left_join(list_price)

lists <- listings |>
  filter(number_of_reviews > 2, minimum_nights <= 7, id != 19249217)

listings <- lists |>
  st_as_sf(coords=c('longitude', 'latitude'), crs = st_crs(4326)) |>
  st_transform(crs = st_crs(7856))

# Border for our analyses (around Sydney CBD)
bb <- st_bbox(c(xmin=328000,
          xmax=350000,
          ymin=6242000,
          ymax=6263000)) |>
  st_set_crs(st_crs(7856)) |>
  st_as_sfc()

# AU map for border
map <- read_sf("stories/airbnb/data/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp")

au <- map |>
  st_transform(crs=st_crs(7856)) |>
  mutate(geometry = st_make_valid(geometry)) |>
  summarise(geometry = st_union(geometry))

border <- au |>
  summarise(geometry = st_union(geometry)) |>
  st_intersection(bb)

plot(border)

# OK, listings should be the smaller region
listings <- listings |> st_intersection(bb)

# train stations
train <- read_csv("stories/airbnb/data/StationEntrances2020_v4.csv")

train_sf <- train |>
  st_as_sf(coords = c('LONG', 'LAT'), crs=st_crs(4326)) |>
  st_transform(crs=st_crs(7856))

# OK, see if we can get the NSW coastline out. Idea is to intersect with
# an expanded NSW border...
buff <- bb |>
  st_geometry() |>
  st_buffer(dist = 3000)

au_line <-  au |> st_geometry() |>
  st_transform(crs = st_crs(buff)) |>
  st_cast(to = "MULTILINESTRING")

coast <- au_line |>
  st_intersection(buff)

plot(coast)

# generate a grid and intersect the grid with the border
syd_grid <- border |>
  st_bbox() |>
  st_as_sfc() |>
  st_make_grid(
    n = 100,
    what = "centers"
  ) |>
  st_intersection(border)

# OK, now compute how far each point is to the coast.
train_locs <- train_sf |> st_geometry() |> st_combine()

gridded <- tibble(dist_to_coast =
                    st_distance(syd_grid, coast),
                  dist_to_rail =
                    st_distance(syd_grid, train_locs)) |>
  st_set_geometry(syd_grid) |>
  mutate(across(starts_with('dist'), as.numeric))

# OK, now generate our dataset with the same measures

# OK, now let's do some modelling...
house_apart <- listings |> filter(room_type == "Entire home/apt") |>
  mutate(bedrooms = case_when(is.na(bedrooms) ~ floor((accommodates + 1)/2),
         TRUE ~ bedrooms)) |>
  filter(bedrooms < 6) |>
  mutate(dist_to_coast = st_distance(geometry, coast),
         dist_to_rail  = st_distance(geometry, train_locs)) |>
  mutate(across(starts_with('dist'), as.numeric))

ggplot(house_apart) +
  geom_point(aes(x=bedrooms, y=clean_price)) +
  scale_y_log10()

ggplot(house_apart) +
  geom_point(aes(x=dist_to_coast, y=clean_price)) +
  scale_x_log10() +
  scale_y_log10()

ggplot(house_apart |> mutate(coast = cut(dist_to_coast, breaks = c(0, 100, 1000, 10000, 100000)))) +
  geom_point(aes(x=dist_to_rail, y=clean_price)) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(vars(bedrooms)) +
  geom_smooth(aes(x=dist_to_rail, y=clean_price))

#m <- lm(log(clean_price) ~ bedrooms + log(dist_to_coast) +
#          log(dist_to_rail), data=house_apart)
m <- lm(log(clean_price) ~ bedrooms + log(dist_to_coast), data=house_apart)
summary(m)
visreg::visreg(m)

# OK, just bedrooms and distance to coast...

# plot them?
if (0) {
ggplot(neighbours) +
  geom_sf() +
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4)

ggplot(border) +
  geom_sf()+
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4)

ggplot(neighbours) +
  geom_sf() +
  geom_point(data=lists |> arrange(price), aes(x=longitude, y=latitude, col=price), alpha=0.4) +
  geom_point(data=train, aes(x=LONG, y=LAT), col='red')

plot(coast)
plot(syd_grid, add=TRUE)

gridded |>
  ggplot() +
  geom_sf(data=coast) +
  geom_sf(aes(col=dist_to_rail, alpha=log(dist_to_coast)))
}


# OK, now do some kriging stuff with listings and gridded for price

library(tidyverse)
library(sf)
library(gstat)
library(stars)

# start by plotting the data. Hmm, we may need to downsample!
set.seed(5)
house_apart_sample <- house_apart |> slice_sample(n=500)

# check locations aren't duplicated
house_apart_sample |>
  st_geometry() |>
  as_Spatial() |>
  sp::zerodist()

ggplot() +
  geom_sf(data=coast) +
  geom_sf(data=house_apart_sample, aes(size=clean_price), shape=21,
          col='steelblue')
  
# zoom in on the cbd
ggplot() +
  geom_sf(data=coast) +
  geom_sf(data=house_apart_sample, aes(size=clean_price), shape=21,
          col='steelblue') +
  xlim(c(320000, 340000)) + ylim(c(6240000, 6260000))

# it stands to reason that the price correlates with number of
# bedrooms and distance to the coast.
ggplot(house_apart_sample) +
  geom_point(aes(x=dist_to_coast, y=clean_price)) +
  scale_y_log10() +
  scale_x_log10()

ggplot(house_apart_sample) +
  geom_boxplot(aes(x=factor(bedrooms), y=clean_price)) +
  scale_y_log10()

ggplot(house_apart) +
  geom_sf(aes(col = log(clean_price/bedrooms)), alpha=0.5)

# Hmm, not much going on. Shall we fuck about with the data a bit to add some correlation?
house_apart_mm <- house_apart |>
  mutate(clean_price = clean_price * exp(-0.3*log(dist_to_coast)) / bedrooms)

# Hmm, way better would be some sort of covariance smoothing. How to do that?
# Add some correlation around train stations?
set.seed(4)
house_apart_mm <- house_apart |>
  filter(bedrooms == 2) |>
  mutate(clean_price = signif(clean_price * exp(-0.3*log(dist_to_rail) - 0.2*log(dist_to_coast)) * 40, 2)) |>
  filter(!duplicated(geometry)) #|>
#  slice_sample(n=1000)

house_apart_mm |>
  select(id, name, price = clean_price, dist_to_coast) |>
  mutate(x = st_coordinates(geometry) |> as.data.frame()) |>
  unnest(x) |>
  st_set_geometry(NULL) |>
  write_csv("data/airbnb/sydney.csv")

gridded |> select(dist_to_coast) |>
  mutate(x = st_coordinates(geometry) |> as.data.frame()) |>
  unnest(x) |>
  st_set_geometry(NULL) |>
  write_csv("data/airbnb/sydney_grid.csv")

coast |>
  write_sf("data/airbnb/sydney_coast.geojson")

# OK, read the data in
sydney <- read_csv("data/airbnb/sydney.csv") |>
  st_as_sf(coords = c('X', 'Y'), crs=7856)

sydney_grid <- read_csv("data/airbnb/sydney_grid.csv") |>
  st_as_sf(coords = c('X', 'Y'), crs=7856)

coast <- read_sf("data/airbnb/sydney_coast.geojson")

sydney <- read_csv("https://www.massey.ac.nz/~jcmarsha/233214/data/airbnb/sydney.csv") |>
  st_as_sf(coords = c('X', 'Y'), crs=7856)

sydney_grid <- read_csv("https://www.massey.ac.nz/~jcmarsha/233214/data/airbnb/sydney_grid.csv") |>
  st_as_sf(coords = c('X', 'Y'), crs=7856)

coast <- read_sf("https://www.massey.ac.nz/~jcmarsha/233214/data/airbnb/sydney_coast.geojson")

ggplot(sydney) +
  geom_sf(aes(col = price), alpha=0.5) +
  geom_sf(data=coast) +
  scale_colour_continuous(trans=scales::log10_trans())

# OK, now do our investigation etc

lm(log(price) ~ log(dist_to_coast), data=sydney) |> summary()

# try some kriging
library(gstat)
const.vgm <- variogram(log(price) ~ 1, data=sydney)
plot(const.vgm)

const.mod <- fit.variogram(const.vgm, model=vgm(0.5, "Sph", 5000, 0.1))
plot(const.vgm, const.mod)

const.krig <- krige(log(price) ~ 1, locations=sydney,
                    newdata=sydney_grid, model=const.mod)

ggplot(const.krig) +
  geom_sf(aes(col=exp(var1.pred))) +
  geom_sf(data=coast) +
  scale_color_viridis_c(trans = 'log10', limits=c(80, 3000)) +
  labs(col = 'Price per night')

ggplot(const.krig) +
  geom_sf(aes(col=var1.var)) +
  geom_sf(data=coast) +
  scale_color_viridis_c(option="B")

#  geom_sf(data=house_apart_mm, mapping=aes(size=log(clean_price)),
 #         shape=21, alpha=0.1)

covar.vgm <- variogram(log(price) ~ log(dist_to_coast), data=sydney)
plot(covar.vgm)

covar.mod <- fit.variogram(covar.vgm, model=vgm(0.5, "Sph", 5000, 0.1))
plot(covar.vgm, covar.mod)

covar.krig <- krige(log(price) ~ log(dist_to_coast), locations=sydney,
                    newdata=sydney_grid, model=covar.mod)

ggplot(covar.krig) +
  geom_sf(aes(col=exp(var1.pred))) +
  geom_sf(data=coast) +
  scale_color_viridis_c(trans = 'log10', limits=c(80,3000)) +
  labs(col = 'Price per night')
#  geom_sf(data=house_apart_mm, mapping=aes(size=log(clean_price)),
#          shape=21, alpha=0.3)

ggplot(covar.krig) +
  geom_sf(aes(col=var1.var)) +
  geom_sf(data=coast) +
  scale_color_viridis_c(option="B")


ggplot(house_apart_mm) +
  geom_point(aes(x=dist_to_coast, y=clean_price)) +
  scale_x_log10() +
  scale_y_log10()
#+
#  geom_sf(data=house_apart_mm, mapping=aes(size=log(clean_price)),
#          shape=21, alpha=0.1)

# add bedrooms
bedroom.vgm <- variogram(log(clean_price) ~ log(dist_to_coast) + bedrooms, data=house_apart_sample, cutoff=500)
plot(bedroom.vgm)

bedroom.mod <- fit.variogram(bedroom.vgm, model=vgm(0.3, "Exp", 20, 0.1), debug.level = 4)
plot(bedroom.vgm, bedroom.mod)

grid_brm <- gridded |> mutate(bedrooms = 2)
bedroom.krig <- krige(log(clean_price) ~ log(dist_to_coast) + bedrooms, locations=house_apart_sample,
                    newdata=grid_brm, model=bedroom.mod)
bedroom.krig

ggplot(bedroom.krig) +
  geom_sf(aes(col=var1.pred)) +
  geom_sf(data=coast) +
  geom_sf(data=house_apart_sample, mapping=aes(size=log(clean_price)),
          shape=21, alpha=0.3)

ggplot(bedroom.krig) +
  geom_sf(aes(col=var1.var)) +
  geom_sf(data=coast) +
  geom_sf(data=house_apart_sample, mapping=aes(size=log(clean_price)),
          shape=21, alpha=0.3)

d <- st_distance(house_apart_sample)
gamma <- dist(log(house_apart_sample$clean_price)) |> as.matrix()
gamma_vs_d <- data.frame(gamma = as.numeric(gamma)/2, dist = as.numeric(d)) |>
  filter(dist < 10000)
ggplot(gamma_vs_d) +
  geom_point(aes(x=dist, y=gamma), alpha=0.3) +
  theme_minimal()

bands <- gamma_vs_d |> mutate(band = cut(dist, breaks=15)) |>
  group_by(band) |> summarise(dist=mean(dist), gamma=mean(gamma))

ggplot(gamma_vs_d) +
  geom_point(aes(x=dist, y=gamma), alpha=0.1) +
  geom_point(data=bands, mapping=aes(x=dist, y=gamma), size=4, col='dark red') +
  theme_minimal()

foo <- variogram(log(clean_price) ~ 1, data=house_apart_sample |> st_transform(crs=st_crs(7856)), cutoff=1500)
foo.vgm <- fit.variogram(foo, model=vgm(0.6, "Sph", 500, 0.1))
plot(foo, foo.vgm)

foo <- variogram(log(clean_price) ~ log(dist_to_coast), data=house_apart_sample |> st_transform(crs=st_crs(7856)), cutoff=1500)
foo.vgm <- fit.variogram(foo, model=vgm(0.6, "Sph", 500, 0.1))
plot(foo, foo.vgm)

foo <- variogram(log(clean_price) ~ log(dist_to_coast) + bedrooms, data=house_apart_sample |> st_transform(crs=st_crs(7856)), cutoff=1500)
foo.vgm <- fit.variogram(foo, model=vgm(0.6, "Sph", 500, 0.1))
plot(foo, foo.vgm)

ggplot() +
  geom_sf(data=house_apart_sample |> st_transform(crs=st_crs(7856)), mapping=aes(size=log(clean_price)),
          shape=21, alpha=0.3) +
  geom_sf(data=coast)

  