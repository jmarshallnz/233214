library(tidyverse)
#library(raster)
library(sf)
library(gstat)

####  Mass  analysis ####

if (0) {
  ### Path to covariates
  base_path <- "../winTor_aux/data/"
  ## Co-variates
  env.names <- c("NA_dem", "NA_northing", "NA_nFrostyDays",
                 "NA_nonGrowingDays", "NA_nDaysFreeze", "NA_OG1k")
  env.stk <- raster::subset(raster::stack(list.files(base_path, pattern = "NA_*", full.names = T)), env.names)
  
  env.stk
  test <- aggregate(env.stk, fact=100)
  ## ammending to have co-variate data
  mass <- read.csv("data/massDataReferenced.csv")
  coordinates(mass) <- ~ Long + Lat
  proj4string(mass) <- proj4string(env.stk)
  mass@data <- as.data.frame(cbind(mass, raster::extract(env.stk, mass)))
  
  ## try to determine if there are any points that'd fall within the
  ## same raster block (over inflate spatial AC)
  mass.cellList <- cellFromXY(env.stk$NA_dem, mass)
  anyDuplicated(mass.cellList) 
  which(duplicated(mass.cellList)) ## 22 25 41 42
  ##examine
  which(mass.cellList==mass.cellList[[22]]) ## 21 & 22
  which(mass.cellList==mass.cellList[[25]]) ## 24 & 25
  which(mass.cellList==mass.cellList[[41]]) ## 40 & 41
  which(mass.cellList==mass.cellList[[42]]) ## 37 & 42
  ## average the 2 and create new entries
  mass@data[c(21,22),] ## average
  mass@data$avgMass[21] <- (mass@data$avgMass[21] + mass@data$avgMass[22])/2 
  mass@data[c(24,25),] ## average
  mass@data$avgMass[24] <- (mass@data$avgMass[24] + mass@data$avgMass[25])/2 
  mass@data[c(40,41),]
  mass@data$avgMass[41] <- (mass@data$avgMass[40] + mass@data$avgMass[41])/2 
  mass@data[c(37,42),]
  mass@data$avgMass[37] <- (mass@data$avgMass[37] + mass@data$avgMass[42])/2 
  
  nrow(mass)
  mass <- mass[-c(22,25,40,42),]
  nrow(mass)
  ## checked and fixed

  # OK, write out stuff here? What type of object should we use. Why not sf?
  library(sf)
  # remove the 14.5 value
  mass_sf <- st_as_sf(mass)
  fs::file_delete("data/bat_mass.csv")
  write_sf(mass_sf %>% filter(avgMass != 14.5), "data/bat_mass.csv")

  # Can we downsample our raster images?

  library(stars)
  test.stars <- st_as_stars(test)
  test.sf <- st_as_sf(test.stars)

  test.sf %>% select(elevation = NA_dem,
                     northing = NA_northing,
                     frost_days = NA_nFrostyDays,
                     freeze_days = NA_nDaysFreeze,
                     growing_days = NA_nonGrowingDays) %>%
    mutate(growing_days = 365 - growing_days) %>%
    mutate(across(-geometry, ~round(., digits=2))) %>%
    write_sf("data/NA_info.sqlite")
  
  # downsample our boundary file
  NA_boundary <- read_sf("data/nam_high_res.shp")
  NA_boundary %>% mutate(Country = ifelse(nchar(STATE) == 2, "UNITED STATES", STATE),
                         Country = ifelse(Country %in% c("UNITED STATES", "MEXICO"), Country, "CANADA")) %>%
    group_by(Country) %>%
    summarise(geometry = st_union(geometry)) %>%
    st_simplify(dTolerance=2000) %>%
    st_transform(proj4string(env.stk)) -> countries
  write_sf(countries, "data/NA_boundary.sqlite")
}

# Read in NA covariate data and boundary
NA_info <- read_sf("stories/bats/data/NA_info.sqlite")
NA_boundary <- read_sf("stories/bats/data/NA_boundary.sqlite")

# Read in MASS data
mass.csv = read_csv("stories/bats/data/bat_mass.csv")
mass = st_as_sf(mass.csv, coords=c("Long","Lat"), crs=st_crs(NA_info))

proj4_map <- "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +datum=WGS84 +units=m"

# plot our area
ggplot() +
  geom_sf(data=NA_boundary) +
  geom_sf(data=mass, aes(size=avgMass), shape=21) +
  coord_sf(crs = proj4_map)

# NICE! Now we are getting somewhere.
ggplot() +
  geom_sf(data=NA_info, aes(fill=elevation), lwd=0) +
  geom_sf(data=NA_boundary, fill=NA) +
  geom_sf(data=mass, aes(size=avgMass), shape=21) +
  coord_sf(crs = proj4_map)

# OK, now let's try the kriging bullshit

# To do kriging, we really need a set of points that we
# can interpolate to. If we give the routines polygons,
# it's just going to discretise them anyway. Given we have
# a nice grid already, let's use the centroid of each tile.
# We'll keep the tile geometry in the 'tile' column:

NA_points <- NA_info %>%
  mutate(tile = GEOMETRY,
         GEOMETRY = st_centroid(tile))

st_crs(NA_points)
st_crs(mass)

# We'll start with inverse distance weighting
mass.idw <- idw(avgMass ~ 1, locations=mass, newdata=NA_points)

# Join back to our NA_points data so we get the tile column for plotting
NA.idw <- mass.idw %>%
  st_join(NA_points)

# Plot
ggplot() +
  geom_sf(data=NA_boundary, fill=NA) +
  geom_sf(data=NA.idw, aes(geometry=tile, fill=var1.pred), col=NA) +
  geom_sf(data=mass, aes(size=avgMass), shape=21) +
  coord_sf(crs = proj4_map) +
  scale_fill_viridis_c()

# The alternate is kriging: we fit a model to the variogram (spatial variation)
const.vgm <- variogram(avgMass ~ 1, data=mass, cutoff=5000)
const.mod <- fit.variogram(const.vgm, model=vgm(1, "Exp", 3000, 1))
plot(const.vgm, const.mod)

# A linear variogram fits here. Let's run with it.

# This isn't the best fit - we don't have huge numbers of points to play with
# but let's run with it.

kout <- krige(avgMass ~ 1, locations=mass, newdata=NA_points, model=const.mod)
kout <- kout %>% st_join(NA_points)

ggplot() +
  geom_sf(data=kout, aes(geometry=tile, fill=var1.pred), col=NA) +
  geom_sf(data=NA_boundary, fill=NA) +
  geom_sf(data=mass, aes(size=avgMass), shape=21) +
  coord_sf(crs = proj4_map) +
  scale_fill_viridis_c()

# OK, how about modelling?
covar.vgm <- variogram(avgMass ~ NA_northing + NA_nDaysFreeze, locations=mass, cutoff=5000)
covar.mod <- fit.variogram(covar.vgm, model=vgm(1, "Exp", 3000, 1))
plot(covar.vgm, covar.mod)

kout <- krige(avgMass ~ NA_northing + NA_nDaysFreeze, locations=mass, newdata=NA_points, model=covar.mod)
kout <- kout %>% st_join(NA_points)

ggplot() +
  geom_sf(data=kout, aes(geometry=tile, fill=var1.pred), col=NA) +
  geom_sf(data=countries, fill=NA) +
  geom_sf(data=mass_sf, aes(size=avgMass), shape=21) +
  coord_sf(crs = proj4_map) +
  scale_fill_viridis_c()


# Example linear model for trend:
mod <- lm(avgMass ~ NA_northing + NA_nDaysFreeze, data=mass_sf)
mod %>% summary()

fit <- mod %>% broom::augment(newdata=mass_sf) %>%
  st_as_sf()

# What we notice is spatial clustering of similar colours.
# This isn't really what we want.
ggplot() +
  geom_sf(data=countries, fill=NA) +
  geom_sf(data=fit, aes(col=.resid))

# We could try and fit a variogram to the underlying data.
# and then to the data with various covariates
mass_vgm <- variogram(avgMass ~ NA_northing + NA_nDaysFreeze, mass_sf)
mass_mod <- fit.variogram(mass_vgm, model=vgm(1, "Exp", 1000, 1))
plot(mass_vgm, mass_mod)

# This isn't the best fit - we don't have huge numbers of points to play with
# but let's run with it.

kout <- krige(avgMass ~ NA_northing + NA_nDaysFreeze, mass_sf, newdata=test2, model=mass_mod)
kout <- kout %>% st_join(test2)

ggplot() +
  geom_sf(data=kout, aes(geometry=tile, fill=var1.pred), col=NA) +
  geom_sf(data=countries, fill=NA) +
  geom_sf(data=mass_sf, aes(size=avgMass), shape=21) +
  coord_sf(crs = proj4_map) +
  scale_fill_viridis_c()

# uncertainty
ggplot() +
  geom_sf(data=kout, aes(geometry=tile, fill=var1.var), col=NA) +
  geom_sf(data=countries, fill=NA) +
  geom_sf(data=mass_sf, aes(size=avgMass), shape=21) +
  coord_sf(crs = proj4_map) +
  scale_fill_viridis_c()

# Ok, this is sort-of working!



vgm()
show.vgms()

?stars
(test)

crs(NA_boundary)

points(mass)

## re-create distance matrix since an obs was removed
knn <- knearneigh(mass, k = 6) ## k set as the sqrt of nrow(mass) rounded
nneighbor <- knn2nb(knn)
## plot to check things out
plot(mass)
plot(nneighbor, coords = coordinates(mass), col = "red", add = T)
#convert to weights list
nn.list <- nb2listw(nneighbor, style = "W")
## test
moran.mc(mass$avgMass, nn.list, nsim = 999)

## test residuals for spatial correlation
moran.mc(mass$RESID, listw = nn.list, nsim = 999) 
## lagrange multiplier - another test for spatial effects
lm.LMtests(mass.mods[[10]], listw = nn.list, test = "all")

## significant auto-correlation among the residuals
#### Spatial GLM ####
## Methodology adapted from: https://cran.r-project.org/web/packages/glmmfields/vignettes/spatial-glms.html
library(glmmfields)
options(mc.cores = parallel::detectCores())  
system.time(
mass_spatial <- glmmfields(mass.mod,
                          data = mass@data,
                          lat = "Lat", lon = "Long",
                          nknots = 6, iter = 10000, chains = 5,
                          prior_intercept = student_t(3, 0, 10), 
                          prior_beta = student_t(3, 0, 3),
                          prior_sigma = half_t(3, 0, 3),
                          prior_gp_theta = half_t(3, 0, 10),
                          prior_gp_sigma = half_t(3, 0, 3),
                          seed = 123,
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 15)# passed to rstan::sampling()
))

## Summary
summary(mass_spatial)

plot(mass_spatial, type = "spatial-residual", link = TRUE) +
  geom_point(size = 3)
## looks much better
plot(mass_spatial, type = "residual-vs-fitted")
plot(mass_spatial, type = "prediction", link = FALSE) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)

#### Predictions and confidence bounds ####
## create layers for lat and long
env.stk$Long <- xFromCell(env.stk[[1]], cell = 1:ncell(env.stk))
env.stk$Lat <- yFromCell(env.stk[[1]],  cell = 1:ncell(env.stk))

## our prediction function - returns matrix of 3 columns
conffun <- function(model, data = NULL, iter='all') {
  v <- predict(object = model,
               newdata = data,
               type = "response",
               interval = "prediction",
               iter = iter)
  return(as.matrix(v))  
}

## change raster options to have way lower chunksize to help
## fit into memory. Can maybe increase this?
raster::rasterOptions(chunksize=1e5)
rasterOptions(memfrac = .3); rasterOptions(maxmemory = 1e+08)
library(cluster)
library(parallel)


## Use for prediction to full
beginCluster(2) ##can hand way too long depending on the number of cores
system.time({
  r.prob.Cluster<-clusterR(env.stk, fun=predict, args=list(model=mass_spatial,
                                                         fun=conffun,
                                                         progress='text',
                                                         index=1:3,
                                                         iter=1000))
})
endCluster() #delete the cluster


names(r.prob.Cluster) <- c("p", "lwr", "upr")
writeRaster(r.prob.Cluster,
            filename = file.path(win.res, "mass.tif"),
            format = "GTiff",
            bylayer = T,
            suffix = "names",
            overwrite = T)

#### Clean up script items ####
env.post <- ls()
to.remove <- env.post[env.post %!in% env.prior]
rm(list=to.remove); rm(env.post, to.remove)
