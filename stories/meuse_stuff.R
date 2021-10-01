library(tidyverse)
library(sf)
library(gstat)


data(meuse.all)
data(meuse)
data(meuse.grid)
meuse

meuse.sf <- st_as_sf(meuse, coords = c('x', 'y'))
meuse.grid.sf <- st_as_sf(meuse.grid, coords = c('x', 'y'))

ggplot(meuse.sf) +
  geom_sf(aes(col=zinc))

data(meuse.riv)

st_as_sf(meuse.area, coords=c('x', 'y'))

# convert river to a line string for plotting with ggplot
river <- st_linestring(meuse.riv)

ggplot(meuse.sf) +
  geom_sf(data=river) +
  geom_sf(aes(col=zinc))

# ok, now do some messing about with kriging
meuse.idw <- idw(zinc ~ 1, locations=meuse.sf, newdata=meuse.grid.sf)

# Now, we notice that areas close to the river have high zinc counts
ggplot(meuse.sf) +
  geom_point(aes(x=dist, y=zinc))

# We want to fit a model to these data, so let's linearize by taking
# a log of the zinc variable:
ggplot(meuse.sf) +
  geom_point(aes(x=dist, y=zinc)) +
  scale_y_log10()

# We can do even better by switching to sqrt(dist)
ggplot(meuse.sf) +
  geom_point(aes(x=sqrt(dist), y=zinc)) +
  scale_y_log10()

# The idea here is that perhaps zinc inflitrates the area from the river
# and the infiltration would be in proportion to sqrt(dist)?

# so perhaps modelling log(zinc) in terms of distance from river might be useful?
m1 <- lm(log(zinc) ~ sqrt(dist), data=meuse.sf)
# clearly something going on here

# Now we could check residuals:
fit <- st_as_sf(broom::augment(m1, newdata=meuse.sf))
ggplot(fit) +
  geom_sf(aes(col=.resid))

# It seems, on the face of it that these are not random in space, right?
# there is clear grouping in space, suggesting that things are not random
# so things aren't independent.

# So the modelling assumptions don't hold, and any predictions we make will
# be ignoring the spatial correlation. Let's look at those predictions
pred <- broom::augment(m1, newdata=meuse.grid.sf) %>% st_as_sf()
ggplot(pred) +
  geom_sf(aes(col=exp(.fitted)))

# The predictions from the linear model are exactly what we'd expect - it
# basically just predicts high counts near the river. This is because we've
# assumed no spatial structure at all in the model.

# We can add some spatial structure by Kriging. Essentially this models
# the spatial dependency between observations with a variogram which can
# then be applied across the whole area to help interpolate the surface


# Before we try this with the distance data, let's see how it works if
# our model was just constant:

# We're going to fit a 'spherical' variogram model.
zn.vgm <- variogram(log(zinc) ~ 1, meuse.sf) # calculates sample variogram values 
zn.fit <- fit.variogram(zn.vgm, model=vgm(1, "Sph", 1000, 1)) # fit model
plot(zn.vgm, zn.fit)

# Once we have this, we can now do Kriging:
krig.const <- krige(log(zinc) ~ 1, meuse.sf, meuse.grid.sf, model=zn.fit)
krig.const.df <- krig.const %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

# This looks way better than predicting from our linear model that
# contained more information in it.
ggplot(krig.const.df) +
  geom_sf(data=river) +
  geom_tile(aes(x=x, y=y, fill=var1.pred)) +
  geom_sf(data=meuse.sf, shape='open circle') +
  scale_fill_viridis_c()

# Uncertainty
ggplot(krig.const.df) +
  geom_sf(data=river) +
  geom_tile(aes(x=x, y=y, fill=var1.var)) +
  geom_sf(data=meuse.sf, shape='open circle') +
  scale_fill_viridis_c()

# We could also do this with dist as well:
zn.vgm <- variogram(log(zinc) ~ sqrt(dist), meuse.sf) # calculates sample variogram values 
zn.fit <- fit.variogram(zn.vgm, model=vgm(1, "Sph", 1000, 1)) # fit model
#plot(zn.vgm, zn.fit)

krig.dist <- krige(log(zinc) ~ sqrt(dist), meuse.sf, meuse.grid.sf, model=zn.fit)
krig.dist.df <- krig.dist %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

ggplot(krig.dist.df) +
  geom_sf(data=river) +
  geom_tile(aes(x=x, y=y, fill=var1.pred)) +
  geom_sf(data=meuse.sf, shape='open circle') +
  scale_fill_gradientn(colors = scales::viridis_pal()(9), limits=c(3.5,7.5), 
                       na.value = "#FDE725FF")

# uncertainty. Neat!
ggplot(krig.dist.df) +
  geom_sf(data=river) +
  geom_tile(aes(x=x, y=y, fill=var1.var)) +
  geom_sf(data=meuse.sf, shape='open circle') +
  scale_fill_viridis_c()

# OK, what about elevation? Can we add this in from meuse.alt?
dim(meuse.grid)
dim(meuse.alt)
range(meuse.alt$x)
range(meuse.grid$x)
range(meuse.alt$y)
range(meuse.grid$y)

# Seems like we might be able to. What would we need to do? Convert to a
# stars raster then align and render? Then fill in from that?
meuse.alt.sf <- st_as_sf(meuse.alt, coords=c('x', 'y'))

test <- raster(meuse.grid.sf)

test <- st_rasterize(meuse.alt.sf)
plot(test)

do <- st_as_stars(st_bbox(meuse.grid.sf), values=NA_real_)
val=st_rasterize(meuse.alt.sf, template=do)
plot(val)

# We have to fucking krig it you dumbarse!
var <- variogram(alt ~ 1, meuse.alt.sf)
mod <- fit.variogram(var, model=vgm(1, "Exp", 1000, 1)) # fit model
plot(var, mod)

# Looks OK - now krige
o <- krige(alt ~ 1, meuse.alt.sf, newdata=meuse.grid.sf, model=mod)

# OK, check this:
alt.df <- o %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

ggplot(alt.df) +
  geom_sf(data=river) +
  geom_tile(aes(x=x, y=y, fill=var1.pred)) +
  geom_sf(data=meuse.sf, shape='open circle') +
  scale_fill_viridis_c()

gridded <- meuse.grid.sf %>% st_join(o %>% rename(alt=var1.pred))

# What about fitting to our original data?
oo <- krige(alt ~ 1, meuse.alt.sf, newdata=meuse.sf, model=mod)
meuse_alt <- meuse.sf %>% st_join(oo %>% rename(alt=var1.pred))

# NOW, do this with dist and alt...

# We could also do this with dist AND alt as well:
zn.vgm <- variogram(log(zinc) ~ sqrt(dist) + alt, meuse_alt) # calculates sample variogram values 
zn.fit <- fit.variogram(zn.vgm, model=vgm(1, "Sph", 1000, 1)) # fit model
plot(zn.vgm, zn.fit)

krig.dist <- krige(log(zinc) ~ sqrt(dist) + alt, meuse_alt, gridded, model=zn.fit)
krig.dist.df <- krig.dist %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

ggplot(meuse_alt) +
  geom_point(aes(x=alt, y=log(zinc))) # not a strong relationship. What a fucking
# waste of time...

ggplot(meuse_alt) +
  geom_boxplot(aes(x=ffreq, y=log(zinc))) # not a strong relationship. What a fucking

lm(log(zinc) ~ sqrt(dist) + ffreq + soil, data=meuse_alt) %>%
  summary()

zn.vgm <- variogram(log(zinc) ~ sqrt(dist) + ffreq + soil, meuse_alt) # calculates sample variogram values 
zn.fit <- fit.variogram(zn.vgm, model=vgm(1, "Sph", 1000, 1)) # fit model
#plot(zn.vgm, zn.fit)

krig.dist <- krige(log(zinc) ~ sqrt(dist) + ffreq + soil, meuse_alt, gridded, model=zn.fit)
krig.dist.df <- krig.dist %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

ggplot(krig.dist.df) +
  geom_sf(data=river) +
  geom_tile(aes(x=x, y=y, fill=var1.pred)) +
  geom_sf(data=meuse.sf, shape='open circle') +
  scale_fill_gradientn(colors = scales::viridis_pal()(9), limits=c(3.5,7.5), 
                       na.value = "#FDE725FF")

# uncertainty. Neat!
ggplot(krig.dist.df) +
  geom_sf(data=river) +
  geom_tile(aes(x=x, y=y, fill=var1.var)) +
  geom_sf(data=meuse.sf, shape='open circle') +
  scale_fill_viridis_c()

# less uncertainty?
ggplot(meuse.grid.sf) +
  geom_sf(aes(col=ffreq)) +
  geom_sf(data=meuse.sf, aes(fill=log(zinc)), shape=21)

# NOW LET'S TRY ORIGINAL AGAIN

# Notice that this does a little strange thing right next to the river, where dist=0
