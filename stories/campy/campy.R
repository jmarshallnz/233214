library(tidyverse)
library(sf)
library(sparr)

# Now let's see if we can do something with campy cases maybe?
if (0) {
  shp <- read_sf("stories/campy/data/midcentral_phu.shp")
  owin <- st_union(shp) # boundary file
  write_sf(owin, "stores/campy/data/midcentral_phu_boundary.shp")
  
  cases <- read_csv("data/campy/data/campy_midcentral_phu.csv")
  cases %>% filter(ReportDate >= '2005-01-01' &
                     ReportDate <= '2007-12-31') %>% select(x, y) %>%
    write_csv("data/campy/midcentral_phu_cases.csv")
  
  set.seed(2021)
  shp <- read_sf("stories/campy/data/midcentral_phu.shp")
  controls <- shp %>% slice_sample(n=1000, weight_by=POPULATION, replace = TRUE)
  
  # and generate their centroid?
  controls_cent <- controls %>% st_centroid() %>%
    st_coordinates() %>% as.data.frame() %>% set_names(nm=c("x","y"))
  controls_cent %>% write_csv("stories/campy/data/midcentral_phu_controls.csv")
}

# PALMY WORKS BETTER!
if (0) {
  shp <- read_sf("stories/campy/data/midcentral_phu.shp")
  plot(shp["SDI"], xlim = c(2725000,2740000), ylim=c(6085000, 6097000))

  shp_palmy <- st_crop(shp, xmin=2725000, ymin=6085000, xmax=2740000, ymax=6097000)
#  write_sf(shp_palmy %>% st_union(), "data/campy/palmy_boundary.shp")
  
  # Ok, find our cases inside palmy
  cases <- read_csv("stories/campy/data/campy_midcentral_phu.csv")
  cases_palmy <- cases %>% filter(ReportDate >= '2005-01-01' &
                                    ReportDate <= '2007-12-31') %>%
    semi_join(shp_palmy %>% dplyr::select(Meshblock06 = MB06)) %>%
    dplyr::select(x, y)
  
  cases.df <- as.ppp(cases, W = shp_palmy %>% st_union()) %>%
    as.data.frame()
  cases.df %>% write_csv("data/campy/palmy_cases.csv")
  
  # noice
  set.seed(2021)
  controls <- shp_palmy %>% slice_sample(n=1000, weight_by=POPULATION, replace = TRUE)
  
  # and generate their centroid?
  controls_cent <- controls %>% st_sample(rep(1, nrow(controls))) %>%
    st_coordinates() %>% as.data.frame() %>% set_names(nm=c("x","y"))
  
  controls_cent %>% write_csv("data/campy/palmy_controls.csv")
}

# OK, now we have x and y already
cases <- read_csv("data/campy/palmy_cases.csv")
controls <- read_csv("data/campy/palmy_controls.csv")

palmy.owin <- owin(xrange=c(2725000, 2740000), yrange=c(6085000, 6097000))
cases.ppp <- as.ppp(cases, W = palmy.owin)
controls.ppp <- as.ppp(controls, W = palmy.owin)

# Plot using spatstat
plot(palmy.owin, main="Palmerston North")
points(controls.ppp, col='black')
points(cases.ppp, col='red')

# Plot using ggplot
points <- bind_rows(control=controls, case=cases, .id = 'type')
ggplot(points) +
  geom_point(mapping=aes(x=x, y=y, col=type), alpha=0.5) +
  coord_fixed() +
  theme_void()

# Risk
rr <- risk(cases.ppp, controls.ppp, adapt=TRUE, pilot.symmetry='pooled')
plot(rr)

# Better colours. Most of the data is in -1 to 0.5. Make it symmetric
col_breaks = c(-2,seq(-1, 1, by=0.1),2)
col_breaks
num_cols = length(col_breaks) - 1
num_cols
cols = hcl.colors(num_cols, "blue-red")
plot(rr, breaks=col_breaks, col = cols)
# Ok, so places of increased risk on the eastern side of Palmerston North (generally wealthier areas!)

# Let's try modelling cases. We could use the controls, alongside perhaps a social deprivation surface

# What bandwidth was used above?
rr$f$h0

# Ok, let's use h0=800
control.bd <- bivariate.density(controls.ppp, h0 = 800, adapt=TRUE, trim=10)
plot(control.bd)

# pull out the surface and fit a model for cases in terms of the control distribution
control.im <- control.bd$z
out <- ppm(cases.ppp ~ log(control.im))
out

# As expected, where people are matter for where cases of campylobacteriosis are.

# Now social deprivation.
# raster for SDI
if (0) {
  library(stars)
  template <- st_as_stars(control.im)
  shp <- read_sf("stories/campy/data/midcentral_phu.shp")
  shp_palmy <- st_crop(shp, xmin=2725000, ymin=6085000, xmax=2740000, ymax=6097000)
  test.st = st_rasterize(shp_palmy["SDI"], template = template)
  test.rast <- as(test.st, "Raster")
  writeRaster(test.rast, "data/campy/palmy_SDI.tif", overwrite=TRUE)
}

sdi.rast = raster("data/campy/palmy_SDI.tif")
sdi.im <- as.im(sdi.rast)
plot(sdi.im)

# maybe that is enough. Yes!
out <- ppm(cases.ppp ~ log(control.im) + log(sdi.im))
out

# Yay!