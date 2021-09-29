library(tidyverse)
library(sparr)
library(raster)
library(maptools)

# Idea: Lung cancer data from UK.
#       Campy data from NZ.

# Fit with sparr
# Fit via PPM to assess effect of rurality or socdep?

# Would need some sort of raster for ppm

# Use chorley ribble
data(chorley)

# Take a look at the data structure. It is a ppp (planer point pattern)
class(chorley)
chorley

# Note it's a marked point pattern - it has a mark showing larynx vs lung
# <detail dataset>

# There's a polygonal boundary
# Let's plot it
plot(chorley)

# It's a bit hard to see what is going on here. Let's see if we can improve
# the plot a bit by splitting into larynx and lung datasets and plotting
# them separately
plot(chorley, chars=c(2,1), cols=c('red', 'blue'))

larynx <- split(chorley)$larynx
lung   <- split(chorley)$lung

# Note that larynx and lung are no longer marked point patterns:
larynx

plot(chorley$window, main="Chorley-Ribble lung data")
points(lung, pch=2, col='blue')
points(larynx, pch=19, col='red')

# Alternatively we could do this with ggplot by converting the various
# objects to data.frames
chorley.df <- as.data.frame(chorley) %>%
  arrange(desc(marks))
chorley.win <- as.data.frame(chorley$window)

ggplot(data=chorley.df) +
  geom_point(aes(x=x, y=y, col=marks)) +
  geom_polygon(data=chorley.win, aes(x=x, y=y), fill='transparent', col='black') +
  coord_fixed()

# OK, now let's look at the distribution of larynx cases by using
# a kernel density estimate. We could do this with spatstat like was
# done in Week 10?
larynx.dens <- density(larynx)
plot(larynx.dens)
attr(larynx.dens, "sigma") # smoothing amount

# try a bit less smoothing
larynx.dens <- density(larynx, sigma=1)
plot(larynx.dens)

#  let's look at lungs as well with same smoothing amount
lung.dens <- density(lung, sigma=1)
plot(lung.dens)
points(lung)

# Note that the points are highly clustered. In some areas we'd want
# a very fine bandwidth (eg. in the built up areas of Bamber bridge (top blob),
# Leyland and Chorley), while in the more rural areas we'll want a
# large bandwidth.

# The sparr package can do that. We'll start by reproducing our last result
lung.dens2 <- bivariate.density(lung, h0 = 1)
plot(lung.dens2)

# Now we might want to try adaptive smoothing
lung.dens2 <- bivariate.density(lung, h0 = 1, adapt=TRUE)
plot(lung.dens2)
# This focuses the smoothing around the built up areas, giving a little more
# detail

# We can also do this for the larynx data
larynx.dens2 <- bivariate.density(larynx, h0 = 1, adapt=TRUE)
plot(larynx.dens2)

# A quick look shows that things are a little different with the larynx data
# but could this just be random? After all there aren't many cases of larynx
# cancer.

# One thing we could try is a relative risk surface. This will be the ratio
# of the lung and larynx densities:
rr <- risk(larynx.dens2, lung.dens2)
plot(rr)

# By default the colours here suck. Let's try and get better colours by centering
# things on zero (a relative risk of 1).
# Most of the data is in the range 2 to -2.
col_breaks = c(-10, seq(-2, 2, by=0.1),10)
col_breaks
num_cols = length(col_breaks) - 1
num_cols
cols = hcl.colors(num_cols, "blue-red")
plot(rr, breaks=col_breaks, col = cols)

# Let's add some tolerance curves. These are like p-value usrfaces. It shows
# there seems to be a couple of areas with increased risk of larynx cancer
# comapred to lung cancer

tol <- tolerance(rr)
tol.contour(tol, levels=c(0.1, 0.05, 0.01), add=TRUE)

# OK, now one thing we haven't revealed is that there is an incinerator.
# let's add that to our plot
points(chorley.extra$incin, pch=17, cex=1.5)

# Hmm!

# Ok, say we wanted to have some statistical measure of whether risk is
# increased due to the incinerator. What might we do?

# We could fit a model for the point process `larynx` with covariates
# for where we expect cancer cases to be (using `lung`) and a covariate
# for the distance from the incinerator.

# When we look at ppm, the covariates are of the form of a pixel image,
# which we can get from the `z` part of the density object:
lung.im <- lung.dens2$z

# we start by fitting a model using only the lung.im stuff.
# the first thing to note is that the model runs on the log scale. i.e.
# we're estimating the  log(lambda) of the Poisson point process. It
# thus makes sense to use log(lung.im) for this, as ofcourse lung.im
# is not on the log scale!
out <- ppm(larynx ~ log(lung.im))
out

# we can see that the log(lung.im) covariate is highly significant (Ztest has
# lots of stars.

# what about the incinerator? This involves quite a lot more work!

# To do this, we first create a ppp with all the pixels from the lung
# cancer image (i.e. all the pixels we're running out model over)

pts <- expand.grid(x=lung.im$xcol, y=lung.im$yrow)
dist <- 1/crossdist(pts$x, pts$y, chorley.extra$incin$x,
                    chorley.extra$incin$y)^2
pts <- pts %>% mutate(z = dist)

# chop out the bits we dont want by converting to a ppp then
# back to a data.frame. Theres probably a better way to do this
pts <- as.data.frame(as.ppp(pts, W = larynx$window))

# and then convert to a pixel image
pts.im <- as.im(pts)
plot(log(pts.im))

out <- ppm(larynx ~ log(lung.im) + log(pts.im))
out
# Interestingly, it's not super important, but the effect is in the
# expected direction (+ve coefficient). This might in part be due
# to the other area of increased risk in the west, and the relatively
# low sample size - there's only 4 points right near the incinerator

# Ok, that's interesting enough I guess. Some evidence of an effect, but not lots
# with this very simple model and few data points
