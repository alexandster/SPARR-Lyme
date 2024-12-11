#################################
## Copy of Alex's read-in code ##
#################################
library(sparr)
library(dplyr)
library(sf)
library(gt)
library(gtExtras)
library(raster)
library(tmap)
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#read Baltimore County boundary
geom <- st_read("CONUS_counties.shp") %>%
  subset(., GEOID == "24005") %>%
  st_transform(., crs = 2248) %>%
  st_coordinates(.) %>%
  as.data.frame(.) %>%
  subset(., L2 == 4)
#read table
df <- read.csv("lyme_alex.csv") %>%
  st_as_sf(., coords = c("long", "Lat")) %>%
  st_set_crs(4326) %>%
  st_transform(., crs = 2248)
#frequency
table(df$Type)
#separate 
df_0 <- subset(df, Type == 0) %>% #controls
  st_coordinates(.)
df_1 <- subset(df, Type == 1) %>% #cases
  st_coordinates(.)
#ppp object window
w <- owin(poly=list(x=rev(geom$X),y=rev(geom$Y))) #window
glyme <- ppp(x = df$x, y = df$y, window = w, marks = df$class)
#ppp
df_0_ppp <- ppp(df_0[, "X"], df_0[, "Y"], window = w)
df_1_ppp <- ppp(df_1[, "X"], df_1[, "Y"], window = w)
#################################
#################################
#################################



# Tilman's scratchpad

lyppp <- superimpose(df_0_ppp,df_1_ppp)
marks(lyppp)<- factor(rep(c("control","case"),c(npoints(df_0_ppp),npoints(df_1_ppp))))
plot(lyppp)
summary(lyppp)

OS(lyppp)
h <- LSCV.risk(lyppp,method="kelsall-diggle");h 
h <- LSCV.risk(lyppp,method="davies");h 
# SLIK.adapt(df_0_ppp) ## This gives h around 3980. Too small.

lyrr1 <- risk(lyppp,h0=7200,tolerate=TRUE) # using roughly the "Davies" bandwidth. In general, I don't recommend using a bandwidth smaller than 7000 for these data
plot(lyrr1)
points(df_0_ppp,pch=3,col="peachpuff4")
points(df_1_ppp,pch=19,col="seagreen3") ## Just to see where the cases fall. It's important for interpretation of tolerance contours. For example, if we get tolerance contours occurring due simply to one point, this is not particularly reliable!


par(mfrow=c(1,2))  #comparing asymmetric and symmetric (pooled) surfaced
lyrr2 <- risk(lyppp,adapt=TRUE,h0=7200,tolerate=TRUE,pilot.symmetry="none") 
plot(lyrr2,main="adaptive asymmetric, h0=7200")
points(df_0_ppp,pch=3,col="peachpuff4")
points(df_1_ppp,pch=19,col="seagreen3")

lyrr3 <- risk(lyppp,adapt=TRUE,h0=7200,tolerate=TRUE,pilot.symmetry="pooled") 
plot(lyrr3,main="adaptive symmetric (pooled), h0=7200")
points(df_0_ppp,pch=3,col="peachpuff4")
points(df_1_ppp,pch=19,col="seagreen3")


# Testing out small bandwidth on fixed bandwidth rr surface -- this is too small -- prefer adaptive estimation for this pattern
ffix <- bivariate.density(df_1_ppp,h0=3980)
gfix <- bivariate.density(df_0_ppp,h0=3980)
plot(ffix)
plot(gfix)
rfix <- risk(ffix,gfix,tolerate=TRUE)
plot(rfix)
points(df_1_ppp,pch=19,col="seagreen3")

# Comparing 'rfix' with an adaptive estimate using the same small bandwidth. Better, but still 'noisy' -- notice that a couple of significant contours are due solely to single data points.
fada <- bivariate.density(df_1_ppp,h0=3980,adapt=TRUE)
gada <- bivariate.density(df_0_ppp,h0=3980,adapt=TRUE)
plot(fada)
plot(gada)
rada <- risk(fada,gada,tolerate=TRUE,adapt=T)
plot(rada)
points(df_1_ppp,pch=19,col="seagreen3")


### Conclusion: Use an adaptive estimate here with a bandwidth of at least 7000. Depending on how sensitive you want the risk surface to be to the specific/individual case locations, you could choose either an asymmetric (more sensitive) or symmetric-pooled (slightly more conservative) estimator.
