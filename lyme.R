#install.packages("gt")
#install.packages("gtExtras")
library(sparr)
library(dplyr)
library(sf)
library(gt)
library(gtExtras)
library(raster)
library(tmap)
library(ggplot2)
library(spatstat)
library(spatstat.utils)

# set workspace----
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


#kde 0: control
con <- bivariate.density(pp=df_0_ppp, h0=OS(df_0_ppp)/2, adapt=FALSE, resolution=750, verbose=TRUE, parallelise = 7)

#kde 1: case
cas <- bivariate.density(pp=df_1_ppp, h0=OS(df_0_ppp)/2, adapt=FALSE, resolution=750, verbose=TRUE, parallelise = 7)

#risk
rho <- risk(cas, con, tolerate = TRUE)

plot(rho, tol.show = TRUE)

#classify points
rho.class <- tol.classify(rho, cutoff = 0.05)

#plot
plot(rho)
points(rho.class$fin,col=2)
points(rho.class$fout)

#cluster report table
ID <- 1:length(rho.class[["finsplit"]])       #cluster identifier
Cases <- lengths(rho.class[["finsplit"]])     #case count
Controls <- lengths(rho.class[["ginsplit"]])  #control count
N <- Cases + Controls                         #point count
Risk <- Cases/N                               #

#contours to sf
pcpolys <- rho.class$pcpolys %>%
  lapply(., FUN = st_as_sf) %>%
  do.call(rbind, .) %>%
  st_set_crs(26912)
pcpolys$ID <- 1:length(rho.class[["finsplit"]])
st_write(pcpolys,"pcpolys.shp")

Area <- st_area(pcpolys) #Take care of units
Case_density <- Cases/Area                    #

#style it
df_res <- data.frame(ID, N, Cases, Controls, Risk, Case_density, Area) %>%
  gt() %>%
  gt_theme_nytimes() %>%
  tab_header(title = "Clusters of Lyme disease in Baltimore, MD")
df_res
