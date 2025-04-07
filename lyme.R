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

# Tilman's scratchpad
lyppp <- superimpose(df_0_ppp,df_1_ppp)
marks(lyppp)<- factor(rep(c("control","case"),c(npoints(df_0_ppp),npoints(df_1_ppp))))
plot(lyppp)
summary(lyppp)

#bandwidth selection
OS(lyppp)
h <- LSCV.risk(lyppp,method="kelsall-diggle");h 
h <- LSCV.risk(lyppp,method="davies");h 

#compute SPARR
t1 <- Sys.time()
#lyrr <- risk(lyppp, adapt=TRUE, h0=7200, resolution = 450, tolerate=TRUE, pilot.symmetry="none") 
lyrr <- risk(lyppp, adapt=TRUE, h0=7200, resolution = 100, tolerate=TRUE, pilot.symmetry="none") 
t2 <- Sys.time()
t2-t1 

plot(lyrr,main="adaptive asymmetric, h0=7200")
points(df_0_ppp,pch=3,col="peachpuff4")
points(df_1_ppp,pch=19,col="seagreen3")

#classify points
rho.class <- tol.classify(lyrr, cutoff = 0.05)

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
  st_set_crs(2248)
pcpolys$ID <- 1:length(rho.class[["finsplit"]])
#st_write(pcpolys,"pcpolys.shp", append = FALSE)

Area <- st_area(pcpolys) %>%
  units::set_units(., value = km^2) #Takes care of units #Take care of units
Case_density <- Cases/Area                    #

#cluster report table
df_res <- data.frame(ID, N, Cases, Controls, Risk, Case_density, Area)
#write.csv(df_res, "clusters.csv", row.names = FALSE)

#risk surface to raster
r <- raster(lyrr$rr)
crs(r) <- crs(pcpolys)

#map it
tm <- tm_shape(r) +
  tm_raster(col.scale = tm_scale(style = "quantile", 
                                 values = "brewer.oranges")) +
  tm_shape(pcpolys) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(frame = FALSE, legend.show = TRUE)
tm

#tmap_save(tm, "lyme_clusters.jpg", width = 4, height = 4)

