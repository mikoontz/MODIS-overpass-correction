# Purpose: overpass correction for Aqua / Terra to account for more possible
# overpasses at higher latitudes.

library(tidyverse)
library(lubridate)
library(sf)
library(lwgeom)
library(raster)
library(fasterize)
library(viridis)
library(purrr)
library(furrr)

plan(multiprocess)

# global variables --------------------------------------------------------

global_extent <- 
  tibble(lon = c(-180, 180, 180, -180, -180), lat = c(90, 90, -90, -90, 90)) %>% 
  as.matrix() %>% 
  list() %>% 
  st_polygon() %>% 
  st_sfc(crs = 4326) %>% 
  st_sf()

# raster templates

r_0.25 <- raster::raster("data/data_raw/grid_0_25_degree_vars_modis_D_AFC_num_April_2001.tif")

r_2.5 <- raster::raster("data/data_raw/grid_2_5_degree_vars_modis_D_AFC_num_April_2001.tif")

# build rasterized count of overpasses ------------------------------------

count_overpasses <- function(footprints) {
  
  overpasses <- 
    pmap(footprints, 
         .f = function(satellite, path, year, month, day, yday, GranuleID, StartDateTime, ArchiveSet, OrbitNumber, DayNightFlag, EastBoundingCoord, NorthBoundingCoord, SouthBoundingCoord, WestBoundingCoord, GRingLongitude1, GRingLongitude2, GRingLongitude3, GRingLongitude4, GRingLatitude1, GRingLatitude2, GRingLatitude3, GRingLatitude4, geometry) {
           
           daynight_r <- 
             r_0.25 %>% 
             coordinates() %>% 
             as.data.frame() %>% 
             setNames(c("longitude", "latitude")) %>% 
             dplyr::mutate(solar_offset = longitude / 15,
                           acq_datetime = StartDateTime,
                           acq_hour = hour(acq_datetime),
                           acq_min = minute(acq_datetime),
                           local_datetime = acq_datetime + as.duration(solar_offset * 60 * 60),
                           local_doy = lubridate::yday(local_datetime),
                           local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                           local_solar_hour_decmin_round = round(local_hour_decmin),
                           local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                           h = (local_hour_decmin - 12) * 15 * pi / 180,
                           phi = latitude * pi / 180,
                           delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                           solar_elev_angle = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                           day = ifelse(solar_elev_angle >= 0, yes = 1, no = 0),
                           night = ifelse(solar_elev_angle < 0, yes = 1, no = 0))
           
           day_r <- night_r <- r_0.25
           values(day_r) <- daynight_r$day
           values(night_r) <- daynight_r$night
           
           r <- 
             fasterize(sf = st_sf(st_sfc(geometry, crs = 4326)), 
                       raster = r_0.25, 
                       fun = "count", 
                       background = 0)
           
           day_r <- day_r * r
           night_r <- night_r * r
           
           return(list(day = day_r, night = night_r))
         })
  
  day_sum <- lapply(overpasses, FUN = function(x) x$day) %>% stack() %>% sum()
  night_sum <- lapply(overpasses, FUN = function(x) x$night) %>% stack() %>% sum()
  
  par(mfrow = c(1, 2))
  plot(day_sum, col = viridis::viridis(6), main = "Day overpass count")
  plot(night_sum, col = viridis::viridis(6), main = "Night overpass count")
  
}

footprints_to_process <-
  system2(command = "aws", args = "s3 ls s3://earthlab-mkoontz/MODIS-footprints/", stdout = TRUE) %>% 
  tibble::enframe(name = NULL) %>% 
  setNames("files_on_aws") %>% 
  dplyr::mutate(footprints_processed = substr(files_on_aws, start = 32, stop = 63)) %>% 
  dplyr::mutate(year_month = substr(footprints_processed, start = 1, stop = 7))

this_footprint <- 
  footprints_to_process %>% 
  dplyr::filter(year_month == "2005-06") %>% 
  slice(1)

s3_file <- this_footprint$footprints_processed

if (!file.exists(paste0('data/data_output/MODIS-footprints/', s3_file))) {
system2(command = "aws", args = paste0('s3 cp s3://earthlab-mkoontz/MODIS-footprints/', s3_file, ' data/data_output/MODIS-footprints/', s3_file))
}

footprints <- sf::st_read(paste0('data/data_output/MODIS-footprints/', s3_file), 
                          stringsAsFactors = FALSE)

plot(footprints[1:100, ] %>% st_geometry())

unique(footprints$satellite)

nrow(footprints)
length(unique(footprints$StartDateTime))




satellite <- footprints$satellite[3]
path <- footprints$path[3]
year <- footprints$year[3]
month <- footprints$month[3]
day <- footprints$day[3]
yday <- footprints$yday[3]
GranuleID <- footprints$GranuleID[3]
StartDateTime <- footprints$StartDateTime[3]
ArchiveSet <- footprints$ArchiveSet[3]
OrbitNumber <- footprints$OrbitNumber[3]
DayNightFlag <- footprints$DayNightFlag[3]
EastBoundingCoord <- footprints$EastBoundingCoord[3]
NorthBoundingCoord <- footprints$NorthBoundingCoord[3]
SouthBoundingCoord <- footprints$SouthBoundingCoord[3]
WestBoundingCoord <- footprints$WestBoundingCoord[3]
GRingLongitude1 <- footprints$GRingLongitude1[3]
GRingLongitude2 <- footprints$GRingLongitude2[3]
GRingLongitude3 <- footprints$GRingLongitude3[3]
GRingLongitude4 <- footprints$GRingLongitude4[3]
GRingLatitude1 <- footprints$GRingLatitude1[3]
GRingLatitude2 <- footprints$GRingLatitude2[3]
GRingLatitude3 <- footprints$GRingLatitude3[3]
GRingLatitude4 <- footprints$GRingLatitude4[3]
geometry <- footprints$geometry[3]