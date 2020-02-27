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
library(data.table)
# We want to use run length encoding to dramatically reduce the memory usage while iterating
# through the steps here
if (!require("BiocManager")) {
  install.packages("BiocManager")
  BiocManager::install("S4Vectors")
}
library(S4Vectors)

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

count_overpasses <- function(footprints, raster_template) {
  
  # Ensure geometry sfc column is called "geometry"
  footprints <- st_sf(st_drop_geometry(footprints), geometry = st_geometry(footprints))
  
  overpasses <- 
    pmap(footprints, 
         .f = function(satellite, path, year, month, day, yday, GranuleID, StartDateTime, ArchiveSet, OrbitNumber, DayNightFlag, EastBoundingCoord, NorthBoundingCoord, SouthBoundingCoord, WestBoundingCoord, GRingLongitude1, GRingLongitude2, GRingLongitude3, GRingLongitude4, GRingLatitude1, GRingLatitude2, GRingLatitude3, GRingLatitude4, geometry) {
           
           daynight_r <- 
             raster_template %>% 
             coordinates() %>% 
             as.data.table() %>% 
             setNames(c("longitude", "latitude"))
           
           daynight_r[, `:=`(solar_offset = longitude / 15,
                             StartDateTime = StartDateTime,
                             acq_datetime = StartDateTime + as.duration(2.5 * 60))]
           
           daynight_r[, `:=`(acq_hour = lubridate::hour(acq_datetime),
                             acq_min = lubridate::minute(acq_datetime),
                             local_datetime = acq_datetime + as.duration(solar_offset * 60 * 60))]
           
           daynight_r[, `:=`(local_doy = lubridate::yday(local_datetime),
                             local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24)]
           
           daynight_r[, `:=`(h = (local_hour_decmin - 12) * 15 * pi / 180,
                             phi = latitude * pi / 180,
                             delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))))]
           
           daynight_r[, solar_elev_angle := (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi]
           
           
           daynight_r[, `:=`(day = ifelse(solar_elev_angle >= 0, yes = 1, no = 0),
                             night = ifelse(solar_elev_angle < 0, yes = 1, no = 0))]
           
           r <- 
             fasterize(sf = st_sf(st_sfc(geometry, crs = 4326)), 
                       raster = raster_template, 
                       fun = "count", 
                       background = 0)
           
           daynight_r[, overpass := raster::values(r)]
           daynight_r[, `:=`(day_overpass = day * overpass,
                             night_overpass = night * overpass)]
           
           if(sum(daynight_r$day_overpass) > 0) {
             day_overpass <- S4Vectors::Rle(daynight_r$day_overpass)
           } else {day_overpass <- NULL}
           
           if(sum(daynight_r$night_overpass) > 0) {
             night_overpass <- S4Vectors::Rle(daynight_r$night_overpass)
           } else {night_overpass <- NULL}
           
           
           return(list(day_overpass = day_overpass, night_overpass = night_overpass))
         })
  
  day_sum <- 
    lapply(overpasses, FUN = function(x) x$day_overpass)
  
  day_sum <-
    day_sum[!sapply(day_sum, is.null)] %>% 
    do.call("+", .) %>% 
    S4Vectors::decode()
  
  night_sum <- 
    lapply(overpasses, FUN = function(x) x$night_overpass)
  
  night_sum <-
    night_sum[!sapply(night_sum, is.null)] %>% 
    do.call("+", .) %>% 
    S4Vectors::decode()
  
  day_r <- night_r <- raster_template
  values(day_r) <- day_sum
  values(night_r) <- night_sum
  
  this_year <- unique(footprints$year)
  this_month <- unique(footprints$month)
  this_satellite <- unique(footprints$satellite)
  
  night_file <- paste0(this_year, "-", this_month, "_", this_satellite, "_night-overpass-count_", res(raster_template)[1], ".tif")
  
  day_file <- paste0(this_year, "-", this_month, "_", this_satellite, "_day-overpass-count_", res(raster_template)[1], ".tif")
  
  night_path <- file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1]), night_file)
  
  day_path <- file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1]), day_file)
  
  raster::writeRaster(x = night_r, filename = night_path)
  raster::writeRaster(x = day_r, filename = day_path)
  
  rm(night_r)
  rm(day_r)
  
  system2(command = "aws", 
          args = paste0('s3 cp ', night_path, ' s3://earthlab-mkoontz/MODIS-overpass-counts_', res(raster_template)[1], '/', night_file))
  
  system2(command = "aws", 
          args = paste0('s3 cp ', day_path, ' s3://earthlab-mkoontz/MODIS-overpass-counts_', res(raster_template)[1], '/', day_file))
  
  return(invisible())
}


# Figure out what has been processed already

raster_template <- r_0.25
  
overpasses_processed <-
  system2(command = "aws", 
          args = paste0("s3 ls s3://earthlab-mkoontz/MODIS-overpass-counts_", res(raster_template)[1], "/"), 
          stdout = TRUE) %>% 
  tibble::enframe(name = NULL) %>% 
  setNames("files_on_aws") %>% 
  dplyr::mutate(satellite = ifelse(grepl(pattern = "Terra", x = files_on_aws), 
                                   yes = "Terra", 
                                   no = "Aqua")) %>% 
  dplyr::mutate(night_or_day = ifelse(grepl(pattern = "night", x = files_on_aws), 
                                      yes = "night", 
                                      no = "day")) %>% 
  dplyr::mutate(overpasses_processed = substr(files_on_aws, 
                                              start = 32, 
                                              stop = nchar(files_on_aws) - 4)) %>% 
  dplyr::mutate(year_month_sat = ifelse(satellite == "Terra", 
                                        yes = substr(overpasses_processed, start = 1, stop = 13),
                                        no = substr(overpasses_processed, start = 1, stop = 12)))

dir.create(paste0("data/data_output/MODIS-overpass-counts_", res(raster_template)[1], "/"), recursive = TRUE)

overpasses_to_process <- 
  system2(command = "aws", 
          args = "s3 ls s3://earthlab-mkoontz/MODIS-footprints/", 
          stdout = TRUE) %>% 
  tibble::enframe(name = NULL) %>% 
  setNames("files_on_aws") %>% 
  dplyr::mutate(satellite = ifelse(grepl(pattern = "Terra", x = files_on_aws), 
                                   yes = "Terra", 
                                   no = "Aqua")) %>% 
  dplyr::mutate(footprints_processed = substr(files_on_aws, 
                                              start = 32, 
                                              stop = nchar(files_on_aws))) %>% 
  dplyr::mutate(year_month_sat = ifelse(satellite == "Terra", 
                                        yes = substr(footprints_processed, start = 1, stop = 13),
                                        no = substr(footprints_processed, start = 1, stop = 12))) %>% 
  dplyr::filter(!(year_month_sat %in% unique(overpasses_processed$year_month_sat))) %>% 
  dplyr::pull(footprints_processed)


plan(multiprocess, workers = 2)

furrr::future_map(overpasses_to_process, .f = function(this_footprint) {
  
  local_path <- paste0("data/data_output/MODIS-footprints/", this_footprint)
  
  if(!file.exists(local_path)) {
    system2(command = "aws", args = paste0('s3 cp s3://earthlab-mkoontz/MODIS-footprints/', this_footprint, ' ', local_path))
  }
  
  footprints <- sf::st_read(local_path)
  
  this_overpass <- count_overpasses(footprints, raster_template)
  
  # Don't delete local file just yet.
  # unlink(local_path)
  
})