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

# NOTE: Must use latest versions of GDAL, GEOS, and PROJ for this to work!
# This SO answer was helpful for making sure {sf} linked to the latest
# versions: https://stackoverflow.com/questions/44973639/trouble-installing-sf-due-to-gdal

aqua_years <- 2002:2019
terra_years <- 2000:2019
geoMeta_source <- "s3" # could also be "laads-daac" or "local"

if(geoMeta_source != "local") {
  if(geoMeta_source == "laads-daac") {
    api_keys <- read_csv("data/data_raw/LAADS-DAAC_api-keys.csv")
    
    if(!dir.exists(paths = file.path("data", "data_output", "MODIS-geoMeta61", "TERRA"))) {
      
      dir.create(path = file.path("data", "data_output", "MODIS-geoMeta61", "TERRA"))
      
      lapply(terra_years, FUN = function(this_year) {
        system2(command = "wget", args = paste0('-e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov/archive/geoMeta/61/TERRA/', this_year, '/" --header "Authorization: Bearer ', dplyr::filter(api_keys, satellite == "Terra") %>% dplyr::pull(key), '" -P data/data_output/MODIS-geoMeta61/'))
      })
      
    }
    if(!dir.exists(paths = file.path("data", "data_output", "MODIS-geoMeta61", "AQUA"))) {
      
      dir.create(path = file.path("data", "data_output", "MODIS-geoMeta61", "AQUA"))
      
      lapply(aqua_years, FUN = function(this_year) {
        system2(command = "wget", args = paste0('-e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov/archive/geoMeta/61/AQUA/', this_year, '/" --header "Authorization: Bearer ',  dplyr::filter(api_keys, satellite == "Aqua") %>% dplyr::pull(key), '" -P data/data_output/MODIS-geoMeta61/'))
      })
      
    }
  } else if(geoMeta_source == "s3") {
    system2(command = "aws", args = paste0('s3 sync s3://earthlab-mkoontz/MODIS-geoMeta61/ data/data_output/MODIS-geoMeta61/'))
  }
}

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

# function to write polygons (a slow step) to disk ------------------------------------

get_MODIS_footprints <- function(geoMeta) {
  
  geoMeta <- tidyr::unnest(geoMeta, cols = geoMeta)
  
  footprints <-
    pmap(geoMeta, .f = function(satellite, path, year, month, day, yday, GranuleID, StartDateTime, ArchiveSet, OrbitNumber, DayNightFlag, EastBoundingCoord, NorthBoundingCoord, SouthBoundingCoord, WestBoundingCoord, GRingLongitude1, GRingLongitude2, GRingLongitude3, GRingLongitude4, GRingLatitude1, GRingLatitude2, GRingLatitude3, GRingLatitude4) {
      
      pt1 <- c(GRingLongitude1, GRingLatitude1)
      pt2 <- c(GRingLongitude2, GRingLatitude2)
      pt3 <- c(GRingLongitude3, GRingLatitude3)
      pt4 <- c(GRingLongitude4, GRingLatitude4)
      
      mat <- st_sfc(st_multipoint(matrix(rbind(pt1, pt2, pt3, pt4, pt1), ncol = 2)), crs = 4326)
      
      footprint_line <-
        mat %>% 
        st_cast("LINESTRING") %>%
        st_segmentize(dfMaxLength = units::set_units(10, km))
      
      if(SouthBoundingCoord <= -89.9 | NorthBoundingCoord >= 89.9) {
        
        blade <- 
          footprint_line %>% 
          st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) %>% 
          st_cast("POINT") %>% 
          sf::st_coordinates()
        
        blade <- blade[order(blade[, "X"]), ]
        
        blade <-
          blade %>% 
          st_linestring() %>% 
          st_sfc(crs = 4326)
        
        split_globe <- lwgeom::st_split(x = global_extent, y = blade)
        
        suppressWarnings({
          polys <- 
            st_collection_extract(x = split_globe, type = "POLYGON") %>% 
            dplyr::mutate(centroid_lat = st_coordinates(st_centroid(.))[, 2])
        })
        
        if(SouthBoundingCoord > 0) {
          footprint <- 
            polys %>% 
            filter(centroid_lat == max(centroid_lat)) %>% 
            dplyr::select(-centroid_lat)
          
        } else {
          footprint <- 
            polys %>% 
            filter(centroid_lat == min(centroid_lat)) %>% 
            dplyr::select(-centroid_lat)
          
        }
      } else {
        footprint <- 
          footprint_line %>% 
          st_cast("POLYGON") %>% 
          st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) %>% 
          st_sfc(crs = 4326) %>% 
          st_sf(geometry = .)
        
      }
      
      # rebuild the sf object row by row
      # cast everything to multipolygon so we can combine the resulting 
      # list elements with data.table::rbindlist()
      obj <- 
        tibble(satellite, path, year, month, day, yday, GranuleID, StartDateTime, ArchiveSet, OrbitNumber, DayNightFlag, EastBoundingCoord, NorthBoundingCoord, SouthBoundingCoord, WestBoundingCoord, GRingLongitude1, GRingLongitude2, GRingLongitude3, GRingLongitude4, GRingLatitude1, GRingLatitude2, GRingLatitude3, GRingLatitude4, geometry = footprint$geometry) %>% 
        st_sf() %>% 
        st_cast("MULTIPOLYGON") 
      
      return(obj)
      
    })
  
  footprints <- do.call(what = "rbind", args = footprints)
  
  this_year <- unique(footprints$year)
  this_month <- unique(footprints$month)
  this_file <- paste0(this_year, "-", this_month, "_5-minute-footprints.gpkg")
  
  this_path <- file.path("data", "data_output", "MODIS-footprints", this_file)
  
  sf::st_write(obj = footprints, dsn = this_path, delete_dsn = TRUE)
  
  system2(command = "aws", args = paste0('s3 cp ', this_path, ' s3://earthlab-mkoontz/MODIS-footprints/', this_file))
  
  return(footprints)
  
}

# read all the daily geoMeta data ----------------------------------------

aqua <- 
  list.files("data/data_output/MODIS-geoMeta61/AQUA", recursive = TRUE, full.names = TRUE) %>% 
  tibble(satellite = "Aqua",
         path = .,
         year = substr(path, start = 39, stop = 42),
         month = substr(path, start = 55, stop = 56),
         day = substr(path, start = 58, stop = 59),
         yday = lubridate::yday(ymd(substr(path, start = 50, stop = 59)))) %>% 
  dplyr::mutate(geoMeta = purrr::map(path, .f = function(this_path) {
    read_delim(this_path, 
               delim = ",", 
               skip = 2,
               col_types = cols(
                 `# GranuleID` = col_character(),
                 StartDateTime = col_datetime(format = ""),
                 ArchiveSet = col_double(),
                 OrbitNumber = col_double(),
                 DayNightFlag = col_character(),
                 EastBoundingCoord = col_double(),
                 NorthBoundingCoord = col_double(),
                 SouthBoundingCoord = col_double(),
                 WestBoundingCoord = col_double(),
                 GRingLongitude1 = col_double(),
                 GRingLongitude2 = col_double(),
                 GRingLongitude3 = col_double(),
                 GRingLongitude4 = col_double(),
                 GRingLatitude1 = col_double(),
                 GRingLatitude2 = col_double(),
                 GRingLatitude3 = col_double(),
                 GRingLatitude4 = col_double()
               )) %>%
      dplyr::rename(GranuleID = `# GranuleID`)
    
  }))

terra <- 
  list.files("data/data_output/MODIS-geoMeta61/TERRA", recursive = TRUE, full.names = TRUE) %>% 
  tibble(satellite = "Terra",
         path = .,
         year = substr(path, start = 40, stop = 43),
         month = substr(path, start = 56, stop = 57),
         day = substr(path, start = 59, stop = 60),
         yday = lubridate::yday(ymd(substr(path, start = 51, stop = 60)))) %>% 
  dplyr::mutate(geoMeta = purrr::map(path, .f = function(this_path) {
    read_delim(this_path, 
               delim = ",", 
               skip = 2,
               col_types = cols(
                 `# GranuleID` = col_character(),
                 StartDateTime = col_datetime(format = ""),
                 ArchiveSet = col_double(),
                 OrbitNumber = col_double(),
                 DayNightFlag = col_character(),
                 EastBoundingCoord = col_double(),
                 NorthBoundingCoord = col_double(),
                 SouthBoundingCoord = col_double(),
                 WestBoundingCoord = col_double(),
                 GRingLongitude1 = col_double(),
                 GRingLongitude2 = col_double(),
                 GRingLongitude3 = col_double(),
                 GRingLongitude4 = col_double(),
                 GRingLatitude1 = col_double(),
                 GRingLatitude2 = col_double(),
                 GRingLatitude3 = col_double(),
                 GRingLatitude4 = col_double()
               )) %>%
      dplyr::rename(GranuleID = `# GranuleID`)
    
  }))

geoMeta_df <- rbind(aqua, terra)

# Find out what footprints have already been processed and don't redo that work
already_processed <- 
  system2(command = "aws", args = "s3 ls s3://earthlab-mkoontz/MODIS-footprints/", stdout = TRUE) %>% 
  tibble::enframe(name = NULL) %>% 
  setNames("files_on_aws") %>% 
  dplyr::mutate(footprints_processed = substr(files_on_aws, start = 32, stop = 63)) %>% 
  dplyr::mutate(year_month = substr(footprints_processed, start = 1, stop = 7))

geoMeta_list <- 
  geoMeta_df %>% 
  dplyr::mutate(year_month = paste(year, month, sep = "-")) %>% 
  dplyr::filter(!(year_month %in% already_processed$year_month)) %>% 
  dplyr::select(-year_month) %>% 
  dplyr::group_by(year, month) %>% 
  dplyr::group_split()

dir.create("data/data_output/MODIS-footprints", recursive = TRUE)

(start <- Sys.time())
footprints_list <- furrr::future_map(geoMeta_list, get_MODIS_footprints)
(difftime(Sys.time(), start, units = "min"))

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


# data.table implementation of solar elevation angle ----------------------

this_afd[, acq_hour := floor(ACQ_TIME / 100)]
this_afd[, acq_min := ((ACQ_TIME / 100) - acq_hour) * 100]
this_afd[, acq_datetime := as.POSIXct((ACQ_DATE / 1000) + (acq_hour * 3600) + (acq_min * 60), 
                                      origin = "1970-01-01", 
                                      tz = "America/Los_Angeles")]

this_afd[, `:=`(acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = LONGITUDE / 15,
                hemisphere = ifelse(LATITUDE >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"))]

this_afd[, acq_datetime_local := acq_datetime + as.duration(solar_offset * 60 * 60)]

this_afd[, `:=`(local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24)]

this_afd[, `:=`(local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2)]

# https://en.wikipedia.org/wiki/Solar_zenith_angle
# https://en.wikipedia.org/wiki/Position_of_the_Sun#Declination_of_the_Sun_as_seen_from_Earth
# https://en.wikipedia.org/wiki/Hour_angle

this_afd[, `:=`(h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = LATITUDE * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))))]

this_afd[, solar_elev_ang := (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi]





























# # write to disk
# write.csv(overpass_corrections_0.25, file = "analyses/analyses_output/aqua-terra-overpass-corrections-table_0.25-degree-grid.csv", row.names = FALSE)
# write.csv(overpass_corrections_2.5, file = "analyses/analyses_output/aqua-terra-overpass-corrections-table_2.5-degree-grid.csv", row.names = FALSE)
# 
# # save the visualization to disk
# png("figures/aqua-terra-overpass-corrections-map_0.25-degree-grid.png")
# plot(orbit_overlap_0.25, col = viridis(30))
# plot(st_as_sf(ne_coastline()) %>% st_geometry(), add = TRUE)
# dev.off()
# 
# png("figures/aqua-terra-overpass-corrections-map_2.5-degree-grid.png")
# plot(orbit_overlap_2.5, col = viridis(30))
# plot(st_as_sf(ne_coastline()) %>% st_geometry(), add = TRUE)
# dev.off()
# 
# # save the empirical model plot to disk
# ggplot(overpass_corrections_0.25 %>% filter(lat %in% c(seq(-83.5, -70, by = 0.25), -69:69, seq(70, 83.5, by = 0.25))), aes(x = lat, y = mean_overpasses)) +
#   geom_point(cex = 0.3) +
#   theme_bw() +
#   geom_ribbon(aes(ymin = min_overpasses, ymax = max_overpasses), fill = "red", alpha = 0.1)
# ggsave("figures/aqua-terra-overpass-corrections-function_0.25-degree-grid.png")
# 
# ggplot(overpass_corrections_2.5 %>% filter(lat %in% c(seq(-83.5, -70, by = 0.25), -69:69, seq(70, 83.5, by = 0.25))), aes(x = lat, y = mean_overpasses)) +
#   geom_point(cex = 0.3) +
#   theme_bw() +
#   geom_ribbon(aes(ymin = min_overpasses, ymax = max_overpasses), fill = "red", alpha = 0.1)
# ggsave("figures/aqua-terra-overpass-corrections-function_2.5-degree-grid.png")

