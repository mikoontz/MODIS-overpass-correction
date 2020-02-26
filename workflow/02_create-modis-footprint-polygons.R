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
library(geosphere)

plan(multiprocess)


# global variables --------------------------------------------------------

global_extent <- 
  tibble(lon = c(-180, 180, 180, -180, -180), lat = c(90, 90, -90, -90, 90)) %>% 
  as.matrix() %>% 
  list() %>% 
  st_polygon() %>% 
  st_sfc(crs = 4326) %>% 
  st_sf()

# function to write polygons (a slow step) to disk ------------------------------------

get_MODIS_granule_footprints <- function(satellite, path, year, month, day, yday, GranuleID, StartDateTime, ArchiveSet, OrbitNumber, DayNightFlag, EastBoundingCoord, NorthBoundingCoord, SouthBoundingCoord, WestBoundingCoord, GRingLongitude1, GRingLongitude2, GRingLongitude3, GRingLongitude4, GRingLatitude1, GRingLatitude2, GRingLatitude3, GRingLatitude4) {
  
  pt1 <- c(GRingLongitude1, GRingLatitude1)
  pt2 <- c(GRingLongitude2, GRingLatitude2)
  pt3 <- c(GRingLongitude3, GRingLatitude3)
  pt4 <- c(GRingLongitude4, GRingLatitude4)
  
  linestring_resolution <- units::set_units(10, km)
  # If the footprint overlaps one of the poles, then the WestBoundingCoord
  # and the EastBoundingCoord will be -180 and 180, respectively
  # This requires a special accomodation to ensure the right polygon is made
  if(WestBoundingCoord == -180 | EastBoundingCoord == 180) {
    
    mat <- st_sfc(st_multipoint(matrix(rbind(pt1, pt2, pt3, pt4, pt1), ncol = 2)), crs = 4326)
    
    # Create a single linestring that spans the -180 to 180 global extent
    # which will be used to split a polygon of the global extent
    blade <- 
      mat %>% 
      st_cast("LINESTRING") %>%
      st_segmentize(dfMaxLength = linestring_resolution) %>% 
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
    
    # We know when the footprint crosses the dateline because the
    # WestBoundingCoord is actually greater than the EastBoundingCoord
    # First, we'll deal with the simple case where the footprint doesn't
    # cross the dateline
    
    if (WestBoundingCoord < EastBoundingCoord) {
      
      mat <- st_sfc(st_multipoint(matrix(rbind(pt1, pt2, pt3, pt4, pt1), ncol = 2)), crs = 4326)
      
      footprint <- 
        mat %>% 
        st_cast("LINESTRING") %>%
        st_segmentize(dfMaxLength = linestring_resolution) %>% 
        st_cast("POLYGON") %>% 
        st_sfc(crs = 4326) %>% 
        st_sf()
      
    } else 
      # The condition when the footprint crosses the dateline
      # Rotate the frame of reference so the dateline is now at 0 degrees
      # Build the polygon (that now won't cross the 'dateline' of the 
      # modified projection)
      # Split the geometry at the 0 meridian in the rotated projection
      # (which is the dateline in the original projection)
      # Translate the split polygons to their appropriate places
      if (WestBoundingCoord >= EastBoundingCoord) {
        
        offset <- 180
        
        pt1 <- pt1 + c(offset, 0)
        pt2 <- pt2 + c(offset, 0)
        pt3 <- pt3 + c(offset, 0)
        pt4 <- pt4 + c(offset, 0)
        
        
        pt1 <- pt1 + c(ifelse(pt1[1] > 180, yes = -360, no = 0), 0)
        pt2 <- pt2 + c(ifelse(pt2[1] > 180, yes = -360, no = 0), 0)
        pt3 <- pt3 + c(ifelse(pt3[1] > 180, yes = -360, no = 0), 0)
        pt4 <- pt4 + c(ifelse(pt4[1] > 180, yes = -360, no = 0), 0)
        
        mat <- sf::st_sfc(st_multipoint(matrix(rbind(pt1, pt2, pt3, pt4, pt1), ncol = 2)), crs = 4326)
        
        footprint <- 
          mat %>% 
          sf::st_cast("LINESTRING") %>%
          sf::st_segmentize(dfMaxLength = linestring_resolution) %>% 
          sf::st_cast("POLYGON") %>% 
          sf::st_sfc(crs = 4326) %>% 
          sf::st_sf()
        
        blade <- 
          sf::st_linestring(x = matrix(rbind(c(0, 90), c(0, -90)), ncol = 2)) %>% 
          sf::st_sfc(crs = 4326)
        
        split_poly <- lwgeom::st_split(x = footprint, y = blade)
        
        suppressWarnings({
          polys <- 
            sf::st_collection_extract(x = split_poly, type = "POLYGON") %>% 
            dplyr::mutate(centroid_lon = st_coordinates(st_centroid(.))[, 1])
        })
        
        east_poly <- 
          polys %>% 
          dplyr::filter(centroid_lon < 0) %>% 
          dplyr::summarize() %>% 
          sf::st_geometry()
        
        east_poly <- east_poly + sf::st_sfc(st_point(c(offset, 0)))
        east_poly <- 
          east_poly %>% 
          sf::st_set_crs(4326)
        
        west_poly <-
          polys %>% 
          dplyr::filter(centroid_lon >= 0) %>% 
          dplyr::summarize() %>% 
          sf::st_geometry()
        
        west_poly <- west_poly + sf::st_sfc(st_point(c(-offset, 0)))
        west_poly <- 
          west_poly %>% 
          sf::st_set_crs(4326)
        
        footprint <- 
          lapply(list(east_poly, west_poly), sf::st_as_sf) %>% 
          do.call("rbind", .) %>% 
          dplyr::rename(geometry = x) %>% 
          dplyr::summarize()
        
      } 
  }
  
  # rebuild the sf object row by row
  # cast everything to multipolygon so we can combine the resulting 
  # list elements with data.table::rbindlist()
  obj <- 
    tibble(satellite, path, year, month, day, yday, GranuleID, StartDateTime, ArchiveSet, OrbitNumber, DayNightFlag, EastBoundingCoord, NorthBoundingCoord, SouthBoundingCoord, WestBoundingCoord, GRingLongitude1, GRingLongitude2, GRingLongitude3, GRingLongitude4, GRingLatitude1, GRingLatitude2, GRingLatitude3, GRingLatitude4, geometry = footprint$geometry) %>% 
    st_sf() %>% 
    st_cast("MULTIPOLYGON") 
  
  return(obj)
  
}

get_MODIS_footprints <- function(geoMeta, overwrite = FALSE) {
  
  # First do a check to see whether the file exists locally already. If it does,
  # (and the user doesn't want to overwrite it), then upload that local file
  # to the S3 bucket
  # 
  this_year <- unique(geoMeta$year)
  this_month <- unique(geoMeta$month)
  this_satellite <- unique(geoMeta$satellite)
  
  this_file <- paste0(this_year, "-", this_month, "_", this_satellite, "_5-minute-footprints.gpkg")
  
  this_path <- file.path("data", "data_output", "MODIS-footprints", this_file)
  
  if(!file.exists(this_path) | overwrite) {
    geoMeta <- tidyr::unnest(geoMeta, cols = geoMeta)
    
    footprints <-
      pmap(geoMeta, .f = get_MODIS_granule_footprints)
    
    footprints <- do.call(what = "rbind", args = footprints)
    
    sf::st_write(obj = footprints, dsn = this_path, delete_dsn = TRUE)
  }
  
  system2(command = "aws", args = paste0('s3 cp ', this_path, ' s3://earthlab-mkoontz/MODIS-footprints/', this_file))
  
  # Don't return anything
  return(invisible())
  
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


# get terra geometa61 -----------------------------------------------------
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

# test small snippets of geoMeta for polygon creation ---------------------

# geoMeta_target <- 
#   geoMeta_df %>% 
#   dplyr::filter(path == "data/data_output/MODIS-geoMeta61/AQUA/2005/MYD03_2005-06-02.txt")
# 
# footprints <- get_MODIS_footprints(geoMeta_target)
# plot(footprints[sample(1:nrow(footprints), size = 20), ] %>% st_geometry(), axes = TRUE)

# Find out what footprints have already been processed and don't redo that work
already_processed <- 
  system2(command = "aws", args = "s3 ls s3://earthlab-mkoontz/MODIS-footprints/", stdout = TRUE) %>% 
  tibble::enframe(name = NULL) %>% 
  setNames("files_on_aws") %>% 
  dplyr::mutate(satellite = ifelse(grepl(pattern = "Terra", x = files_on_aws), yes = "Terra", no = "Aqua")) %>% 
  dplyr::mutate(footprints_processed = substr(files_on_aws, start = 32, stop = nchar(files_on_aws))) %>% 
  dplyr::mutate(year_month_sat = ifelse(satellite == "Terra", 
                                        yes = substr(footprints_processed, start = 1, stop = 13),
                                        no = substr(footprints_processed, start = 1, stop = 12)))

geoMeta_list <- 
  geoMeta_df %>% 
  dplyr::mutate(year_month_sat = paste0(year, "-", month, "_", satellite)) %>% 
  dplyr::filter(!(year_month_sat %in% already_processed$year_month_sat)) %>% 
  dplyr::select(-year_month_sat) %>% 
  dplyr::group_by(year, month, satellite) %>% 
  dplyr::group_split()

dir.create("data/data_output/MODIS-footprints", recursive = TRUE)

(start <- Sys.time())
footprints_list <- furrr::future_map(geoMeta_list, get_MODIS_footprints)
(difftime(Sys.time(), start, units = "min"))