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


# dateline split ----------------------------------------------------------

.breakAtDateLine <- function(x) {
  r <- range(x[,1]) 
  r <- r[2] - r[1]
  # if (r > 200) {
    dif <- abs(x[-nrow(x),1] - x[-1,1])
    tr <- which(dif==max(dif))
    x1 <- x[1:tr, ,drop=FALSE]
    x2 <- x[(tr+1):nrow(x), ,drop=FALSE]
    if (x1[tr,1] < 0) { 
      x1[tr,1] <- -180 
      x2[1,1] <- 180
    } else { 
      x1[tr,1] <- 180 
      x2[1,1] <- -180 
    }
    if (nrow(x1) <= 1) {
      res <- x2
    } else if (nrow(x2) <= 1) {
      res <- x1
    } else {
      res <- list(x1, x2)
    }
    return(res)
  # } 
  # return(x)
}

# global variables --------------------------------------------------------

global_extent <- 
  tibble(lon = c(-180, 180, 180, -180, -180), lat = c(90, 90, -90, -90, 90)) %>% 
  as.matrix() %>% 
  list() %>% 
  st_polygon() %>% 
  st_sfc(crs = 4326) %>% 
  st_sf()

# function to write polygons (a slow step) to disk ------------------------------------

get_MODIS_footprints <- function(geoMeta) {
  
  geoMeta <- tidyr::unnest(geoMeta, cols = geoMeta)
  
  footprints <-
    pmap(geoMeta, .f = function(satellite, path, year, month, day, yday, GranuleID, StartDateTime, ArchiveSet, OrbitNumber, DayNightFlag, EastBoundingCoord, NorthBoundingCoord, SouthBoundingCoord, WestBoundingCoord, GRingLongitude1, GRingLongitude2, GRingLongitude3, GRingLongitude4, GRingLatitude1, GRingLatitude2, GRingLatitude3, GRingLatitude4) {
      
      pt1 <- c(GRingLongitude1, GRingLatitude1)
      pt2 <- c(GRingLongitude2, GRingLatitude2)
      pt3 <- c(GRingLongitude3, GRingLatitude3)
      pt4 <- c(GRingLongitude4, GRingLatitude4)
      
      # We know when the footprint crosses the dateline because the
      # WestBoundingCoord is actually greater than the EastBoundingCoord
      if (WestBoundingCoord > EastBoundingCoord) {
        
        line1 <- 
          st_linestring(matrix(rbind(pt1, pt2), ncol = 2)) %>% 
          st_sfc(crs = 4326) %>% 
          st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=100")) %>% 
          st_segmentize(dfMaxLength = units::set_units(10, km)) %>% 
          st_wrap_dateline()
        
        line2 <- 
          st_linestring(matrix(rbind(pt2, pt3), ncol = 2)) %>% 
          st_sfc(crs = 4326) %>% 
          st_segmentize(dfMaxLength = units::set_units(10, km))

        line3 <- 
          st_sfc(st_linestring(rbind(pt3, pt4)), crs = 4326) %>% 
          st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=190"))   %>%
          st_segmentize(units::set_units(10, km)) %>% 
          st_wrap_dateline()

        line4 <- 
          st_linestring(matrix(rbind(pt4, pt1), ncol = 2)) %>% 
          st_sfc(crs = 4326)
        
        footprint_line <- do.call("rbind", lapply(list(line1, line2, line3, line4), st_as_sf))
        
      } else {
        
        mat <- st_sfc(st_multipoint(matrix(rbind(pt1, pt2, pt3, pt4, pt1), ncol = 2)), crs = 4326)
        
        # 
        # line1 <- geosphere::gcIntermediate(pt1, pt2, breakAtDateLine = TRUE, sp = TRUE) %>% st_as_sf()
        # 
        # line2 <- geosphere::gcIntermediate(pt2, pt3, breakAtDateLine = TRUE, sp = TRUE) %>% st_as_sf()
        # line3 <- geosphere::gcIntermediate(pt3, pt4)
        # test1 <- .breakAtDateLine(line3)[[1]] %>% as.data.frame() %>% st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>% summarize() %>% st_cast("LINESTRING")
        # 
        # test2 <- .breakAtDateLine(line3)[[2]] %>% as.data.frame() %>% st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>% summarize() %>% st_cast("LINESTRING")
        # 
        # line3 <- do.call("rbind", list(test1, test2))
        # line4 <- geosphere::gcIntermediate(pt4, pt1, breakAtDateLine = TRUE, sp = TRUE) %>% st_as_sf()
        
        plot(line1, axes = TRUE)
        plot(line2, add = TRUE)
        plot(line3, add = TRUE)
        plot(line4, add = TRUE)
        
        footprint_line <-
          mat %>% 
          st_cast("MULTILINESTRING") %>%
          st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=190")) %>%
          st_segmentize(dfMaxLength = units::set_units(10, km)) %>% 
          st_wrap_dateline()
      }

      if(WestBoundingCoord == -180 | EastBoundingCoord == 180) {
        
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

geoMeta <- geoMeta_df %>% filter(path == "data/data_output/MODIS-geoMeta61/AQUA/2005/MYD03_2005-06-02.txt")

geoMeta <- tidyr::unnest(geoMeta, cols = geoMeta)

# test <- geoMeta %>% filter(GranuleID %in% c("MYD03.A2005153.0030.061.2018022195730.hdf", 
#                                             "MYD03.A2005153.0035.061.2018022195452.hdf"))
# 
# geoMeta <- geoMeta %>% filter(GranuleID == "MYD03.A2005153.0030.061.2018022195730.hdf")

test <-
  pmap(geoMeta[8:10, ], .f = function(satellite, path, year, month, day, yday, GranuleID, StartDateTime, ArchiveSet, OrbitNumber, DayNightFlag, EastBoundingCoord, NorthBoundingCoord, SouthBoundingCoord, WestBoundingCoord, GRingLongitude1, GRingLongitude2, GRingLongitude3, GRingLongitude4, GRingLatitude1, GRingLatitude2, GRingLatitude3, GRingLatitude4) {
    
    pt1 <- c(GRingLongitude1, GRingLatitude1)
    pt2 <- c(GRingLongitude2, GRingLatitude2)
    pt3 <- c(GRingLongitude3, GRingLatitude3)
    pt4 <- c(GRingLongitude4, GRingLatitude4)
    
    mat <- st_sfc(st_multipoint(matrix(rbind(pt1, pt2, pt3, pt4, pt1), ncol = 2)), crs = 4326) %>% 
      st_sf() %>% 
      st_cast("POINT") %>% 
      rowid_to_column() %>% 
      dplyr::mutate(rowid = ifelse(rowid == 5, yes = 1, no = rowid))
    
    
  }) %>% do.call("rbind", .)

plot(test, axes = TRUE)

test2 <- 
  geoMeta %>% 
  dplyr::mutate(north_pole_overlap = 
                  NorthBoundingCoord - GRingLatitude1 > 0 |
                  NorthBoundingCoord - GRingLatitude2 > 0 |
                  NorthBoundingCoord - GRingLatitude3 > 0 |
                  NorthBoundingCoord - GRingLatitude4 > 0)

abs(-17.4247071) - abs(c(-34.8902211, -38.665270, -20.383383, -17.3517469))

satellite <- geoMeta$satellite[7]
path <- geoMeta$path[7]
year <- geoMeta$year[7]
month <- geoMeta$month[7]
day <- geoMeta$day[7]
yday <- geoMeta$yday[7]
GranuleID <- geoMeta$GranuleID[7]
StartDateTime <- geoMeta$StartDateTime[7]
ArchiveSet <- geoMeta$ArchiveSet[7]
OrbitNumber <- geoMeta$OrbitNumber[7]
DayNightFlag <- geoMeta$DayNightFlag[7]
EastBoundingCoord <- geoMeta$EastBoundingCoord[7]
NorthBoundingCoord <- geoMeta$NorthBoundingCoord[7]
SouthBoundingCoord <- geoMeta$SouthBoundingCoord[7]
WestBoundingCoord <- geoMeta$WestBoundingCoord[7]
GRingLongitude1 <- geoMeta$GRingLongitude1[7]
GRingLongitude2 <- geoMeta$GRingLongitude2[7]
GRingLongitude3 <- geoMeta$GRingLongitude3[7]
GRingLongitude4 <- geoMeta$GRingLongitude4[7]
GRingLatitude1 <- geoMeta$GRingLatitude1[7]
GRingLatitude2 <- geoMeta$GRingLatitude2[7]
GRingLatitude3 <- geoMeta$GRingLatitude3[7]
GRingLatitude4 <- geoMeta$GRingLatitude4[7]






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