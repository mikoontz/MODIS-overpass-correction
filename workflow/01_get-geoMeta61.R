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

# NOTE: Must use latest versions of GDAL, GEOS, and PROJ for this to work!
# This SO answer was helpful for making sure {sf} linked to the latest
# versions: https://stackoverflow.com/questions/44973639/trouble-installing-sf-due-to-gdal

aqua_years <- 2020
terra_years <- 2020
geoMeta_source <- "laads-daac" # could also be "laads-daac" or "local"

if(geoMeta_source != "local") {
  if(geoMeta_source == "laads-daac") {
    api_keys <- read_csv("data/data_raw/LAADS-DAAC_api-keys.csv")
    
    if(!dir.exists(paths = file.path("data", "data_output", "MODIS-geoMeta61", "TERRA"))) {
      dir.create(path = file.path("data", "data_output", "MODIS-geoMeta61", "TERRA"))
    }
    
    lapply(terra_years, FUN = function(this_year) {
      system2(command = "wget", args = paste0('-e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov/archive/geoMeta/61/TERRA/', this_year, '/" --header "Authorization: Bearer ', dplyr::filter(api_keys, satellite == "Terra") %>% dplyr::pull(key), '" -P data/data_output/MODIS-geoMeta61/'))
    })
    
    
    if(!dir.exists(paths = file.path("data", "data_output", "MODIS-geoMeta61", "AQUA"))) {
      dir.create(path = file.path("data", "data_output", "MODIS-geoMeta61", "AQUA"))
    }
    
    lapply(aqua_years, FUN = function(this_year) {
      system2(command = "wget", args = paste0('-e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov/archive/geoMeta/61/AQUA/', this_year, '/" --header "Authorization: Bearer ',  dplyr::filter(api_keys, satellite == "Aqua") %>% dplyr::pull(key), '" -P data/data_output/MODIS-geoMeta61/'))
    })
    
  } else if(geoMeta_source == "s3") {
    system2(command = "aws", args = paste0('s3 sync s3://earthlab-mkoontz/MODIS-geoMeta61/ data/data_output/MODIS-geoMeta61/'))
  }
}

system2(command = "aws", args = paste0('s3 sync data/data_output/MODIS-geoMeta61/ s3://earthlab-mkoontz/MODIS-geoMeta61/'))
