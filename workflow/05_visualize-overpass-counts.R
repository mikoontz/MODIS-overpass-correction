# Animate the monthly overpass counts

library(tidyverse)
library(raster)

# raster_template <- raster::raster("data/data_raw/grid_2_5_degree_vars_modis_D_AFC_num_April_2001.tif")
# raster_template <- raster::raster("data/data_raw/grid_0_25_degree_vars_modis_D_AFC_num_April_2001.tif")
raster_template <- raster::raster("data/data_raw/AGG_1deg_SOLARELEV_nocorn_perc_day_total.tif")

# Get all the rasters to local disk (shouldn't take too long, as they are pretty small)
system2(command = "aws", 
        args = paste0("s3 sync s3://earthlab-mkoontz/MODIS-overpass-counts_", res(raster_template)[1], "/ data/data_output/MODIS-overpass-counts_", res(raster_template)[1], "/"))

raster_files <- 
  list.files(paste0("data/data_output/MODIS-overpass-counts_", res(raster_template)[1], "/"), full.names = ) %>% 
  enframe(name = NULL) %>% 
  setNames("filename") %>% 
  dplyr::mutate(filepath = file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1]), filename),
                daynight = ifelse(grepl(pattern = "day", x = filename), yes = "day", no = "night"), 
                yearmonth = substr(filename, start = 1, stop = 7)) %>% 
  tidyr::separate(col = yearmonth, into = c("year", "month"), sep = "-",remove = FALSE)


rasters_as_df <- 
  raster_files %>% 
  dplyr::filter(year == 2015) %>% 
  dplyr::mutate(r = map(filepath, .f = function(path) {
    
    r <- 
      raster::raster(path) %>% 
      as.data.frame(xy = TRUE) %>% 
      setNames(c("lon", "lat", "overpass_count"))
    
    return(r)
    })) %>% 
  tidyr::unnest(cols = r)

ggplot(rasters_as_df,
       aes(x = lon, y = lat, color = overpass_count, fill = overpass_count)) +
  geom_raster() +
  facet_grid(rows = vars(yearmonth), cols = vars(daynight)) +
  scale_fill_viridis_c()
