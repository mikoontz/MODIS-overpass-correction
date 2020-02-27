# Combine overpass rasters across satellite into desired temporal aggregations

# day overpass count per year/month 
# night overpass count per year/month
# day overpass count per year
# night overpass count per year
# day overpass count for the whole study
# night overpass count for the whole study

library(tidyverse)
library(raster)

# create new directory to house the analysis ready data
if(!dir.exists(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/")))) {
  dir.create(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/")))
}

# What is the raster template for the rasters being combined?
raster_template <- raster::raster("data/data_raw/grid_2_5_degree_vars_modis_D_AFC_num_April_2001.tif")

# Get all the rasters to local disk (shouldn't take too long, as they are pretty small)
system2(command = "aws", 
        args = paste0("s3 sync s3://earthlab-mkoontz/MODIS-overpass-counts_", res(raster_template)[1], "/ data/data_output/MODIS-overpass-counts_", res(raster_template)[1], "/"))

# all of the rasters on disk, plus some extra metadata extracted from the filename
all_rasters <- 
  list.files(paste0("data/data_output/MODIS-overpass-counts_", res(raster_template)[1], "/")) %>% 
  enframe(name = NULL) %>% 
  setNames("filename") %>% 
  tidyr::separate(col = filename, into = c("year_month", "satellite", "daynight_overpass_count", "resolution"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(year = substr(year_month, start = 1, stop = 4),
                resolution = substr(resolution, start = 1, stop = nchar(resolution) - 4)) %>% 
  dplyr::mutate(daynight = ifelse(str_detect(daynight_overpass_count, pattern = "day"), yes = "day", no = "night"))


# split the rasters by various means --------------------------------------

# Split by year/month/daynight combinations -------------------------------

year_month_combos <-
  all_rasters %>% 
  group_by(year_month, daynight) %>% 
  group_split()

# Create directory to house the year/month/daynight combinations
if(!dir.exists(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/monthly")))) {
  dir.create(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/monthly")))
}

# work through the year/month/daynight combos -----------------------------

plan(multiprocess)

future_map(year_month_combos, .f = function(this_year_month_combo) {
  
  rasters_to_read <- file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1]), this_year_month_combo$filename)
  
  # Combine the rasters in this group by summing them across the stack
  combined_rasters <- map(rasters_to_read, raster::raster) %>% raster::brick() %>% sum()
  
  combined_filename <- paste0(unique(this_year_month_combo$year_month), "_", unique(this_year_month_combo$daynight), "_overpass-count.tif")
  
  combined_filepath <- file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/monthly/", combined_filename))
  
  raster::writeRaster(x = combined_rasters, filename = combined_filepath)
  
  system2(command = "aws", args = paste0('s3 cp ', combined_filepath, ' s3://earthlab-mkoontz/MODIS-overpass-counts_', res(raster_template)[1], '_analysis-ready/monthly/', combined_filename))
})


# Annual aggregations -----------------------------------------------------

# We can use the rasters already created on a monthly timestep to create the ones on an
# annual time step. This should save quite a bit of processing time. Since we are summing,
# the results will be exactly the same as if we started with all the individual
# year/month/satellite/daynight rasters first!
year_combos <-
  list.files(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/monthly"))) %>% 
  enframe(name = NULL) %>% 
  setNames("filename") %>% 
  dplyr::mutate(year = substr(year_month, start = 1, stop = 4),
                resolution = substr(resolution, start = 1, stop = nchar(resolution) - 4)) %>% 
  dplyr::mutate(daynight = ifelse(str_detect(daynight_overpass_count, pattern = "day"), yes = "day", no = "night")) %>% 
  group_by(year, daynight) %>% 
  group_split()

# Create directory to house the year/daynight combinations
if(!dir.exists(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/annual")))) {
  dir.create(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/annual")))
}


# Aggregate across teh whole record ---------------------------------------

# Again, we can use the previously-combined rasters per year to add together and
# make the final raster for the whole record.
daynight_combos <-
  list.files(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/annual"))) %>% 
  enframe(name = NULL) %>% 
  setNames("filename") %>% 
  dplyr::mutate(year = substr(year_month, start = 1, stop = 4),
                resolution = substr(resolution, start = 1, stop = nchar(resolution) - 4)) %>% 
  dplyr::mutate(daynight = ifelse(str_detect(daynight_overpass_count, pattern = "day"), yes = "day", no = "night")) %>% 
  group_by(daynight) %>% 
  group_split()

# Create directory to house the whole record's day/night overpass count
if(!dir.exists(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/whole-record")))) {
  dir.create(file.path("data", "data_output", paste0("MODIS-overpass-counts_", res(raster_template)[1], "_analysis-ready/whole-record")))
}


