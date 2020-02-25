library(sf)
library(dplyr)

sessionInfo()
sf::sf_extSoftVersion()
#  GEOS           GDAL         proj.4 GDAL_with_GEOS     USE_PROJ_H 
# "3.8.0"        "3.0.4"        "6.3.0"         "true"         "true" 

# Example in ?st_wrap_dateline() works!
(ls1 <- st_sfc(st_linestring(rbind(c(-179,0), c(179,0))), crs = 4326))
(ls1.wrapped <- st_wrap_dateline(ls1))

# Doesn't work
(ls2 <- st_sfc(st_linestring(rbind(c(-100,0), c(100,0))), crs = 4326))
(ls2.wrapped <- st_wrap_dateline(ls2))

# Doesn't work
(ls3 <- st_sfc(st_linestring(rbind(c(-100,0), c(100,0))), crs = 4326))
(ls3.wrapped <- st_wrap_dateline(ls3, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=160")))

# Works!
(ls4 <- st_sfc(st_linestring(rbind(c(-100,0), c(100,0))), crs = 4326))
(ls4.wrapped <- st_wrap_dateline(ls4, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=160.000001")))

# Still works!
(ls5 <- st_sfc(st_linestring(rbind(c(-100,0), c(100,0))), crs = 4326))
(ls5.wrapped <- st_wrap_dateline(ls5, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=359.9999")))

# Doesn't work
(ls6 <- st_sfc(st_linestring(rbind(c(-100,0), c(100,0))), crs = 4326))
(ls6.wrapped <- st_wrap_dateline(ls6, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=360")))


###

pt1 <- c(104.50757, 70.95887)
pt2 <- c(-68.18395, 87.56392)

(180 - pt1[1]) + (pt2[1] + 180)

ls <- st_sfc(st_linestring(rbind(pt1, pt2)), crs = 4326)
(ls.wrapped <- st_wrap_dateline(x = ls, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=359.999")))

plot(ls, axes = TRUE)
plot(st_sfc(st_point(pt1), crs = 4326), add = TRUE, col = "black", pch = 19, cex = 2)
plot(st_sfc(st_point(pt2), crs = 4326), add = TRUE, col = "red", pch = 19, cex = 2)
