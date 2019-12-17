library(dplyr)
library(sf)

# Set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))
print(getwd())
dir <- "/Users/curt.whitmire/Documents/GIS/IEA/FishingEffort"

# import feature classes
gdb1 <- paste0(dir, "/", "IEA_FishingEffort_through2018.gdb")
tl <- sf::st_read(dsn=gdb1, layer = "LB0218_TL_Final")
# st_write(tl, dsn=gdb1, layer="LB0218_TL_Final_copy", driver="FileGDB")
str(tl) # returns column types

tl <- tl %>% 
  mutate(TowDate = as.Date(TOWDATE, format="%d-%b-%Y"),
         TowYr = as.numeric(format(TowDate, "%Y"))
         )

minYr = min(tl$TowYr)
maxYr = max(tl$TowYr)

# Create columns with 5-yr intervals
library(lubridate)

tl.summ <- tl %>% 
  filter(TowDate>=as.Date('2002-01-01')) %>% 
  # group_by(by5=cut(TowDate, "5 year")) # create column with 5-yr interval using cut
  # group_by(by5=cut(TowDate, "5 year")) # create column with 5-yr interval using lubridate interval
  mutate(y02_06="y02_06") # create 5-yr fields the old-fashined way, using filter and mutate


# project haul data to TM spatial object
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
tlTM <- st_as_sf(x = tl,                         
                coords = c("longitude_dd", "latitude_dd"),
                crs = projcrs)

