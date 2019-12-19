library(tidyr)
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
# str(tl) # returns column types
# grd10 <- sf::st_read(dsn=gdb1, layer = "grid10km_ETsquare")
grd20 <- sf::st_read(dsn=gdb1, layer = "grid20km_ETsquare")

tl <- tl %>% 
  mutate(
    TowDate = as.Date(TOWDATE, format="%d-%b-%Y"),
    TowYr = as.numeric(format(TowDate, "%Y"))
    )

# minYr = min(tl$TowYr)
# maxYr = max(tl$TowYr)

# # Create columns with 5-yr intervals
# library(lubridate)
# tl.summ <- tl %>% 
#   filter(TowDate>=as.Date('2002-01-01')) %>% 
#   # group_by(by5=cut(TowDate, "5 year")) # create column with 5-yr interval using cut
#   # group_by(by5=cut(TowDate, "5 year")) # create column with 5-yr interval using lubridate interval
#   mutate(y02_06="y02_06") # create 5-yr fields the old-fashined way, using filter and mutate

# project haul data to TM spatial object
projcrs <- "+proj=tmerc +lat_0=31.96 +lon_0=-121.6 +k=1 +x_0=390000 +y_0=0 +datum=WGS84 +units=m +no_defs"
tlTM <- st_transform(x=tl, crs=projcrs)
# st_crs(tlTM)
tlTM <- tlTM %>% 
  mutate(
    origLen = st_length(tlTM)
  )

# Calculate intersection of towlines within Grid Cells
# If necessary because entire dataset is too large, calculate iteratively by year and grid cell size

# i10 <- st_intersection(tlTM, grd10)
i20 <- st_intersection(tlTM, grd20)
n_distinct(i20$DRVID)

# Create summary table by CellID and Year
i20_summ <- as.data.frame(i20) %>% 
  # select(HAUL_ID, DRVID, DUR, TowYr, CellID_20km, origLen) %>% 
  mutate(
    partLen = st_length(i20),
    propLen = partLen/origLen,
    propDur = propLen*DUR
  ) %>% 
  group_by(CellID_20km, TowYr) %>% 
  summarise(totLen=sum(propLen),
            totDur=sum(propDur)
  ) 
# print(paste0("Error Check:  Maximum proportional length = ", max(i20_summ$propLen))) # old ErrCk

# Need to add logic to count unique DRVID values and filter by CellID's with values <3
# use dplyr::n_distinct()

# Use tidyr pivot_wider create wide data frame
i20pivLen <- i20_summ %>% 
  filter(TowYr!=2001) %>% 
  select(-totDur) %>% 
  pivot_wider(names_from = TowYr, 
              values_from = totLen,
              names_prefix = "len")


i20pivDur <- i20_summ %>% 
  filter(TowYr!=2001) %>% 
  select(-totLen) %>% 
  pivot_wider(names_from = TowYr, 
              values_from = totDur,
              names_prefix = "dur")

# Join 3 tables (Vessel Counts, Cell Lengths, Cell Durations)??

# Output to CSV file for joining to Polygon features OR join here and output to FileGDB using st_write

