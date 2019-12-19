ptm <- proc.time() # Start the clock!

library(tidyr)
library(tidyselect)
library(dplyr)
library(sf)
library(smoothr)

# Set working directory
# current_path <- getActiveDocumentContext()$path
# setwd(dirname(current_path))
# print(getwd())
dir <- "/Users/curt.whitmire/Documents/GIS/IEA/FishingEffort"

# Import feature classes
gdb <- paste0(dir, "/", "IEA_FishingEffort_through2018.gdb")
tl_dd <- sf::st_read(dsn=gdb, layer = "LB0218_TL_Final")
tl <- densify(tl_dd, max_distance = 0.001) # add vertices to preserve shape during reprojection
# tlTM <- sf::st_read(dsn=gdb, layer = "LB0218_TL_Final_TM")
# project haul data to TM spatial object
projcrs <- "+proj=tmerc +lat_0=31.96 +lon_0=-121.6 +k=1 +x_0=390000 +y_0=0 +datum=WGS84 +units=m +no_defs"
tlTM <- st_transform(x=tl, crs=projcrs) 

grd <- sf::st_read(dsn=gdb, layer = "grid20km_ETsquare")
grdSz <- substring(colnames(grd)[4],8,9)

# Update attribute table
tlTM <- tlTM %>% 
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

tlTM <- tlTM %>% 
  mutate(
    origLen = st_length(tlTM)
  )

# Calculate intersection of towlines within Grid Cells
# If necessary because entire dataset is too large, calculate iteratively by year and grid cell size

# i10 <- st_intersection(tlTM, grd10)
# ptm <- proc.time() # Start the clock!
isct <- st_intersection(tlTM, grd)
# proc.time() - ptm # Stop the clock
# n_distinct(isct$DRVID)

# Create summary table by CellID, Year, and VesselID
# ptm <- proc.time() # Start the clock!
summ1 <- as.data.frame(isct) %>% 
  # select(HAUL_ID, DRVID, DUR, TowYr, CellID_20km, origLen) %>% 
  mutate(
    partLen = st_length(isct),
    propLen = partLen/origLen,
    propDur = propLen*DUR
  ) %>% 
  group_by(CellID, TowYr, DRVID) %>% 
  summarise(preLen = sum(propLen),
            preDur = sum(propDur),
            preTow = n_distinct(HAUL_ID)
  ) 
# proc.time() - ptm # Stop the clock
# print(paste0("Error Check:  Maximum proportional length = ", max(i20_summ$propLen))) # old ErrCk

# Create summary table by CellID and Year
# Need to add logic to count unique DRVID values and filter by CellID's with values <3
# use dplyr::n_distinct()
summ2 <- summ1 %>% 
  group_by(CellID, TowYr) %>% 
  summarise(totLen = sum(preLen),
            totDur = sum(preDur),
            cntVes = n_distinct(DRVID),
            cntTow = sum(preTow)
            )

# Use tidyr pivot_wider create wide data frame
yrSrt <- 2002 # Set up for user input
yrEnd <- 2018 # Set up for user input

pivLen <- summ2 %>% 
  filter(TowYr >= yrSrt, TowYr <= yrEnd) %>% 
  select(-totDur, -cntVes, -cntTow) %>% 
  pivot_wider(names_from = TowYr, 
              values_from = totLen,
              names_prefix = "len") %>% 
  select(CellID, sort(tidyselect::peek_vars()))

pivDur <- summ2 %>% 
  filter(TowYr >= yrSrt, TowYr <= yrEnd) %>% 
  select(-totLen, -cntVes, -cntTow) %>% 
  pivot_wider(names_from = TowYr, 
              values_from = totDur,
              names_prefix = "dur") %>% 
  select(CellID, sort(tidyselect::peek_vars()))

pivVes <- summ2 %>% 
  filter(TowYr >= yrSrt, TowYr <= yrEnd) %>% 
  select(-totLen, -totDur, -cntTow) %>% 
  pivot_wider(names_from = TowYr, 
              values_from = cntVes,
              names_prefix = "ves") %>% 
  select(CellID, sort(tidyselect::peek_vars()))

# Join 3 tables (Vessel Counts, Cell Lengths, Cell Durations)
jn <- pivVes %>% 
  full_join(pivLen, by = "CellID") %>% 
  full_join(pivDur, by = "CellID") %>% 
  mutate(CellSz = paste0(grdSz,"-km")) %>% 
  select(CellSz, CellID, everything())

# Output to CSV file for joining to Polygon features OR join here and output to FileGDB using st_write
setwd(paste0(dir, "/Data"))
write.csv(jn, 
          file = paste0("LB", substring(yrSrt,3,4), substring(yrEnd,3,4), "_iden", grdSz, "km_1yr.csv"), 
          na = "0")

proc.time() - ptm # Stop the clock
