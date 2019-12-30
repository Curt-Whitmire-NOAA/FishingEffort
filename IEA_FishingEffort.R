library(tictoc)
tic.clearlog()
tic("total") # Start the clock!

# Load required libraries
library(tidyr)
library(tidyselect)
library(dplyr)
library(sf)
library(zoo)

tic("data preparation") # Start the clock!
# Set working directory
dir <- "/Users/curt.whitmire/Documents/GIS/IEA/FishingEffort"
# library(rstudioapi)
# current_path <- getActiveDocumentContext()$path
# setwd(dirname(current_path))

# Import feature classes
gdb <- paste0(dir, "/", "IEA_FishingEffort_through2018.gdb")

if (exists("tlTMfc")) {
} else {
  tlTMfc <- sf::st_read(dsn=gdb, layer = "LB0218_TL_Final_TM")
}
  
# library(smoothr) # for densifying unprojected lines
# tl_dd <- sf::st_read(dsn=gdb, layer = "LB0218_TL_Final") # Read in unprojected line fc
# tl_dens <- smoothr::densify(tl_dd, max_distance = 0.001) # Add vertices to preserve shape during reprojection
# projcrs <- "+proj=tmerc +lat_0=31.96 +lon_0=-121.6 +k=1 +x_0=390000 +y_0=0 +datum=WGS84 +units=m +no_defs"
# tlTMfc <- st_transform(x=tl_dens, crs=projcrs) # Project haul data to TM spatial object

grd <- sf::st_read(dsn=gdb, layer = "grid20km_ETsquare") # Set up for user input
grdSz <- substring(colnames(grd)[4],8,9)

# Filter for desired year range and Update attribute table
yrSrt <- 2013 # Set up for user input
yrEnd <- 2018 # Set up for user input

tlTM <- tlTMfc %>% 
  mutate(
    TowDate = as.Date(TOWDATE, format="%d-%b-%Y"),
    TowYr = as.numeric(format(TowDate, "%Y")),
    origLen = st_length(tlTMfc)
    ) %>% 
  filter(TowYr >= yrSrt, TowYr <= yrEnd)

# # Create columns with 5-yr intervals
# library(lubridate)
# tl.summ <- tl %>% 
#   filter(TowDate>=as.Date('2002-01-01')) %>% 
#   # group_by(by5=cut(TowDate, "5 year")) # create column with 5-yr interval using cut
#   # group_by(by5=cut(TowDate, "5 year")) # create column with 5-yr interval using lubridate interval
#   mutate(y02_06="y02_06") # create 5-yr fields the old-fashined way, using filter and mutate

# Create summary table by Year
summ1 <- as.data.frame(tlTM) %>% 
  group_by(TowYr) %>% 
  summarise(totLen = sum(origLen),
            totDur = sum(DUR),
            cntTow = n_distinct(HAUL_ID),
            numVes = n_distinct(DRVID)
            )

toc(log = TRUE, quiet = TRUE) # End subclock
log.txt <- tic.log(format = TRUE)

tic("data geoprocessing") # Start the clock!
# Calculate intersection of towlines within Grid Cells
isct <- st_intersection(tlTM, grd) %>% 
  mutate(partLen = st_length(.),
         propLen = partLen/origLen,
         partDur = propLen*DUR) %>% 
  select(CellID, CentroidLonDD, CentroidLatDD, TowYr, HAUL_ID, DRVID, partLen, partDur)

toc(log = TRUE, quiet = TRUE) # End subclock
log.txt <- tic.log(format = TRUE)

tic("data summarization") # Start the clock!
# Create summary table by CellID and Year
summ2 <- as.data.frame(isct) %>% 
  group_by(CellID, TowYr) %>% 
  summarise(sumLen = sum(partLen),
            sumDur = sum(partDur),
            cntTow = n_distinct(HAUL_ID),
            cntVes = n_distinct(DRVID)
            )

# Create summary table by Year
summ3 <- as.data.frame(isct) %>% 
  group_by(TowYr) %>% 
  summarise(sumLen = sum(partLen),
            sumDur = sum(partDur),
            cntTow = n_distinct(HAUL_ID),
            cntVes = n_distinct(DRVID),
            cntCell = n_distinct(CellID)
            ) %>% 
  mutate(CellSz = paste0(grdSz,"-km")) %>% 
  select(CellSz, everything())

# Create summary table by 5-yr increments
incr = 5 # Include in function as user input
summ4 <- as.data.frame(isct) %>% 
  group_by(CellID, TowYr) %>% 
  summarise(sumLen = sum(partLen),
            sumDur = sum(partDur),
            cntTow = n_distinct(HAUL_ID),
            cntVes = n_distinct(DRVID) 
            ) %>% 
  mutate(len5sum=rollsum(sumLen,incr,align='right',fill=NA),
         len5mean=rollmean(sumLen,incr,align='right',fill=NA),
         tow5sum=rollsum(cntTow,incr,align='right',fill=NA),
         ves5min=rollapply(cntVes,incr,min,align='right',fill=NA)
  )

# Use tidyr pivot_wider create wide data frame
pivLen <- summ2 %>% 
  select(-sumDur, -cntTow, -cntVes) %>% 
  pivot_wider(names_from = TowYr, 
              values_from = sumLen,
              names_prefix = "len") %>% 
  select(CellID, sort(tidyselect::peek_vars()))

pivDur <- summ2 %>% 
  select(-sumLen, -cntTow, -cntVes) %>% 
  pivot_wider(names_from = TowYr, 
              values_from = sumDur,
              names_prefix = "dur") %>% 
  select(CellID, sort(tidyselect::peek_vars()))

pivVes <- summ2 %>% 
  select(-sumLen, -sumDur, -cntTow) %>% 
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

toc(log = TRUE, quiet = TRUE) # End subclock
log.txt <- tic.log(format = TRUE)

# Output to CSV file for joining to Polygon features OR join here and output to FileGDB using st_write
setwd(paste0(dir, "/Data"))
fp <- paste0("LB", substring(yrSrt,3,4), substring(yrEnd,3,4), "_iden", grdSz, "km_1yr")

write.csv(jn, 
          file = paste0(fp, ".csv"), 
          na = "0")
write.csv(summ3, 
          file = paste0(fp, "_summ.csv"), 
          na = "0")

toc(log = TRUE, quiet = TRUE) # End clock
log.txt <- tic.log(format = TRUE)

# Write elapsed code section times to log file
log <- paste0(fp, "_log.txt")
fileConn <- file(log)
writeLines(unlist(log.txt), con = fileConn)
close(fileConn)

rm(grd,isct,jn,log.txt,pivDur,pivLen,pivVes,summ2,summ3,tlTM,fileConn,fp,grdSz,log,yrEnd,yrSrt)

tic.clearlog()
