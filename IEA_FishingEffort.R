# Script:   IEA_FishingEffort.R
# Author:   Curt Whitmire (NOAA Fisheries)
# Date:     10 Jan 2021
# Purpose:  Summarize fishing effort, represented by line features, into
#           Cartesian grids of specified size.
#           Summarize by 1-yr and 5-yr increments

# Specify the required packages
packages = c("tictoc", "tidyr", "tidyselect", "dplyr", "sf", "zoo", "rstudioapi")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed and loaded
#function provided by Vikram Baliga, July 19, 2015
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Verify all required packages are loaded
search()

tic.clearlog()
tic("total") # Start the clock
tic("data preparation") # Start the subclock

# Set working directory
dir <- "/Users/curt.whitmire/Documents/GIS/IEA/FishingEffort"
# current_path <- getActiveDocumentContext()$path
# setwd(dirname(current_path))
# dir <- dirname(current_path)

# Import feature classes and set input variables
gdb <- paste0(dir, "/", "IEA_FishingEffort_through2019.gdb")

# Check if feature class already exists in memory
if (exists("tlTMfc")) {
} else {
  tlTMfc <- sf::st_read(dsn=gdb, layer = "LB_TL_0219_TM")
}
st_crs(tlTMfc) # Verify projection is "WGS_1984_Transverse_Mercator"
  
# library(smoothr) # for densifying unprojected lines; !!probably best to project in ArcGIS to ensure consistent projection properties
# tl_dd <- sf::st_read(dsn=gdb, layer = "LB_TL_0219") # Read in unprojected line fc
# tl_dens <- smoothr::densify(tl_dd, max_distance = 0.001) # Add vertices to preserve shape during reprojection
# projcrs <- "+proj=tmerc +lat_0=31.96 +lon_0=-121.6 +k=1 +x_0=390000 +y_0=0 +datum=WGS84 +units=m +no_defs"
# tlTMfc <- st_transform(x=tl_dens, crs=projcrs) # Project haul data to TM spatial object
# st_crs(tlTMfc) # Check projection status

grd <- sf::st_read(dsn=gdb, layer = "grid20km_ETsquare") # Include as user input
grdSz <- substring(colnames(grd)[2], 8)

# Filter for desired year range and Update attribute table
yrSrt <- 2010 # Include as user input
yrEnd <- 2019 # Include as user input

tlTM <- tlTMfc %>% 
  mutate(
    TowDate = as.Date(TOWDATE, format="%d-%b-%Y"),
    TowYr = as.numeric(format(TowDate, "%Y")),
    origLen = st_length(tlTMfc)
    ) %>% 
  filter(TowYr >= yrSrt, TowYr <= yrEnd, DUR > 0, Shape_Length > 0)

# Create summary table by Year
summ1 <- as.data.frame(tlTM) %>% 
  group_by(TowYr) %>% 
  summarise(sumLen = sum(origLen), # first check for negative values
            sumDur = sum(DUR), # first check for negative values
            cntTow = n_distinct(HAUL_ID),
            cntVes = n_distinct(DRVID)
            )

toc(log = TRUE, quiet = TRUE) # End subclock
log.txt <- tic.log(format = TRUE)
tic("data geoprocessing") # Start the subclock

# Calculate intersection of towlines within Grid Cells
isct <- st_intersection(tlTM, grd) %>% 
  mutate(partLen = st_length(.),
         propLen = partLen/origLen,
         partDur = propLen*DUR) %>% 
  select(CellID, CentroidLonDD, CentroidLatDD, TowYr, HAUL_ID, DRVID, DRVIDINT, partLen, partDur)

toc(log = TRUE, quiet = TRUE) # End subclock
log.txt <- tic.log(format = TRUE)
tic("data summarization") # Start the subclock

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
  mutate(CellSz = grdSz) %>% 
  select(CellSz, everything())

# Use tidyr pivot_wider create wide data frames for 1-yr summary
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

# Join 3 tables (Vessel Counts, Cell Lengths, Cell Durations) for 1-yr summary
jn1yr <- pivVes %>% 
  full_join(pivLen, by = "CellID") %>% 
  full_join(pivDur, by = "CellID") %>% 
  mutate(CellSz = grdSz) %>% 
  select(CellSz, CellID, everything())

# Create summary table by 5-yr increments
incr = 5 # Include as user input
fun.prop <- function(x) sum(x > 2) / incr
fun.yrSrt <- function(x){ ifelse(x < (yrSrt+incr), yrSrt, x - incr + 1) }

#### Start of OLD method that summarizes number of years where # unique vessels > 2
# summ4 <- as.data.frame(isct) %>% 
#   group_by(CellID, TowYr) %>% 
#   summarise(sumLen = sum(partLen),
#             sumDur = sum(partDur),
#             cntTow = n_distinct(HAUL_ID),
#             cntVes = n_distinct(DRVID) 
#   ) %>% 
#   mutate(TowYrSrt = fun.yrSrt(TowYr),
#          len5sum = rollapply(sumLen, incr, sum, partial=TRUE, align = 'right'),
#          # len5mean = rollapply(sumLen, incr, mean, partial=TRUE, align = 'right'),
#          dur5sum = rollapply(sumDur, incr, sum, partial=TRUE, align = 'right'),
#          # dur5mean = rollapply(sumDur, incr, mean, partial=TRUE, align = 'right'),
#          tow5sum = rollapply(cntTow, incr, sum, partial=TRUE, align = 'right'),
#          ves5yrs = rollapply(cntVes, incr, FUN = fun.prop, partial=TRUE, align = 'right')
#          ) %>%
#   select(-sumLen, -sumDur, -cntTow, -cntVes)
# # %>% 
# #   filter(TowYr - TowYrSrt >= incr - 1)
# 
# # Use tidyr pivot_wider create wide data frames for 5-yr summary
# pivVes5 <- summ4 %>% 
#   select(-len5sum, -dur5sum, -tow5sum) %>% 
#   pivot_wider(names_from = c(TowYrSrt, TowYr),
#               names_sep = "_",
#               values_from = ves5yrs,
#               names_prefix = "ves") %>% 
#   select(CellID, sort(tidyselect::peek_vars())) %>% 
#   select_if(~!all(is.na(.)))
# 
# pivLen5 <- summ4 %>% 
#   select(-dur5sum, -tow5sum, -ves5yrs) %>% 
#   pivot_wider(names_from = c(TowYrSrt, TowYr),
#               names_sep = "_",
#               values_from = len5sum,
#               names_prefix = "len") %>% 
#   select(CellID, sort(tidyselect::peek_vars())) %>%
#   select_if(~!all(is.na(.)))
# 
# pivDur5 <- summ4 %>% 
#   select(-len5sum, -tow5sum, -ves5yrs) %>% 
#   pivot_wider(names_from = c(TowYrSrt, TowYr),
#               names_sep = "_",
#               values_from = dur5sum,
#               names_prefix = "dur") %>% 
#   select(CellID, sort(tidyselect::peek_vars())) %>% 
#   select_if(~!all(is.na(.)))
#### End of OLD method that summarizes number of years where # unique vessels > 2

#### Start of CURRENT method for calculating # of distinct vessels in 5-yr periods
# from:
# https://stackoverflow.com/questions/41183693/r-rolling-sliding-window-and-distinct-count-in-r-for-sliding-number-of-days
summ5 <- as.data.frame(isct) %>% 
  group_by(CellID, TowYr, DRVID) %>% 
  summarise(sumLen = sum(partLen),
            sumDur = sum(partDur),
            cntTow = n_distinct(HAUL_ID)
            )

# Next two steps create a continuous series and fill in missing combinations of TowYr:CellID
all_combs <- expand.grid(CellID=unique(summ5$CellID),
                         DRVID=unique(summ5$DRVID),
                         TowYr=seq(min(summ5$TowYr), max(summ5$TowYr), by=1)
                         )

df <- merge(summ5, all_combs, by=c('CellID', 'DRVID', 'TowYr'), all=TRUE) %>%
  mutate_if(is.numeric, funs(replace_na(., 0))) %>%
  mutate(DRVID = if_else(cntTow != 0, DRVID, NA_character_)) %>% 
  select(CellID, TowYr, DRVID, sumLen, sumDur, cntTow)
  

# Use rollapply funtion ('zoo' package) to calculate rolling sum and # of unique vessel IDs
# wid <- incr * n_distinct(isct$CellID) # Unique TowYr/CellID
# wid <- n_distinct(isct$TowYr) * n_distinct(isct$CellID) # Unique TowYr/CellID
wid <- incr * n_distinct(isct$DRVID)
# wid <- incr
summ6 <- df %>%
  group_by(CellID) %>%
  arrange(CellID, TowYr, DRVID) %>% 
  mutate(
    TowYrSrt = fun.yrSrt(TowYr),
    len5sum = rollapply(sumLen, width=wid, sum, partial=TRUE, align='right'),
    dur5sum = rollapply(sumDur, width=wid, sum, partial=TRUE, align='right'),
    tow5sum = rollapply(cntTow, width=wid, sum, partial=TRUE, align='right'),
    # ves5sum = rollapply(DRVID, width=wid, FUN=function(x) length(unique(x[!is.na(x)])), partial=TRUE, align='right'),
    ves5sum = rollapply(DRVID, width=wid, FUN=function(x) n_distinct(x, na.rm=TRUE), partial=TRUE, align='right')
    ) %>%
  select(CellID, TowYr, TowYrSrt, ves5sum, tow5sum, len5sum, dur5sum)

# Test filter of a single CellID
# summ6b <- summ6 %>% 
#   filter(CellID == 761)

# Return only the last (right-aligned) row of each unique CellID/5-yr combination, using max # of tows
summ7 <- summ6 %>% 
  group_by(CellID, TowYr) %>% 
  filter(row_number()==n()) %>% # return only last row number in series
  filter(TowYr - TowYrSrt >= incr - 1) # filter out first set of summaries that don't encompass specified incr of years
# %>% 
#   filter(tow5sum == max(tow5sum), !is.na(DRVID)) 

# Use tidyr pivot_wider create wide data frames for 5-yr summary
pivVes5 <- summ7 %>% 
  select(CellID, TowYr, TowYrSrt, ves5sum) %>% 
  pivot_wider(names_from = c(TowYrSrt, TowYr),
              values_from = ves5sum,
              names_prefix = "ves") %>% 
  select(CellID, sort(tidyselect::peek_vars())) %>% 
  arrange(CellID) %>% 
  select_if(~!all(is.na(.)))

pivLen5 <- summ7 %>% 
  select(CellID, TowYr, TowYrSrt, len5sum) %>% 
  pivot_wider(names_from = c(TowYrSrt, TowYr),
              names_sep = "_",
              values_from = len5sum,
              names_prefix = "len") %>% 
  select(CellID, sort(tidyselect::peek_vars())) %>%
  arrange(CellID) %>% 
  select_if(~!all(is.na(.)))

pivDur5 <- summ7 %>% 
  select(CellID, TowYr, TowYrSrt, dur5sum) %>% 
  pivot_wider(names_from = c(TowYrSrt, TowYr),
              names_sep = "_",
              values_from = dur5sum,
              names_prefix = "dur") %>% 
  select(CellID, sort(tidyselect::peek_vars())) %>% 
  arrange(CellID) %>% 
  select_if(~!all(is.na(.)))

# Join 3 tables (Vessel Counts, Cell Lengths, Cell Durations) for 5-yr summary
jn5yr <- pivVes5 %>% 
  full_join(pivLen5, by = "CellID") %>% 
  full_join(pivDur5, by = "CellID") %>% 
  mutate(CellSz = grdSz) %>% 
  select(CellSz, CellID, everything()) %>% 
  filter_at(vars(-CellSz, -CellID), any_vars(!is.na(.)))

# Check output for errors in vessel counts
# Only works if 2nd filter on summ7 above is removed
# fld1 <- paste0("ves", yrSrt, "_", yrSrt+incr-1)
# fld2 <- paste0("ves", yrSrt, "_", yrSrt+incr-2)
# fld3 <- paste0("ves", yrSrt, "_", yrSrt+incr-3)
# fld4 <- paste0("ves", yrSrt, "_", yrSrt+incr-4)
# fld5 <- paste0("ves", yrSrt, "_", yrSrt+incr-5)
# erck5 <- jn5yr %>%
#   filter(fld1 < fld2
#          | fld2 < fld3
#          | fld3 < fld4
#          | fld4 < fld5)

# Import strata and habitat polygon feature class
# strata <- sf::st_read(dsn=gdb, layer = "Strata_Final_TM_WGS84_5") # Include as user input
hab <- sf::st_read(dsn=gdb, layer = "Strata_idenSGH_v4_flat_Bizzarro_dissPhyTyp_SubRegion_HabTyp_multi") # Include as user input

# Clean up topology errors
library(rgeos)
hab_cln <- st_buffer(hab, 0)

# Calculate intersection of towlines within Grid Cells
# Depending on the # of years chosen, this process can take several hours
# Suggest separating processing this using an apply function on each year
isctHab <- st_intersection(tlTM, hab_cln) %>% 
  mutate(partLen = st_length(.),
         propLen = partLen/origLen,
         partDur = propLen*DUR) %>% 
  select(TowYr, HAUL_ID, DRVID, DRVIDINT, partLen, partDur, SubRegion, PhysType, HabType)

# Summarize by TowYr, sub-region, depth strata and habitat type, and tally total effort (km, hrs) and # unique vessels
summ8 <- as.data.frame(isctHab) %>% 
  group_by(TowYr, SubRegion, PhysType, HabType) %>% 
  summarise(Length_km = sum(partLen) / 1000,
            Durat_hr = sum(partDur),
            numVessels = n_distinct(DRVID),
            numHauls = n_distinct(HAUL_ID)
            )

# Output to CSV files for later joining to polygon features
# OR join here and output to FileGDB using st_write
setwd(paste0(dir, "/Data"))
fp <- paste0("LB", substring(yrSrt, 3, 4), substring(yrEnd, 3, 4), "_summ", grdSz)
log <- paste0(fp, "_log.txt")
err_log <- paste0(fp, "_errCk.txt")

write.csv(jn1yr, 
          file = paste0(fp, "_1yr_dat.csv"), 
          na = "0"
          )
write.csv(jn5yr, 
          file = paste0(fp, "_", incr, "yr_dat.csv"), 
          na = "0"
          )
write.csv(summ3, 
          file = paste0(fp, "_1yr_summ.csv"), 
          # na = "0"
          )

write.csv(summ8, 
          file = paste0("LB", substring(yrSrt, 3, 4), substring(yrEnd, 3, 4), "_1yr_hab_dat.csv"), 
          # na = "0"
          )

toc(log = TRUE, quiet = TRUE) # End subclock
log.txt <- tic.log(format = TRUE)
toc(log = TRUE, quiet = TRUE) # End clock
log.txt <- tic.log(format = TRUE)

# Write elapsed code section times to log file
fileConn <- file(log)
writeLines(unlist(log.txt), con = fileConn)
close(fileConn)

# Clear process log
tic.clearlog()

# Delete temporary variables except input lines
rm(list=setdiff(ls(), "tlTMfc"))
