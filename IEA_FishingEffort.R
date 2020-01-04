# Script:   IEA_FishingEffort.R
# Author:   Curt Whitmire (NOAA Fisheries)
# Date:     2 Jan 2020
# Purpose:  Summarize fishing effort, represented by line features, into
#           Cartesian grids of specified size

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

titic.clearlog()
tic("total") # Start the clock
tic("data preparation") # Start the subclock

# Set working directory
dir <- "/Users/curt.whitmire/Documents/GIS/IEA/FishingEffort"
# current_path <- getActiveDocumentContext()$path
# setwd(dirname(current_path))

# Import feature classes and set input variables
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

grd <- sf::st_read(dsn=gdb, layer = "grid20km_ETsquare") # Include as user input
grdSz <- substring(colnames(grd)[4], 8)

# Filter for desired year range and Update attribute table
yrSrt <- 2002 # Include as user input
yrEnd <- 2018 # Include as user input

tlTM <- tlTMfc %>% 
  mutate(
    TowDate = as.Date(TOWDATE, format="%d-%b-%Y"),
    TowYr = as.numeric(format(TowDate, "%Y")),
    origLen = st_length(tlTMfc)
    ) %>% 
  filter(TowYr >= yrSrt, TowYr <= yrEnd)

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
tic("data geoprocessing") # Start the subclock

# Calculate intersection of towlines within Grid Cells
isct <- st_intersection(tlTM, grd) %>% 
  mutate(partLen = st_length(.),
         propLen = partLen/origLen,
         partDur = propLen*DUR) %>% 
  select(CellID, CentroidLonDD, CentroidLatDD, TowYr, HAUL_ID, DRVID, partLen, partDur)

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
fun.cond1 <- function(x){ ifelse(x < (yrSrt+incr), yrSrt, x - incr + 1) }

summ4 <- as.data.frame(isct) %>% 
  group_by(CellID, TowYr) %>% 
  summarise(sumLen = sum(partLen),
            sumDur = sum(partDur),
            cntTow = n_distinct(HAUL_ID),
            cntVes = n_distinct(DRVID) 
  ) %>% 
  mutate(TowYrSrt = fun.cond1(TowYr),
         len5sum = rollapply(sumLen, incr, sum, partial=TRUE, align = 'right', fill = NA),
         len5mean = rollapply(sumLen, incr, mean, partial=TRUE, align = 'right', fill = NA),
         dur5sum = rollapply(sumDur, incr, sum, partial=TRUE, align = 'right', fill = NA),
         dur5mean = rollapply(sumDur, incr, mean, partial=TRUE, align = 'right', fill = NA),
         tow5sum = rollapply(cntTow, incr, sum, partial=TRUE, align = 'right', fill = NA),
         ves5yrs = rollapply(cntVes, incr, FUN = fun.prop, partial=TRUE, align = 'right', fill = NA)
         ) %>% 
  select(-sumLen, -sumDur, -cntTow, -cntVes)

# Use tidyr pivot_wider create wide data frames for 5-yr summary
pivLen5 <- summ4 %>% 
  select(-len5mean, -dur5sum, -dur5mean, -tow5sum, -ves5yrs) %>% 
  pivot_wider(names_from = c(TowYrSrt, TowYr),
              names_sep = "_",
              values_from = len5sum,
              names_prefix = "len") %>% 
  select(CellID, sort(tidyselect::peek_vars())) %>%
  select_if(~!all(is.na(.)))

pivDur5 <- summ4 %>% 
  select(-len5sum, -len5mean, -dur5mean, -tow5sum, -ves5yrs) %>% 
  pivot_wider(names_from = c(TowYrSrt, TowYr),
              names_sep = "_",
              values_from = dur5sum,
              names_prefix = "dur") %>% 
  select(CellID, sort(tidyselect::peek_vars())) %>% 
  select_if(~!all(is.na(.)))

pivVes5 <- summ4 %>% 
  select(-len5sum, -len5mean, -dur5sum, -dur5mean, -tow5sum) %>% 
  pivot_wider(names_from = c(TowYrSrt, TowYr),
              names_sep = "_",
              values_from = ves5yrs,
              names_prefix = "ves") %>% 
  select(CellID, sort(tidyselect::peek_vars())) %>% 
  select_if(~!all(is.na(.)))

# Join 3 tables (Vessel Counts, Cell Lengths, Cell Durations) for 5-yr summary
fun.cond2 <- function(x){ as.numeric(substr(y,9,12)) - as.numeric(substr(y,4,7)) >= (incr-1) }

jn5yr <- pivVes5 %>% 
  full_join(pivLen5, by = "CellID") %>% 
  full_join(pivDur5, by = "CellID") %>% 
  mutate(CellSz = grdSz) %>% 
  select(CellSz, CellID, everything()) %>% 
  filter_at(vars(-CellSz, -CellID), any_vars(!is.na(.))) 

# %>% 
#   select_if(fun.cond2)

# Output to CSV files for later joining to polygon features
# OR join here and output to FileGDB using st_write
setwd(paste0(dir, "/Data"))
fp <- paste0("LB", substring(yrSrt, 3, 4), substring(yrEnd, 3, 4), "_summ", grdSz)

write.csv(jn1yr, 
          file = paste0(fp, "_1yr_dat.csv"), 
          na = "0")
write.csv(jn5yr, 
          file = paste0(fp, "_", incr, "yr_dat.csv"), 
          na = "0")
write.csv(summ3, 
          file = paste0(fp, "_1yr_summ.csv"), 
          na = "0")

toc(log = TRUE, quiet = TRUE) # End subclock
log.txt <- tic.log(format = TRUE)
toc(log = TRUE, quiet = TRUE) # End clock
log.txt <- tic.log(format = TRUE)

# Write elapsed code section times to log file
log <- paste0(fp, "_log.txt")
fileConn <- file(log)
writeLines(unlist(log.txt), con = fileConn)
close(fileConn)

# Clear process log
tic.clearlog()

# Delete temporary variables except input lines
rm(list=setdiff(ls(), "tlTMfc"))
