# This script must be run locally on alicja's computer because farm has outdated terra package and I still couldn't get this updated...
# But it runs quite quickly 

load_pkgs <- function(pkg, bioconductor = FALSE) {
  for (p in pkg) {
    library(p, character.only = TRUE)
  }
}
# CRAN packages
cran_pkgs <- c(
  "dplyr", "cluster", "reshape", "reshape2", "stringdist",
  "ggiraph", "e1071", "gridExtra", "colorspace", "purrr",
  "tidyverse", "RColorBrewer", "scales", "kableExtra",
  "knitr", "patchwork", "rnaturalearth", "rnaturalearthdata", 
  "ggplot2", "tidyr", "stringr", "terra", "dismo",
  "parallel", "raster", "ncdf4", "seqinr", "vegan", "reshape2", "remotes",
  "phangorn", "shiny", "sf", "tibble", "forcats", "lubridate", "viridis", "maps"
)
# Load CRAN packages
load_pkgs(cran_pkgs, bioconductor = FALSE)

# Get sample-level data (csv output of 01_subset_data.R)
working_set <- read.csv("/Users/aw43/farm22/100k_paper/output/BIOSCAN_100k_samples_corrected2025-11-20.csv")
# Select needed columns and remove repeated entries
working_set <- working_set %>% dplyr::select(sts_CATCH_LOT, trap_name, sts_latitude, sts_longitude, sts_col_date, day, month, year) %>% unique()

# Prepare file maps for each variable
make_file_map <- function(folder, var) {
  files <- list.files(folder, pattern = "\\.nc$", full.names = TRUE)
  tibble(
    file = files,
    month = as.integer(str_sub(basename(files), -7, -6)), # e.g. 2023
    year = as.integer(str_sub(basename(files), -11, -8)), # e.g. 05
    var = var
  )
}
# 2021-2024 & provisional 2025 until July [these are present in Alicja's dropbox until the terra update happens]
tasmin_map <- make_file_map("~/Dropbox/1_BIOSCAN2024/weather/tasmin_dump", "tasmin") 
tasmax_map <- make_file_map("~/Dropbox/1_BIOSCAN2024/weather/tasmax_dump", "tasmax")
rain_map   <- make_file_map("~/Dropbox/1_BIOSCAN2024/weather/rain_dump",   "rainfall")

# Map the dates using file names
file_map <- bind_rows(tasmin_map, tasmax_map, rain_map)

# Function to extract values for one record and one variable
extract_vals <- function(var, lon, lat, date) {
  pt <- vect(data.frame(lon = lon, lat = lat), crs = "EPSG:4326")
  
  #### SELECT DAYS ####
  # Days of interest (the sampling day and two days prior - can be adjusted here!!!)
  dates_needed <- c(date, date - 1, date - 2)
  
  out <- sapply(dates_needed, function(d) {
    yr <- year(d); mo <- month(d)
    f <- file_map %>% filter(var == !!var, year == yr, month == mo) %>% pull(file)
    if (length(f) == 0) return(NA_real_)
    r <- rast(f)
    dates <- as.Date(time(r))
    idx <- match(d, dates)
    if (is.na(idx)) return(NA_real_)
    pt_grid <- project(pt, crs(r))
    terra::extract(r[[idx]], pt_grid, method = "bilinear")[,2]
  })
  return(out)
}

# Apply
results <- working_set %>%
  rowwise() %>%
  mutate(
    tasmin_vals = list(extract_vals("tasmin", sts_longitude, sts_latitude, as.Date(sts_col_date))),
    tasmax_vals = list(extract_vals("tasmax", sts_longitude, sts_latitude, as.Date(sts_col_date))),
    rain_vals   = list(extract_vals("rainfall", sts_longitude, sts_latitude, as.Date(sts_col_date)))
  ) %>%
  mutate(
    tasmin_sampling = tasmin_vals[[1]], tasmin2 = tasmin_vals[[2]], tasmin3 = tasmin_vals[[3]],
    tasmax_sampling = tasmax_vals[[1]], tasmax2 = tasmax_vals[[2]], tasmax3 = tasmax_vals[[3]],
    rain_sampling   = rain_vals[[1]],   rain2   = rain_vals[[2]],   rain3   = rain_vals[[3]]
  ) %>%
  dplyr::select(-tasmin_vals, -tasmax_vals, -rain_vals) %>%
  ungroup()

# table(is.na(results))

# Filter results and calculate medians 
results_summary <- results %>% dplyr::select(sts_CATCH_LOT, tasmin_sampling, tasmin2, tasmin3, tasmax_sampling, tasmax2, tasmax3, rain_sampling, rain2, rain3) %>%
  rowwise() %>%
  mutate(
    temp_med_sampling = median(c(tasmin_sampling, tasmax_sampling), na.rm = TRUE),
    temp_med_2 = median(c(tasmin2, tasmax2), na.rm = TRUE),
    temp_med_3 = median(c(tasmin3, tasmax3), na.rm = TRUE),
    tasmin_med = median(c(tasmin_sampling, tasmin2, tasmin3), na.rm = TRUE),
    tasmax_med = median(c(tasmax_sampling, tasmax2, tasmax3), na.rm = TRUE),
    rain_med   = median(c(rain2,  rain2), na.rm = TRUE),
    # temp_med = median(c(tasmin_sampling, tasmin2, tasmin3, tasmax_sampling, tasmax2, tasmax3), na.rm = TRUE),
    temp_med = median(c(temp_med_sampling, temp_med_2, temp_med_3), na.rm = TRUE)
  ) %>%
  ungroup()

# Save
outfile <- paste0("/Users/aw43/farm22/100k_paper/output/weather_data_", format(Sys.Date(), "%Y-%m-%d"), ".csv")
write.csv(results_summary, outfile, row.names = FALSE)

table(results_summary$sts_CATCH_LOT %in% working_set$sts_CATCH_LOT)
table(working_set$sts_CATCH_LOT %in% working_set$sts_CATCH_LOT)

message("Saved: ", outfile)
