# This script takes the output of:
# ~/bioscan/processing/code/manifest_fetch.sh
# And subsets the data for further processing 

# Load required libraries; All libraries should be automatically installed in the environment
load_pkgs <- function(pkg, bioconductor = FALSE) {
  for (p in pkg) {
    library(p, character.only = TRUE)
  }
}
# CRAN packages
cran_pkgs <- c(
  "dplyr", "cluster", "reshape", "reshape2", "stringdist", "pander",
  "ggiraph", "e1071", "gridExtra", "colorspace", "purrr",
  "tidyverse", "RColorBrewer", "scales", "kableExtra",
  "here", "knitr", "patchwork", "rnaturalearth", "rnaturalearthdata", 
  "ggplot2", "tidyr", "stringr", "terra", "dismo",
  "parallel", "bigmemory", "raster", "ncdf4", "seqinr", "vegan", "reshape2", "remotes",
  "phangorn", "shiny", "sf", "textshape", "tibble", "forcats", "lubridate", "viridis", "maps"
)
# Bioconductor packages
bioconductor_pkgs <- c(
  "biomaRt", "Biostrings", "msa", "ape"
)
# Load CRAN packages
load_pkgs(cran_pkgs, bioconductor = FALSE)
# Load Bioconductor packages
load_pkgs(bioconductor_pkgs, bioconductor = TRUE)

# Load the most up to date file created in 01_subset_data.R

pattern  <- "^BIOSCAN_100k_samples_(\\d{4}-\\d{2}-\\d{2})\\.csv$"   
dir_in   <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/"
files <- list.files(dir_in, pattern = pattern, full.names = TRUE)
if (length(files) == 0) {
  stop("No dated BIOSCAN_100k_samples_*.csv files found in ", dir_in)
}
latest_file <- files[
  which.max(
    as.Date(sub(pattern, "\\1", basename(files)))
  )
]
message("Loading ", basename(latest_file))
manifests <- read.csv(latest_file, stringsAsFactors = FALSE)

manifests_ni <- manifests %>% filter(partner == "AFBN")
manifests_uk <- manifests %>% filter(partner != "AFBN")

# Get habitat types at 500m radius (I think what's 1km around a trap is the most reasonable to include atm - larger radius often catches other nearby traps and the trap identity is a strong predictor in all models so the local microhabitats are likely important)
# Load the map
land_dataUK <- terra::rast("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/habitat_complexity/ceh_maps/gblcm2024_10m.tif")
land_dataNI <- terra::rast("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/habitat_complexity/ceh_maps/nilcm2024_10m.tif")
message("Maps loaded")
# Select layer that holds the LCM-2020 classes & convert land cover to categorical factor
land_lcmUK <- land_dataUK[[1]]
land_lcmNI <- land_dataNI[[1]]
message("Map subset created")
# all.equal(land_dataUK[], land_dataNI[])
# cat("\nCombining datasets")
# land_data <- mosaic(land_dataUK, land_dataNI, fun = "first")
# cat("\nUK and NI succesfuly combined")
# Build the RAT (Raster Attribute Table)
lut <- data.frame(
  value    = 1:21,  
  category = c(
    "Broadleaf woodland", "Coniferous woodland",
    "Arable & horticulture", "Improved grassland",
    "Natural grassland", "Calcareous grassland",
    "Acid grassland", "Fen, marsh & swamp",
    "Heather", "Heather grassland", "Bog",
    "Inland rock", "Saltwater", "Freshwater",
    "Supralittoral rock", "Supralittoral sediment",
    "Littoral rock", "Littoral sediment", "Saltmarsh",
    "Urban", "Suburban")
)

# Process NI first

# Attach the table — an alternative is to use as.factor(land_data) but super computationally expensive on the 10m rasters
levels(land_lcmNI) <- lut  
# Convert occurrence points to SpatVector
occ_ni <- manifests_ni 
occ_points_ni <- vect(occ_ni, geom = c("sts_longitude", "sts_latitude"), crs = "EPSG:4326")
# Reproject points to match raster CRS (if needed)
occ_points_ni <- terra::project(occ_points_ni, crs(land_lcmNI))

# Process rest of UK second

# Attach the table — an alternative is to use as.factor(land_data) but super computationally expensive on the 10m rasters
levels(land_lcmUK) <- lut  
# Convert occurrence points to SpatVector
occ_uk <- manifests_uk
occ_points_uk <- vect(occ_uk, geom = c("sts_longitude", "sts_latitude"), crs = "EPSG:4326")
# Reproject points to match raster CRS (if needed)
occ_points_uk <- terra::project(occ_points_uk, crs(land_lcmUK))

buffer_width <- 500

cat("\nProcessing NI traps")

occ_points <- occ_points_ni
land_data <- land_dataNI
cat("\nCreating habitat types per trap")
# Establish habitat labels
hab_labels <- c(
  "Broadleaf woodland", "Coniferous woodland",
  "Arable & horticulture", "Improved grassland",
  "Natural grassland", "Calcareous grassland",
  "Acid grassland", "Fen, marsh & swamp",
  "Heather", "Heather grassland", "Bog",
  "Inland rock", "Saltwater", "Freshwater",
  "Supralittoral rock", "Supralittoral sediment",
  "Littoral rock", "Littoral sediment", "Saltmarsh",
  "Urban", "Suburban"
)
# Create buffer polygons around each trap
buf_poly <- terra::buffer(occ_points, width = buffer_width)
# Extract raster values inside each buffer
cat("Extract raster values from buffer\n")
landcover_within_buffers <- terra::extract(
  land_data, buf_poly,
  cells = TRUE, exact = TRUE
)
# Add trap metadata & pixel count
landcover_within_buffers$trap_name <- occ_points$trap_name[landcover_within_buffers$ID]
no_pixels <- table(landcover_within_buffers$trap_name) |>
  as.data.frame() |>
  setNames(c("trap_name", "no_pixels"))
landcover_within_buffers <- landcover_within_buffers |>
  dplyr::left_join(no_pixels, by = "trap_name")
# Translate land-cover codes to habitat labels
cat("Translate land-cover codes to habitat labels\n")
landcover_within_buffers <- landcover_within_buffers |>
  dplyr::mutate(
    habitat_type = factor(
      nilcm2024_10m_1,
      levels = 1:21,
      labels = hab_labels
    )
  )
landcover_within_buffers$habitat_type[is.na(landcover_within_buffers$habitat_type)] <- "Saltwater"
landcover_within_buffers_ni <- landcover_within_buffers
cat("\nFinished")

today_stamp <- format(Sys.Date(), "%Y-%m-%d")
file_out    <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/CEH2024_NI_trap_500m_%s.csv",
  today_stamp
)
write.csv(landcover_within_buffers_ni, file_out, row.names = FALSE)
cat("\nOutput created: NI")


cat("\nProcessing UK traps")
occ_points <- occ_points_uk
land_data <- land_dataUK
cat("\nCreating habitat types per trap")
# Establish habitat labels
hab_labels <- c(
  "Broadleaf woodland", "Coniferous woodland",
  "Arable & horticulture", "Improved grassland",
  "Natural grassland", "Calcareous grassland",
  "Acid grassland", "Fen, marsh & swamp",
  "Heather", "Heather grassland", "Bog",
  "Inland rock", "Saltwater", "Freshwater",
  "Supralittoral rock", "Supralittoral sediment",
  "Littoral rock", "Littoral sediment", "Saltmarsh",
  "Urban", "Suburban"
)
# Create buffer polygons around each trap
buf_poly <- terra::buffer(occ_points, width = buffer_width)
# Extract raster values inside each buffer
cat("Extract raster values from buffer\n")
landcover_within_buffers <- terra::extract(
  land_data, buf_poly,
  cells = TRUE, exact = TRUE
)
# Add trap metadata & pixel count
landcover_within_buffers$trap_name <- occ_points$trap_name[landcover_within_buffers$ID]
no_pixels <- table(landcover_within_buffers$trap_name) |>
  as.data.frame() |>
  setNames(c("trap_name", "no_pixels"))
landcover_within_buffers <- landcover_within_buffers |>
  dplyr::left_join(no_pixels, by = "trap_name")
# Translate land-cover codes to habitat labels
cat("Translate land-cover codes to habitat labels\n")
landcover_within_buffers <- landcover_within_buffers |>
  dplyr::mutate(
    habitat_type = factor(
      uklcm2024_10m_1,
      levels = 1:21,
      labels = hab_labels
    )
  )
landcover_within_buffers$habitat_type[is.na(landcover_within_buffers$habitat_type)] <- "Saltwater"
landcover_within_buffers_uk <- landcover_within_buffers
cat("\nFinished")

today_stamp <- format(Sys.Date(), "%Y-%m-%d")
file_out    <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/CEH2024_UK_trap_500m_%s.csv",
  today_stamp
)
write.csv(landcover_within_buffers_uk, file_out, row.names = FALSE)
cat("\nOutput created: UK")

# Combine
landcover_within_buffers <- rbind(landcover_within_buffers_uk, landcover_within_buffers_ni)

# Save the file
today_stamp <- format(Sys.Date(), "%Y-%m-%d")
file_out    <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/CEH2024_trap_500m_%s.csv",
  today_stamp
)
write.csv(landcover_within_buffers, file_out, row.names = FALSE)
cat("\nOutput created")

# Get the percentages of habitat types per 500m buffer for each trap 

# "Broadleaf woodland" == "Broadleaf woodland"
# "Coniferous woodland" == "Coniferous woodland"
# "Arable & horticulture" == "Arable"
# "Improved grassland" == "Improved grassland"
# "Natural grassland", "Calcareous grassland", "Acid grassland", "Fen, marsh & swamp" == "Semi-natural grassland
# "Heather", "Heather grassland", "Bog", "Inland rock" == "Mountain, heath and bog"
# "Saltwater" == "Saltwater" 
# "Freshwater" == "Freshwater"
# "Supralittoral rock", "Supralittoral sediment", "Littoral rock", "Littoral sediment", "Saltmarsh" == "Coastal"
# "Urban", "Suburban" == "Built-up areas and gardens" 
cat("\Calculating habitat type ratios per trap")
urban_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    Urban = sum(habitat_type %in% c("Urban", "Suburban"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Urban", "Suburban")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = Urban + Other,
    !!paste0("urban_ratio_500") :=
      (Urban) / total_pixels
  ) %>% 
  dplyr::select(trap_name, urban_ratio_500) %>% unique()

coastal_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    Coastal = sum(habitat_type %in% c("Supralittoral rock", "Supralittoral sediment", "Littoral rock", "Littoral sediment", "Saltmarsh"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Supralittoral rock", "Supralittoral sediment", "Littoral rock", "Littoral sediment", "Saltmarsh")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = Coastal + Other,
    !!paste0("coastal_ratio_500") :=
      (Coastal) / total_pixels
  ) %>% 
  dplyr::select(trap_name, coastal_ratio_500) %>% unique()

freshwater_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    Freshwater = sum(habitat_type %in% c("Freshwater"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Freshwater")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = Freshwater + Other,
    !!paste0("freshwater_ratio_500") :=
      (Freshwater) / total_pixels
  ) %>% 
  dplyr::select(trap_name, freshwater_ratio_500) %>% unique()

saltwater_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    Saltwater = sum(habitat_type %in% c("Saltwater"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Saltwater")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = Saltwater + Other,
    !!paste0("saltwater_ratio_500") :=
      (Saltwater) / total_pixels
  ) %>% 
  dplyr::select(trap_name, saltwater_ratio_500) %>% unique()

mountain_heath_bog_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    Bog = sum(habitat_type %in% c("Heather", "Heather grassland", "Bog", "Inland rock"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Heather", "Heather grassland", "Bog", "Inland rock")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = Bog + Other,
    !!paste0("mountain_heath_bog_ratio_500") :=
      (Bog) / total_pixels
  ) %>% 
  dplyr::select(trap_name, mountain_heath_bog_ratio_500) %>% unique()

semi_natural_grassland_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    NatGrass = sum(habitat_type %in% c("Natural grassland", "Calcareous grassland", "Acid grassland", "Fen, marsh & swamp"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Natural grassland", "Calcareous grassland", "Acid grassland", "Fen, marsh & swamp")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = NatGrass + Other,
    !!paste0("semi_natural_grassland_ratio_500") :=
      (NatGrass) / total_pixels
  ) %>% 
  dplyr::select(trap_name, semi_natural_grassland_ratio_500) %>% unique()

improved_grassland_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    ImpGrass = sum(habitat_type %in% c("Improved grassland"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Improved grassland")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = ImpGrass + Other,
    !!paste0("improved_grassland_ratio_500") :=
      (ImpGrass) / total_pixels
  ) %>% 
  dplyr::select(trap_name, improved_grassland_ratio_500) %>% unique()

arable_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    Arable = sum(habitat_type %in% c("Arable & horticulture"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Arable & horticulture")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = Arable + Other,
    !!paste0("agriculture_ratio_500") :=
      (Arable) / total_pixels
  ) %>% 
  dplyr::select(trap_name, agriculture_ratio_500) %>% unique()

broadleaf_woodland_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    Woodland = sum(habitat_type %in% c("Broadleaf woodland"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Broadleaf woodland")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = Woodland + Other,
    !!paste0("broadleaf_woodland_ratio_500") :=
      (Woodland) / total_pixels
  ) %>% 
  dplyr::select(trap_name, broadleaf_woodland_ratio_500) %>% unique()

coniferous_woodland_ratio <- landcover_within_buffers %>% 
  group_by(trap_name) %>% 
  summarise(
    Woodland = sum(habitat_type %in% c("Coniferous woodland"), na.rm = TRUE),
    Other    = sum(!(habitat_type %in% c("Coniferous woodland")), na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  mutate(
    total_pixels = Woodland + Other,
    !!paste0("coniferous_woodland_ratio_500") :=
      (Woodland) / total_pixels
  ) %>% 
  dplyr::select(trap_name, coniferous_woodland_ratio_500) %>% unique()

habitat_ratios <- Reduce(function(x, y) merge(x, y, by = "trap_name", all = TRUE),
                         list(urban_ratio,
                              coastal_ratio,
                              freshwater_ratio,
                              saltwater_ratio,
                              mountain_heath_bog_ratio,
                              semi_natural_grassland_ratio,
                              improved_grassland_ratio,
                              arable_ratio,
                              broadleaf_woodland_ratio,
                              coniferous_woodland_ratio))


# manifests <- merge(manifests, habitat_ratios, by = "trap_name")

# Save the file 
today_stamp <- format(Sys.Date(), "%Y-%m-%d") 
file_out    <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/habitat_ratios_%s.csv",
  today_stamp
)
write.csv(habitat_ratios, file_out, row.names = FALSE)
cat("\nOutput created")


