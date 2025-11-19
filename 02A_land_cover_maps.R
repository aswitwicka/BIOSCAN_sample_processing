# This script takes the list of traps and their locations and creates the habitat type description within a given radius size (a csv file for each radius)
# NOTE: Norhtern Ireland must be processed seperately

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
  "phangorn", "shiny", "sf", "tibble", "forcats", "lubridate", "viridis", "maps"
)
# Bioconductor packages
bioconductor_pkgs <- c(
  "biomaRt", "Biostrings", "msa", "ape"
)
# Load CRAN packages
load_pkgs(cran_pkgs, bioconductor = FALSE)
# Load Bioconductor packages
load_pkgs(bioconductor_pkgs, bioconductor = TRUE)

# Working at the trap level here - subset
working_set <- read.csv("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/required_files/trap_to_partner.csv")

colnames(working_set) <- c("partner", "trap_name", "trap_no", "sts_latitude", "sts_longitude")

cat(paste("Number of traps:", nrow(working_set),
          "\nCheck if anything is repeated:", length(unique(working_set$trap_name))))

# Set date
today_stamp <- format(Sys.Date(), "%Y-%m-%d") 

### Assign habitat types ###

# Get env info 
land_dataGB <- terra::rast("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/maps/gblcm2024_10m.tif")
land_dataNI <- terra::rast("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/maps/nilcm2024_10m.tif")

out_dir <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output"

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

# Specify size 
buffers <- c(25, 50, 100, 500, 1000)

# Process NI first
land_data <- land_dataNI

# Convert land cover to categorical factor
# Select layer that holds the LCM-2020 classes
land_lcm <- land_data[[1]]

# Attach the table — an alternative is to use as.factor(land_data) but super computationally expensive on the 10m rasters and it was crushing... But good for the 1km data 
levels(land_lcm) <- lut

# Convert occurrence points to SpatVector
# Select only NI traps 
occ_df <- working_set %>% filter(partner == "AFBN")
occ_points <- vect(occ_df, geom = c("sts_longitude", "sts_latitude"), crs = "EPSG:4326")

# Reproject points to match raster CRS (if needed)
occ_points <- project(occ_points, crs(land_lcm))

# Extract raster values at each point
land_cover_values <- terra::extract(land_lcm, occ_points)

# Merge with original data
occ_df$land_cover <- land_cover_values[,2]  # column 2 contains the raster values

# Create buffer files

cat("\nStarting the loop")
for (bw in buffers) {
  cat("Make a buffer polygon around every trap")
  buf_poly <- terra::buffer(occ_points, width = bw)
  cat(bw)
  # Extract raster values inside each buffer
  landcover_within_buffers <- terra::extract(
    land_data, buf_poly,
    cells = TRUE, exact = TRUE
  )
  cat("Add trap metadata + pixel count")
  landcover_within_buffers$trap_name    <- occ_points$trap_name[landcover_within_buffers$ID]
  no_pixels <- table(landcover_within_buffers$trap_name) |>
    as.data.frame() |>
    setNames(c("trap_name", "no_pixels"))
  landcover_within_buffers <- landcover_within_buffers |>
    left_join(no_pixels, by = "trap_name")
  cat("Translate land-cover codes to habitat labels")
  landcover_within_buffers <- landcover_within_buffers |>
    mutate(
      habitat_type = factor(
        # gblcm2024_10m_1,
        nilcm2024_10m_1,
        levels = 1:21,
        labels = hab_labels
      )
    )
  
  cat("Write output")
  table(landcover_within_buffers$trap_name)
  file_out <- file.path(out_dir,
                        # sprintf("%d_buffer_2024_%s.csv", bw, today_stamp))
                        sprintf("%d_buffer_NIL_2024_%s.csv", bw, today_stamp))
  write.csv(landcover_within_buffers, file_out, row.names = FALSE)
  message("Saved ", basename(file_out))
}

# Number of pixels per trap 
#table(landcover_within_buffers$trap_name)

# Process GB next
land_data <- land_dataGB

# Convert land cover to categorical factor
# Select layer that holds the LCM-2020 classes
land_lcm <- land_data[[1]]

# Attach the table — an alternative is to use as.factor(land_data) but super computationally expensive on the 10m rasters and it was crushing... But good for the 1km data 
levels(land_lcm) <- lut

# Convert occurrence points to SpatVector
# Select only NI traps 
occ_df <- working_set %>% filter(partner != "AFBN")
occ_points <- vect(occ_df, geom = c("sts_longitude", "sts_latitude"), crs = "EPSG:4326")

# Reproject points to match raster CRS (if needed)
occ_points <- project(occ_points, crs(land_lcm))

# Extract raster values at each point
land_cover_values <- terra::extract(land_lcm, occ_points)

# Merge with original data
occ_df$land_cover <- land_cover_values[,2]  # column 2 contains the raster values

# Create buffer files

cat("\nStarting the loop")
for (bw in buffers) {
  cat("Make a buffer polygon around every trap")
  buf_poly <- terra::buffer(occ_points, width = bw)
  cat(bw)
  # Extract raster values inside each buffer
  landcover_within_buffers <- terra::extract(
    land_data, buf_poly,
    cells = TRUE, exact = TRUE
  )
  cat("Add trap metadata + pixel count")
  landcover_within_buffers$trap_name    <- occ_points$trap_name[landcover_within_buffers$ID]
  no_pixels <- table(landcover_within_buffers$trap_name) |>
    as.data.frame() |>
    setNames(c("trap_name", "no_pixels"))
  landcover_within_buffers <- landcover_within_buffers |>
    left_join(no_pixels, by = "trap_name")
  cat("Translate land-cover codes to habitat labels")
  landcover_within_buffers <- landcover_within_buffers |>
    mutate(
      habitat_type = factor(
        gblcm2024_10m_1,
        # nilcm2024_10m_1,
        levels = 1:21,
        labels = hab_labels
      )
    )
  
  cat("Write output")
  table(landcover_within_buffers$trap_name)
  file_out <- file.path(out_dir,
                        sprintf("%d_buffer_GBL_2024_%s.csv", bw, today_stamp))
                        # sprintf("%d_buffer_NIL_2024_%s.csv", bw, today_stamp))
  write.csv(landcover_within_buffers, file_out, row.names = FALSE)
  message("Saved ", basename(file_out))
}

