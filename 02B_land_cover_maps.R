# This script uses csv files from 02A_land_cover_maps.R
# It combines matching NI and GB files and creates a summarised .rds file to use in subsequent analyses 

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

# Load the most up-to-date input file
dir_in   <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output"

pattern_gb  <- "*buffer_GBL_2024_(\\d{4}-\\d{2}-\\d{2})\\.csv$"   # YYYY-MM-DD
files_gb <- list.files(dir_in, pattern = pattern_gb, full.names = TRUE)
pattern_n  <- "*buffer_NIL_2024_(\\d{4}-\\d{2}-\\d{2})\\.csv$"   # YYYY-MM-DD
files_n <- list.files(dir_in, pattern = pattern_n, full.names = TRUE)
files <- c(files_gb, files_n)

if (length(files) == 0) {
  stop("No *_buffer_*.csv files found in ", dir_in)
}
message("Loading ", length(files), " files …")
working_sets <- lapply(files, read.csv, stringsAsFactors = FALSE)

extract_meta <- function(fname) {
  bn <- basename(fname)
  buffer_val <- suppressWarnings(as.integer(sub("_buffer.*", "", bn)))
  date_str <- sub(".*_buffer_(GBL|NIL)_2024_(\\d{4}-\\d{2}-\\d{2})\\.csv$", "\\2", bn)
  if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", date_str)) {
    date_val <- NA
  } else {
    date_val <- as.Date(date_str)
  }
  c(buffer = buffer_val, date_stamp = as.numeric(date_val))
}

meta <- t(vapply(files, extract_meta, numeric(2)))   
meta <- data.frame(
  buffer     = meta[, "buffer"],
  date_stamp = as.Date(meta[, "date_stamp"], origin = "1970-01-01")
)

cat("Radius files used:\n")
print(meta)

names(working_sets) <- rownames(meta)

names(working_sets) <- names(working_sets) |>
  basename() |> # drop path
  sub("_\\d{4}-\\d{2}-\\d{2}\\.csv$", "", x = _)  

cat("Saving the rds file")
# Save as R object 
saveRDS(working_sets, 
        file = "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02B_working_sets_radius.rds")
cat("Finished and saved")
