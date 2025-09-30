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

# Load the data from sts 
manifests <- read.table("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/sts_manifests_22092025.tsv", 
                        sep = "\t",
                        header = TRUE,
                        fill = TRUE,       # still useful for short rows
                        quote = ""
)

# General check
cat(paste("Are any samples duplicated:", !(length(unique(manifests$sts_specimen.id)) == nrow(manifests))))

### PROCESSING AND CLEANING ###

# Create plate column
manifests$plate <- gsub(".{3}$", "", manifests$sts_specimen.id)
manifests$plate <- sub("_$", "", manifests$plate)
# Replace MOZZ samples 
mozz_to_partner <- read.csv("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/required_files/mozz_to_partner.csv")
manifests <- manifests %>%
  mutate(plate = ifelse(plate %in% mozz_to_partner$Label,
                        mozz_to_partner$partner[match(plate, mozz_to_partner$Label)],
                        plate))
cat(paste("Number of samples where plate name standarisation did not work:", 
          manifests %>% filter(str_detect("MOZ", plate)) %>% nrow()))

# Create partner column 
manifests <- manifests %>% 
  mutate(partner = sub("_.*", "", plate))

# Remove non bioscan manifests from the list 
partner_patterns <- c("BGEP", "SNST", "JARO", "BGKU", "BSN", "BGEG", "OXHP", "POMS", "AYDI", "BGPT")
manifests <- manifests %>% filter(!(partner %in% partner_patterns))
manifests <- manifests %>% filter(!(str_detect(partner, "BGE")))
manifests <- manifests %>% filter(!(str_detect(partner, "TOL-"))) 

# Remove controls 
manifests <- manifests[!(grepl("_G12|_H12", manifests$sts_specimen.id)),]

# Correct mistakes (not sure id these has been corrected systematically)
# Correct the JHGS partner code 
manifests <- manifests %>%
  mutate(partner = gsub("JHGS-00[1-3]", "JHGS", partner))
# Correct coordinates in NOTT: 
manifests <- manifests %>% 
  mutate(
    sts_longitude = if_else(partner == "NOTT", "-0.5479", sts_longitude),
    sts_latitude  = if_else(partner == "NOTT", "52.8581", sts_latitude)
  )
# Correct coordinates in NNWT: 
manifests <- manifests %>% 
  mutate(
    sts_longitude = if_else(partner == "NWWT", "-4.3328", sts_longitude),
    sts_latitude  = if_else(partner == "NWWT", "52.9998", sts_latitude)
  )
# Correct coordinates in OUMK: 
manifests <- manifests %>%
  mutate(
    sts_longitude = if_else(
      partner == "OUMK" & !str_starts(sts_longitude, "-"),
      paste0("-", sts_longitude),
      sts_longitude
    )
  )

# Get columns: $month $day $year 
manifests$year <- substr(manifests$sts_col_date, 1, 4)
manifests$month <- substr(manifests$sts_col_date, 6, 7)
manifests$day <- substr(manifests$sts_col_date, 9, 10)
manifests$day[manifests$day == ""] <- "None"
manifests$month[manifests$month == ""] <- "None"

# Get mBRAVE data, get FAILED samples, assess plate and catch lot success rates, get info about batches and failed plates 
# Find all files ending with "_sample_stats.txt" in mBRAVE directories 
files <- list.files(path = "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/bioscan_qc/mbrave_batch_data", 
                    pattern = "_sample_stats.txt$", 
                    recursive = TRUE, 
                    full.names = TRUE)
# Combine all files into one long data frame
all_seq_samples_table <- files %>%
  lapply(function(f) {
    df <- tryCatch({
      read.delim(f, header = TRUE, stringsAsFactors = FALSE, fill = TRUE) %>%
        mutate(batch = str_extract(basename(f), "(?<=umi\\.)[^_]+.*(?=_sample_stats\\.txt)"))
    }, error = function(e) {
      message("Error reading file: ", f)
      return(NULL)
    })
    return(df)
  }) %>%
  bind_rows()

all_seq_samples <- unique(all_seq_samples_table$Label) 
all_seq_samples <- str_replace(all_seq_samples, "-SDC", "")
# Are all these samples in the manifest?
manifests %>% filter(!(sts_specimen.id %in% all_seq_samples)) %>% filter(bold_nuc != "None") %>% pull(plate) %>% table()
# Update QC output
manifests <- manifests %>%
  mutate(
    qc_res_updated = case_when(
      bioscan_qc_sanger_qc_result %in% c("PASS", "ON_HOLD") ~ bioscan_qc_sanger_qc_result,
      bioscan_qc_sanger_qc_result == "None" & sts_specimen.id %in% all_seq_samples ~ "FAILED",
      bioscan_qc_sanger_qc_result == "None" ~ "Not in mBRAVE"
    )
  )

# Add comuns with sequencing success
# Catch lot level
manifests <- manifests %>%
  group_by(sts_CATCH_LOT) %>%
  mutate(sequenced_lot_level = if_else(
    any(sts_specimen.id %in% all_seq_samples, na.rm = TRUE),
    "TRUE", "FALSE"
  )) %>%
  ungroup()
# Plate level
manifests <- manifests %>%
  group_by(plate) %>%
  mutate(sequenced_plate_level = if_else(
    any(sts_specimen.id %in% all_seq_samples, na.rm = TRUE),
    "TRUE", "FALSE"
  )) %>%
  ungroup()

# Remove empty negative controls and partial plate samples 
# Get sequenced plates (samples in mBRAVE)
seq_plates <- manifests %>% filter(sts_specimen.id %in% all_seq_samples) %>% pull(plate) %>% unique()
cat(paste("\nNumber of sequenced plates:", length(seq_plates)))
empty_neg_cont <- manifests %>% filter(plate %in% seq_plates) %>% filter(qc_res_updated == "Not in mBRAVE") %>% pull(sts_specimen.id)
manifests <- manifests %>% filter(!(sts_specimen.id %in% empty_neg_cont))

# Upload sampling duration codes 
duration_codes <- read.csv("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/required_files/collection_time_codes.csv")
duration_codes <- duration_codes %>% dplyr::select(sts_DURATION_OF_COLLECTION, actual_duration, duration24) 
# Is there anything new?
missing_values <- setdiff(unique(manifests$sts_DURATION_OF_COLLECTION), unique(duration_codes$sts_DURATION_OF_COLLECTION))
match_table <- table(unique(manifests$sts_DURATION_OF_COLLECTION) %in% unique(duration_codes$sts_DURATION_OF_COLLECTION))
# Display if any values are missing
if ("FALSE" %in% names(match_table) && match_table["FALSE"] > 0) {
  cat("Missing:\n", paste(missing_values, collapse = "\n"))
}
# Change 
manifests <- manifests %>%                 
  left_join(                        
    duration_codes,
    by = "sts_DURATION_OF_COLLECTION"
  )

# Subset true BIOSCAN samples
# Subset 2021-2025 only
manifests <- manifests %>% filter(year %in% c("2021", "2022", "2023", "2024", "2025", "2026"))
# Subset Malaise traps
manifests <- manifests %>% filter(sts_COLLECTION_METHOD == "MALAISE_TRAP")
# Subset 24h only
duration_histogram <- ggplot(manifests %>% filter(actual_duration < 100), aes(x = actual_duration)) +
  geom_histogram(binwidth = 1, boundary = 0,
                 fill = "#d69ee8") +
  geom_vline(xintercept = c(21, 26), linetype = "dashed",
             linewidth = 0.2, colour = "#351440") +
  scale_x_continuous(
    name   = "Sampling duration (hours)",
    breaks = seq(0, 100, by = 6),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(y = "Number of samples") +
  theme_classic(base_size = 10)
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/duration_histogram.pdf",
  plot   = duration_histogram,
  device = "pdf",
  width  = 10, height = 6, units = "cm"
)
manifests <- manifests %>% filter(duration24 == "YES")
cat("\nManifest subsets created")

# Add trap info
trap_per_partner <- read.csv("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/required_files/trap_to_partner.csv")
# Establish trap numbers
# Limit latitude and longitude to 4 decimal places
manifests$sts_latitude <- as.numeric(manifests$sts_latitude)
manifests$sts_longitude <- as.numeric(manifests$sts_longitude)
#
coordinates_vs_traps <- manifests
coordinates_vs_traps$sts_latitude <- round(coordinates_vs_traps$sts_latitude, 4)
coordinates_vs_traps$sts_longitude <- round(coordinates_vs_traps$sts_longitude, 4)
# Are there any new partners that are not in the trap file?
new_partners <- unique(coordinates_vs_traps$partner[!(coordinates_vs_traps$partner %in% trap_per_partner$partner)])
cat(paste("Partners not included in the trap list:", paste(new_partners, collapse = ", ")))
# Function to calculate
# The haversine [great circle] formula determines the great-circle distance between two points on a sphere given their longitudes and latitudes.
haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371 # Earth's radius in kilometers
  delta_lat <- (lat2 - lat1) * pi / 180
  delta_lon <- (lon2 - lon1) * pi / 180
  a <- sin(delta_lat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(delta_lon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R * c
}
# Match traps 
assign_trap_info <- function(coordinates_vs_traps, trap_per_partner, radius = 5) {
  coordinates_vs_traps$trap_no <- NA
  coordinates_vs_traps$trap_name <- NA
  # Iterate through rows
  for (i in seq_len(nrow(coordinates_vs_traps))) {
    # Get partner and coordinates info
    partner <- coordinates_vs_traps$partner[i]
    lat <- coordinates_vs_traps$sts_latitude[i]
    lon <- coordinates_vs_traps$sts_longitude[i]
    
    partner_traps <- subset(trap_per_partner, partner == partner)
    
    if (nrow(partner_traps) == 0) next # Skip if no traps found
    
    if (all(partner_traps$trap_no == 1)) {
      coordinates_vs_traps$trap_no[i] <- partner_traps$trap_no[1]
      coordinates_vs_traps$trap_name[i] <- partner_traps$trap_name[1]
    } else {
      # Calculate distance
      distances <- sapply(1:nrow(partner_traps), function(j) {
        haversine(lat, lon,
                  partner_traps$DECIMAL_LATITUDE[j],
                  partner_traps$DECIMAL_LONGITUDE[j])
      })
      
      # Find the closest trap within the radius
      closest_idx <- which.min(distances)
      if (distances[closest_idx] <= radius) {
        coordinates_vs_traps$trap_no[i] <- partner_traps$trap_no[closest_idx]
        coordinates_vs_traps$trap_name[i] <- partner_traps$trap_name[closest_idx]
      } else {
        # If no trap found place NA
        coordinates_vs_traps$trap_no[i] <- NA
        coordinates_vs_traps$trap_name[i] <- NA
      }
    }
  }
  
  return(coordinates_vs_traps)
}
# Apply
cat("\nCreating trap lists")
coordinates_vs_traps <- assign_trap_info(coordinates_vs_traps, trap_per_partner, radius = 5)
# coordinates_vs_traps <- coordinates_vs_traps %>%
#   mutate(trap_name = ifelse(is.na(trap_name), partner, trap_name))
cat(paste("All detected traps:", paste(unique(coordinates_vs_traps$trap_name), collapse = ", "),
          "\n\nNo. traps:", length(unique(coordinates_vs_traps$trap_name))))
cat(paste("Traps not detected:",
          paste(trap_per_partner$trap_name[
            !(trap_per_partner$trap_name %in%
                unique(coordinates_vs_traps$trap_name))
          ], collapse = ", "))
)
cat("Problematic samples that have wrong coordinates or the traps has been significantly relocated:")
coordinates_vs_traps[is.na(coordinates_vs_traps$trap_name),] %>% nrow()
coordinates_vs_traps[is.na(coordinates_vs_traps$trap_name),] %>% dplyr::select(partner, sts_latitude, sts_longitude, trap_name) %>% unique()
manifests <- coordinates_vs_traps
cat("\nTrap names assigned")

# Save the file 
today_stamp <- format(Sys.Date(), "%Y-%m-%d") 
file_out    <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/BIOSCAN_100k_samples_corrected%s.csv",
  today_stamp
)
write.csv(manifests, file_out, row.names = FALSE)
cat("\nOutput created")





