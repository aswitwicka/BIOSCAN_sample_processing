# This script takes the output of 00_manifest_fetch.sh
# and subsets the data for further processing 
# removes: control samples, all non-BIOSCAN partners
# retaines: malaise trap only, 2021 onwards only
# additional columns: trap, partner, region, plate, day, month, year, 24h sampling selection, if any samples from a given plate / catch lot have been sequenced, 
# The script also fetches all BOLD data and fills any BIN or sequence gaps that may not have synced

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

# Load the most recent data downloaded from sts [output of 00_manifest_fetch.sh]
dir_path <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output"
files <- list.files(dir_path, pattern = "^sts_manifests_[0-9]{8}\\.tsv$", full.names = TRUE)
dates <- as.Date(sub(".*_(\\d{8})\\.tsv", "\\1", files), format = "%Y%m%d")
latest_file <- files[which.max(dates)]
manifests <- read.table(latest_file,
                        sep = "\t",
                        header = TRUE,
                        fill = TRUE,
                        quote = ""
)
cat("Loaded file:", latest_file, "\n")
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

# Get mBRAVE data to assess if all sequenced samples have a QC category assigned 
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
cat("Do all the sequenced (present in mBRAVE) samples have a QC category assigned (anything that has status None)")
manifests %>% filter(sts_specimen.id %in% all_seq_samples) %>% pull(bioscan_qc_sanger_qc_result) %>% table()
cat("Do any samples that are not in mBRAVE have a QC category assigned for some reason?")
manifests %>% filter(!(sts_specimen.id %in% all_seq_samples)) %>% pull(bioscan_qc_sanger_qc_result) %>% table()

# Add sequenced column
manifests <- manifests %>%
  mutate(
    sequenced = case_when(
      bioscan_qc_sanger_qc_result %in% c("PASS", "ON_HOLD", "FAILED") ~ TRUE,
      bioscan_qc_sanger_qc_result == "None" ~ FALSE
    )
  )

# Remove empty negative controls and partial plate samples 
# Get sequenced plates (samples in mBRAVE)
seq_plates <- manifests %>% filter(sts_specimen.id %in% all_seq_samples) %>% pull(plate) %>% unique()
cat(paste("\nNumber of sequenced plates:", length(seq_plates)))
empty_neg_cont <- manifests %>% filter(plate %in% seq_plates) %>% filter(sequenced == FALSE) %>% pull(sts_specimen.id)
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
# The haversine [great circle] formula determines the great-circle distance between two points on a sphere given their longitudes and latitudes
haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371 # earth's radius in kilometers
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

# Add region column
# Load UK regions
Sys.setenv(SHAPE_RESTORE_SHX = "YES")
uk <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name == "United Kingdom")
shapefile_path <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/maps/ne_10m_admin_1_states_provinces.shp"
uk_regions <- st_read(shapefile_path)
# uk_regions <- st_read(shapefile_path, options = "ENCODING=UTF-8", stringsAsFactors = FALSE)
# uk_regions <- rnaturalearth::ne_states(
#   country = "United Kingdom",
#   returnclass = "sf"
# )
# Combine into wider regions 
uk_regions <- uk_regions[uk_regions$admin == "United Kingdom", ]
regions <- uk_regions %>%
  mutate(region = case_when(
    # South-east
    name %in% c("Greater London", "Kent", "Essex", "Surrey", "Sussex", "Hampshire", 
                "Oxfordshire", "Buckinghamshire", "Isle of Wight", "Medway", "Thurrock", 
                "Southend-on-Sea", "Brighton and Hove", "East Sussex", "West Sussex", 
                "Portsmouth", "Southampton", "Bromley", "Croydon", "Kingston upon Thames", 
                "Sutton", "Hillingdon", "Hertfordshire", "Hounslow", "Bexley", "Lambeth", "Hackney",
                "Ealing", "Barnet", "Brent", "Southwark", "Barking and Dagenham", "Wandsworth", "Hammersmith and Fulham",
                "Waltham Forest", "Harrow", "Lewisham", "Westminster", "Tower Hamlets", "Islington", "Camden", 
                "Richmond upon Thames", "Greenwich", "Cambridgeshire", "Havering", "Kensington and Chelsea", "Enfield",
                "Merton", "Redbridge", "Haringey", "Newham", "West Berkshire",
                "Norfolk", "Suffolk", "Peterborough", "Bedford", "Central Bedfordshire", "Luton",
                "Royal Borough of Windsor and Maidenhead", "Slough", "Bracknell Forest",
                "Reading", "Wokingham") ~ "South-east",
    # South-west
    name %in% c("Cornwall", "Devon", "Dorset", "Somerset", "Bristol", "Wiltshire", 
                "Gloucestershire", "Bath and North East Somerset", "South Gloucestershire", 
                "Plymouth", "Torbay", "Poole", "Bournemouth", "Isles of Scilly", 
                "Swindon", "North Somerset") ~ "South-west",
    # Midlands
    name %in% c("West Midlands", "East Midlands", "Warwickshire", "Staffordshire", 
                "Derbyshire", "Leicestershire", "Nottinghamshire", "Northamptonshire", 
                "Herefordshire", "Shropshire", "Worcestershire", "Telford and Wrekin", 
                "Stoke-on-Trent", "Rutland", "Nottingham", "Leicester", "Derby", 
                "Solihull", "Coventry", "Sandwell", "Dudley", "Walsall", "Wolverhampton", 
                "Milton Keynes", "Birmingham",
                "North Lincolnshire", "North East Lincolnshire", "Lincolnshire") ~ "Midlands",
    # North England
    name %in% c("Yorkshire", "Lancashire", "Cumbria", "Northumberland", "Durham", 
                "Merseyside", "Tyne and Wear", "Cheshire East", "Cheshire West and Chester", 
                "Manchester", "Liverpool", "Sheffield", "Leeds", "Bradford", "Wakefield", 
                "Barnsley", "Doncaster", "Rotherham", "Blackburn with Darwen", 
                "Calderdale", "Kirklees", "Stockton-on-Tees", "Darlington", 
                "Middlesbrough", "Gateshead", "North Tyneside", "South Tyneside", 
                "Sunderland", "Hartlepool", "Redcar and Cleveland", "York", 
                "Salford", "Bolton", "Trafford", "Oldham", "Rochdale", "Tameside", 
                "Bury", "Wigan", "Halton", "Knowsley", "Sefton", "Blackpool", 
                "City", "Newcastle upon Tyne", "Warrington", "North Yorkshire", "East Riding of Yorkshire",
                "Stockport", "Kingston upon Hull") ~ "North England",
    # Wales
    name %in% c("Powys", "Gwynedd", "Cardiff", "Swansea", "Ceredigion", "Pembrokeshire", 
                "Conwy", "Flintshire", "Wrexham", "Denbighshire", "Anglesey", 
                "Neath Port Talbot", "Vale of Glamorgan", "Rhondda, Cynon, Taff", 
                "Monmouthshire", "Blaenau Gwent", "Caerphilly", "Bridgend", "Torfaen", 
                "Merthyr Tydfil", "Newport", "Carmarthenshire") ~ "Wales",
    # Lowland Scotland
    name %in% c("Glasgow", "Edinburgh", "Fife", "Dundee", "Stirling", "Aberdeenshire", 
                "Argyll and Bute", "South Lanarkshire", "North Lanarkshire", "Aberdeen",
                "Renfrewshire", "East Ayrshire", "Clackmannanshire", "East Renfrewshire", 
                "West Lothian", "Falkirk", "Scottish Borders", "Dumfries and Galloway", 
                "East Dunbartonshire", "Midlothian", "East Lothian", "Angus", "Perthshire and Kinross",
                "Inverclyde", "West Dunbartonshire", "North Ayshire", "South Ayrshire") ~ "Lowland Scotland",
    # Highland Scotland
    name %in% c("Highland", "Orkney Islands", "Shetland Islands", "Western Isles", 
                "Moray", "Eilean Siar", "Orkney") ~ "Highland Scotland",
    # Northern Ireland
    name %in% c("Derry", "Strabane", "Fermanagh", "Dungannon", "Armagh", "Newry and Mourne", 
                "Antrim", "Lisburn", "Ballymoney", "Ballymena", "Magherafelt", "Craigavon", 
                "Banbridge", "Coleraine", "Moyle", "Larne", "Carrickfergus", "Newtownabbey", 
                "North Down", "Ards", "Down", "Belfast", "Castlereagh", "Omagh", 
                "Mid Ulster", "Limavady") ~ "Northern Ireland",
    # East 
    name %in% c() ~ "East England",
    # Catch-all for unassigned
    TRUE ~ "Unassigned"
  )) 
# Make points
pts_sf <- st_as_sf(
  manifests,
  coords = c("sts_longitude", "sts_latitude"),
  crs = 4326,          
  remove = FALSE       
)
# Ensure regions is sf and CRS matches
stopifnot(inherits(regions, "sf"), "region" %in% names(regions))
if (is.na(st_crs(regions))) st_crs(regions) <- 4326
pts_sf <- st_transform(pts_sf, st_crs(regions))
# Optional - avoids topology errors
regions <- sf::st_make_valid(regions)
# Spatial join
pts_join <- st_join(pts_sf, regions[, c("region")], join = st_within)  # or st_intersects
# Return to data frame format 
manifests <- st_drop_geometry(pts_join)
# Fix main issues manually for now
manifests <- manifests %>% 
  mutate(
    region = case_when(
      grepl("HLNR1", trap_name) ~ "Highland Scotland",
      grepl("NMCC",  trap_name) ~ "Lowland Scotland",
      TRUE ~ region 
    )
  )
# Make sure
is.na(manifests$region) %>% table()
table(manifests$region)

# Additinal columns

# Remove catch lots that should not be here
manifests <- manifests %>% filter(!(sts_CATCH_LOT %in% c("6M", "5M", "NOT_APPLICABLE")))
# "NOT_APPLICABLE" is from NHS - catches over 24h now excluded

# Add updated catch lot column
manifests$catch_lot <- paste(manifests$trap_name, manifests$year, manifests$month, manifests$day, sep = "_")

cat("Catch lot numbers that appear across dates/partners [ERRORS]:\n")
manifests %>%
  group_by(sts_CATCH_LOT) %>%
  filter(n_distinct(catch_lot) > 1) %>% pull(sts_CATCH_LOT) %>% table

# Add coluns with sequencing success
# Catch lot level
manifests <- manifests %>%
  group_by(catch_lot) %>%
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

# Uniform catch duration across catch lots 

# Find catch lots where duration24 contains both "YES" and "NO" (inconsistently entered by a partner)
lots_with_both <- manifests %>%
  group_by(catch_lot) %>%
  filter(any(duration24 == "YES") & any(duration24 == "NO")) %>%
  ungroup()

cat("Catch lots with inconsistent sampling duration codes:\n")
lots_with_both %>% pull(catch_lot) %>% unique()
lots_with_both %>% pull(sts_CATCH_LOT) %>% unique()

# Correct duration24 within these groups (if half or more have 24h sampling then correct for all )
manifests <- manifests %>%
  group_by(catch_lot) %>%
  mutate(
    prop_yes = mean(duration24 == "YES"),
    duration24 = if_else(prop_yes > 0.5, "YES", duration24)
  ) %>%
  dplyr::select(-prop_yes) %>%
  ungroup()

# Display catch_lots that STILL have any “NO” remaining
remaining_no <- manifests %>%
  filter(catch_lot %in% lots_with_both$catch_lot) %>% 
  group_by(catch_lot) %>%
  filter(any(duration24 == "NO")) %>%
  ungroup()

cat("Catch lots with inconsistent sampling duration codes after correction (more than 50% of samples have not been samples for 24h):\n")
remaining_no %>% pull(catch_lot) %>% unique()
remaining_no %>% pull(sts_CATCH_LOT) %>% unique()

# remaining_no %>%
#   dplyr::select(catch_lot, sts_CATCH_LOT, trap_name, duration24, actual_duration) %>%
#   arrange(catch_lot) %>% View()

summary_table <- remaining_no %>%
  group_by(catch_lot, sts_CATCH_LOT, trap_name) %>%
  summarize(
    how_many_24 = sum(duration24 == "YES", na.rm = TRUE),
    how_many_other = sum(duration24 == "NO",  na.rm = TRUE),
    .groups = "drop"
  )

summary_table

# Update BOLD data
library(BOLDconnectR)
bold.apikey('32D02751-6F67-4211-BB84-4892748D5F95')
# Get the data
project_codes <- c("BSCAN")
BSCAN_data <- bold.fetch(get_by = "project_codes", identifiers = project_codes)
# Remove columns with only NA values
BSCAN_data_filt <- BSCAN_data %>%
  dplyr::select(where(~ !all(is.na(.))))
# Select only samples with sequences
BSCAN_data_sequenced <- BSCAN_data_filt %>% filter(!(is.na(nuc))) %>% unique()
# nrow(BSCAN_data_filt)
# nrow(BSCAN_data_sequenced)
# Add plate
BSCAN_data_sequenced$plate <- gsub(".{3}$", "", BSCAN_data_sequenced$sampleid)
BSCAN_data_sequenced$plate <- sub("_$", "", BSCAN_data_sequenced$plate)

BSCAN_data_sequenced <- BSCAN_data_sequenced %>% dplyr::select(sampleid, bin_created_date, bin_uri, class, collectors, family, genus, kingdom, nuc, 
                                                               nuc_basecount, order, phylum, sequence_run_site, sequence_upload_date, 
                                                               species, subfamily, tribe)
colnames(BSCAN_data_sequenced) <- c("sts_specimen.id", "bold_bin_created_date", "bold_bin_uri", "bold_class", "bold_collectors", "bold_family", 
                                    "bold_genus", "bold_kingdom", "bold_nuc",  "bold_nuc_basecount", "bold_order", "bold_phylum", 
                                    "bold_sequence_run_site", "bold_sequence_upload_date", "bold_species", "bold_subfamily", "bold_tribe")



# Define the columns to replace
cols_to_replace <- c("bold_bin_created_date", "bold_bin_uri", "bold_class", "bold_collectors", 
                     "bold_family", "bold_genus", "bold_kingdom", "bold_nuc", "bold_nuc_basecount", 
                     "bold_order", "bold_phylum", "bold_sequence_run_site", "bold_sequence_upload_date", 
                     "bold_species", "bold_subfamily", "bold_tribe")

# Subset sequnced samples that don't have BOLD records 
data_merge_subset <- manifests %>%
  filter(bold_nuc == "None", sts_specimen.id %in% BSCAN_data_sequenced$sts_specimen.id) 

BSCAN_data_sequenced <- BSCAN_data_sequenced %>% filter(sts_specimen.id %in% data_merge_subset$sts_specimen.id)

data_merge_subset_updated <- data_merge_subset %>%
  dplyr::select(-all_of(cols_to_replace)) %>% 
  left_join(BSCAN_data_sequenced %>% dplyr::select(sts_specimen.id, all_of(cols_to_replace)), by = "sts_specimen.id")

data_merge_rest <- manifests %>%
  filter(!(sts_specimen.id %in% data_merge_subset$sts_specimen.id))

data_merge_final <- rbind(data_merge_rest, data_merge_subset_updated)

# Check
nrow(manifests)
nrow(data_merge_final)

# Summary
table(data_merge_final$duration24)
table(data_merge_final$sts_COLLECTION_METHOD)
table(data_merge_final$sts_organism_part)
table(data_merge_final$bioscan_qc_sanger_qc_result)
table(data_merge_final$sequenced)
table(data_merge_final$sequenced_plate_level)
table(data_merge_final$sequenced_lot_level)
cat("No sts_CATCH_LOT:\n")
length(unique(data_merge_final$sts_CATCH_LOT))
cat("No corrected catch_lot:\n")
length(unique(data_merge_final$catch_lot))

table(data_merge_final$bold_class)
table(data_merge_final$bold_kingdom)

data_merge_final %>% filter(bold_bin_uri == "None") %>% nrow()
manifests %>% filter(bold_bin_uri == "None") %>% nrow()
data_merge_final %>% filter(bold_nuc == "None") %>% nrow()
manifests %>% filter(bold_nuc == "None") %>% nrow()

# Save the file 
today_stamp <- format(Sys.Date(), "%Y-%m-%d") 
file_out    <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/BIOSCAN_100k_samples_corrected%s.csv",
  today_stamp
)
write.csv(data_merge_final, file_out, row.names = FALSE)
cat("\nOutput created")

# This dataset containds all BIOSCAN samples (Malaise traps + 24h (±2h) sampling and all other > 24h samples)
# It includes all samples sampled, sequenced, QC-ed
# It also includes all catch lots regardless of sequencing success 
# and all samples regardless of taxonomy (arthropod or not)
# Therefore is must be processed before using in in any further analyses 
