# This script uses CEH maps to visualise locations of the traps [map cut-outs and pie charts]

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
  "ggplot2", "tidyr", "stringr", "terra", "dismo", "rlang", "viridis",
  "parallel", "bigmemory", "raster", "ncdf4", "seqinr", "vegan", "reshape2", "remotes",
  "phangorn", "shiny", "sf", "textshape", "tibble", "forcats", "lubridate", "maps"
)
# Bioconductor packages
bioconductor_pkgs <- c(
  "biomaRt", "Biostrings", "msa", "ape"
)
# Load CRAN packages
load_pkgs(cran_pkgs, bioconductor = FALSE)
# Load Bioconductor packages
load_pkgs(bioconductor_pkgs, bioconductor = TRUE)

# Load csv file (load the one you need in terms of the radius size)
map2023 <- read.csv("~/bioscan/habitat_complexity/output/intermediary_files/500_buffer_NIL_2024_2025-10-29.csv")
# map2023 <- read.csv("~/bioscan/habitat_complexity/output/intermediary_files/5000_buffer_2021_2025-09-05.csv")

land_data <- terra::rast("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/habitat_complexity/ceh_maps/gblcm2024_10m.tif")
# Convert land cover to categorical factor
# Select layer that holds the LCM-2020 classes
land_lcm <- land_data[[1]]
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
# Attach the table — an alternative is to use as.factor(land_data) but super computationally expensive on the 10m rasters
levels(land_lcm) <- lut

custom_palette <- c(
  #"Coastal" = "#05aae6",
  "Broadleaf woodland" = "#c75a06",
  "Coniferous woodland" = "#eb8f49",
  "Arable & horticulture" = "#dbd21d",
  "Improved grassland" = "#90d1ad",
  "Natural grassland" = "#135901",
  "Calcareous grassland" = "#06c75d",
  "Acid grassland" = "#89f26f",
  "Fen, marsh & swamp" = "#351440",
  "Heather" = "#5d1075",
  "Heather grassland" = "#d69ee8",
  "Bog" = "#e6d1ed",
  "Inland rock" = "#926d9e",
  "Saltwater" = "#05aae6",
  "Freshwater" = "#acd8e8",
  "Supralittoral rock" = "#0c5399",
  "Supralittoral sediment" = "#34699e",
  "Littoral rock" = "#0c2d99",
  "Littoral sediment" = "#3852a8",
  "Saltmarsh" = "#6562fc",
  "Urban" = "#212021",
  "Suburban" = "#b0aeb0"
)


# Function to plot habitat map for a given trap
plot_trap_habitat <- function(trap_id, map_df, palette) {
  # filter to one trap
  trap_data <- map_df %>%
    filter(trap_name == trap_id) %>%
    dplyr::select(cell, habitat_type)
  
  # convert cell ids to xy coordinates (terra::cellFromXY works on raster, but we can decode if map_df already has xy cols)
  # If you have x/y columns in map2023, use them. Otherwise derive from raster:
  coords <- terra::xyFromCell(land_lcm, trap_data$cell)
  trap_df <- cbind(as.data.frame(coords), habitat_type = trap_data$habitat_type)
  
  # plot
  ggplot(trap_df, aes(x = x, y = y, fill = habitat_type)) +
    geom_raster() +
    scale_fill_manual(values = palette) +
    coord_equal() +
    labs(title = paste("Trap:", trap_id),
         fill = "Habitat") +
    theme_classic()
}

# Example: plot for one trap
plot_trap_habitat("RAIM", map2023, custom_palette)

# If you want to save all traps to files:
unique_traps <- unique(map2023$trap_name)

for (t in unique_traps) {
  p <- plot_trap_habitat(t, map2023, custom_palette)
  ggsave(
    filename = sprintf("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/habitat_complexity/output/processing_plots/trap_maps/2024trap_maps/trap_%s_50m_habitat_map.png", t),
    plot = p,
    width = 6, height = 6, dpi = 300
  )
}

# Pie charts for each trap 
outdir <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/habitat_complexity/output/processing_plots/trap_maps/pie_charts_2024/"
# Loop through trap_name
for(trap in unique(map2023$trap_name)) {
  # Summarize % habitat composition for this trap
  trap_data <- map2023 %>%
    filter(trap_name == trap) %>%
    group_by(habitat_type) %>%
    summarise(total = n()) %>%
    mutate(perc = total / sum(total) * 100)
  # Skip if no valid data
  if(nrow(trap_data) == 0) next
  # Pie chart
  p <- ggplot(trap_data, aes(x = "", y = perc, fill = habitat_type)) +
    geom_bar(stat = "identity", width = 0.5, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = custom_palette) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  # Save plot
  outfile <- file.path(outdir, paste0(trap, "_2024ceh_piechart_500mRadius.pdf"))
  ggsave(outfile, p, width = 7, height = 7)
}




# Comparison 2023 and 2021 (this has been done already)
# Calculate indices for 2021:
base_buf_name <- "map2021_5000"

div_list <- map2023 %>% 
    group_by(trap_name, habitat_type) %>% 
    summarise(area = sum(fraction, na.rm = TRUE), .groups = "drop") %>% 
    ungroup() %>%
    group_by(trap_name) %>% 
    summarise(
      !!paste0("shannon_",  base_buf_name) :=
        { p <- area / sum(area); -sum(p * log(p)) },
      !!paste0("simpson_",  base_buf_name) :=
        { p <- area / sum(area); 1 - sum(p ^ 2) },
      .groups = "drop"
    )

# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/habitat_complexity/output/intermediary_files/05_2021COMPARISON_map5000m_diversity_indices_%s.csv",
  today_stamp
)
write.csv(div_list, file_out, row.names = FALSE)
