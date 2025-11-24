# This script creates a "calendar" plot with all catch plots that exist in the data subset
# The plot is in a form of a heatmap, where pink blocks represent catch lots that have been [at least partially] sequenced
# the green blocks represent catch lots that have been submitted in the manifests but not yet sequenced
# Empty spaces represent either catch loths that have been samples but not yet submitted or catch lots that are missing completely [likely because parther joined us later in the project]

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

# Load data [put the right input directory and the file name]
dir_in   <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output"
pattern  <- "^BIOSCAN_100k_samples_corrected(\\d{4}-\\d{2}-\\d{2})\\.csv$" 
files <- list.files(dir_in, pattern = pattern, full.names = TRUE)
if (length(files) == 0) {
  stop("No dated BIOSCAN_working_data_*.csv files found in ", dir_in)
}
latest_file <- files[
  which.max(
    as.Date(sub(pattern, "\\1", basename(files)))
  )
]
message("Loading ", basename(latest_file))
cat("Loading ", basename(latest_file))
working_set <- read.csv(latest_file, stringsAsFactors = FALSE)

# Prepare data
plot_data <- working_set %>%
  dplyr::select(year, month, trap_name, catch_lot, plate, sequenced_lot_level) %>%
  distinct() %>% 
  mutate(# day = as.integer(day),
    month = as.integer(month),
    year = as.integer(year))

# Create plot grid
grid_data <- expand.grid(
  year = unique(plot_data$year),
  month = unique(plot_data$month),
  # day = unique(plot_data$day),
  trap_name = unique(plot_data$trap_name)
)

# Join with the plot data
plot_data_full <- grid_data %>%
  left_join(plot_data, by = c("year", "month", "trap_name"))

# Reorder traps by earliest day of sampling 
trap_order <- plot_data_full %>%
  filter(!is.na(catch_lot)) %>%
  group_by(trap_name) %>%
  summarise(earliest_month = min(year, na.rm = TRUE)) %>%
  arrange(earliest_month) %>%
  pull(trap_name)

plot_data_full <- plot_data_full %>%
  mutate(trap_name = factor(trap_name, levels = trap_order))

plot_data_full$month_factor <- as.factor(plot_data_full$month)

# Plot

completness <- ggplot(
  plot_data_full %>% filter(!(is.na(trap_name))) %>%
    mutate(
      # Only fill where catch_lot exists
      trap_name_fill = ifelse(is.na(catch_lot), NA, as.character(trap_name))
    ),
  aes(x = month_factor, y = trap_name, fill = sequenced_lot_level)
) +
  geom_tile(color = "grey70") +
  facet_grid(~year, switch = "y") +
  labs(
    x = "Month of sampling",
    y = "Trap name"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) + scale_fill_manual(values = c("#548B54", "#FFC0CB"),na.value = "white")

print(completness)

ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/completness2025_updated_catch_lot_sequenced.pdf",
  plot   = completness,
  device = "pdf",
  width  = 40, height = 45, units = "cm"
)
