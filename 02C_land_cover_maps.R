# Continnuation of habitat type description per trap
# The goal is to calculate and compare habitat diversty metrics across different radius sies around the trap 
# The output files are a series of plots and csv files

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
  "knitr", "patchwork", "rnaturalearth", "rnaturalearthdata", 
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

# Load RDS from the previous script 
working_sets <- readRDS("/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02B_working_sets_radius.rds")

# Today stamp
today_stamp <- format(Sys.Date(), "%Y-%m-%d") 

# Combine dsf
standardize_cols <- function(df) {
  # GB 
  if ("gblcm2024_10m_1" %in% names(df)) {
    names(df)[names(df) == "gblcm2024_10m_1"] <- "lcm2024_value"
    names(df)[names(df) == "gblcm2024_10m_2"] <- "lcm2024_fraction"
  }
  # NI 
  if ("nilcm2024_10m_1" %in% names(df)) {
    names(df)[names(df) == "nilcm2024_10m_1"] <- "lcm2024_value"
    names(df)[names(df) == "nilcm2024_10m_2"] <- "lcm2024_fraction"
  }
  df
}
working_sets <- lapply(working_sets, standardize_cols)

obj_names <- names(working_sets)
buffer_sizes <- sub("_buffer.*", "", obj_names)
groups <- split(obj_names, buffer_sizes)

combined <- lapply(groups, function(nm) {
  dplyr::bind_rows(working_sets[nm])
})
names(combined) <- paste0("combined_buffer_", names(combined))

# Extract dfs
list2env(combined, envir = .GlobalEnv)

# Extract unique radius IDs once
radius_ids <- unique(buffer_sizes)

######### Agriculture ratio in all traps ######### 

buffer_objs <- ls(pattern = "combined_buffer_*")  
# Loop
meta_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i] 
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")
  meta_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      Arable = sum(habitat_type %in% c("Arable & horticulture"), na.rm = TRUE),
      Other = sum(!(habitat_type %in% c("Arable & horticulture")), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      total_pixels = Arable + Other,
      !!paste0("agriculture_ratio_", base_buf_name) :=
        (Arable) / total_pixels
    ) %>% 
    dplyr::select(trap_name,
                  starts_with("agriculture_ratio_"))
}
# Join
agriculture_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                           meta_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_agriculture_%s.csv",
  today_stamp
)
write.csv(agriculture_meta, file_out, row.names = FALSE)
# Heatmap
agri_long <- agriculture_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    buffer = sub("^agriculture_ratio_", "", buffer),
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
# Order traps by their mean ratio
agri_long <- agri_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
agri_heatmap <- ggplot(agri_long, aes(x = buffer,
                                      y = fct_rev(trap_name), # puts highest mean at top
                                      fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Arable & horticulture\nratio",
    option = "C",    
    direction = -1,
    limits = c(0, 1) 
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),   # long labels fit
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/agriculture_heatmap.pdf",
  plot = agri_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

# Correlation plot
out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/agriculture_correlation.pdf"
pdf(out_file, width = 10, height = 10)
plot(agriculture_meta %>% dplyr::select(-trap_name))
dev.off() 

######### Natural grasslands ######### 

# Loop
meta_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i] 
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")
  meta_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      Grassland_natural = sum(habitat_type %in% 
                                c("Acid grassland", "Natural grassland", "Calcareous grassland", "Fen, marsh & swamp"), na.rm = TRUE),
      Other = sum(!(habitat_type %in% 
                         c("Acid grassland", "Natural grassland", "Calcareous grassland", "Fen, marsh & swamp")), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      total_pixels = Grassland_natural + Other,
      !!paste0("natural_grasslapn_ratio_", base_buf_name) :=
        (Grassland_natural) / total_pixels
    ) %>% 
    dplyr::select(trap_name,
                  starts_with("natural_grasslapn_ratio_"))
}
# Join
natural_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                       meta_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_natural_grasslands_meta_%s.csv",
  today_stamp
)
write.csv(natural_meta, file_out, row.names = FALSE)
# Heatmap
natural_long <- natural_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    buffer = sub("^natural_grasslapn_ratio_", "", buffer),
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
# Order traps by their mean ratio
natural_long <- natural_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
natural_heatmap <- ggplot(natural_long, aes(x = buffer,
                                            y = fct_rev(trap_name),   
                                            fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Natural grasslands\nratio",
    option = "C",    
    direction = -1,
    limits = c(0, 1) 
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),   # long labels fit
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/natural_grasslands_heatmap.pdf",
  plot = natural_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

# Correlation plot
out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/natural_grasslands_correlation.pdf"
pdf(out_file, width = 10, height = 10)
plot(natural_meta %>% dplyr::select(-trap_name))
dev.off() 

######### Forest ######### 

# Loop
meta_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i] 
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")
  meta_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      Broadleaf = sum(habitat_type == "Broadleaf woodland", na.rm = TRUE),
      Coniferous = sum(habitat_type == "Coniferous woodland",    na.rm = TRUE),
      Other = sum(!habitat_type %in% c("Coniferous woodland", "Broadleaf woodland"), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      total_pixels = Broadleaf + Coniferous + Other,
      !!paste0("forest_ratio_", base_buf_name) :=
        (Broadleaf + Coniferous) / total_pixels
    ) %>% 
    dplyr::select(trap_name,
                  starts_with("forest_ratio_"))
}
# Join
forest_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                      meta_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_forest_%s.csv",
  today_stamp
)
write.csv(forest_meta, file_out, row.names = FALSE)
# Heatmap
forest_long <- forest_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    ## drop the “urbanisation_ratio_” prefix
    buffer = sub("^forest_ratio_", "", buffer),
    ## order buffers numerically (X100, X500, X1000, …)
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
# Order traps by their mean ratio
forest_long <- forest_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
forest_heatmap <- ggplot(forest_long, aes(x = buffer,
                                          y = fct_rev(trap_name),   # puts highest mean at top
                                          fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Forest (Broadleaf & Coniferous)\nratio",
    option = "C",    
    direction = -1,
    limits = c(0, 1) 
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),   # long labels fit
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/forest_heatmap.pdf",
  plot = forest_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

# Correlation plot
out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/forest_correlation.pdf"
pdf(out_file, width = 10, height = 10)
plot(forest_meta %>% dplyr::select(-trap_name))
dev.off()  

######### Urbanisation levels ######### 

# Loop
meta_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i] 
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")
  meta_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      Suburban = sum(habitat_type == "Suburban", na.rm = TRUE),
      Urban = sum(habitat_type == "Urban",    na.rm = TRUE),
      Other = sum(!habitat_type %in% c("Suburban", "Urban"), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      total_pixels = Suburban + Urban + Other,
      !!paste0("urbanisation_ratio_", base_buf_name) :=
        (Suburban + Urban) / total_pixels
    ) %>% 
    dplyr::select(trap_name,
                  starts_with("urbanisation_ratio_"))
}
# Join
urbanisation_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                            meta_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_urbanisation_%s.csv",
  today_stamp
)
write.csv(urbanisation_meta, file_out, row.names = FALSE)
# Heatmap
urban_long <- urbanisation_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    buffer = sub("^urbanisation_ratio_", "", buffer),
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
urban_long <- urban_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
urbanisation_heatmap <- ggplot(urban_long, aes(x = buffer,
                                               y = fct_rev(trap_name),  
                                               fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Urbanisation\nratio",
    option = "C",    
    direction = -1,
    limits = c(0, 1) 
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),   # long labels fit
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/urbanisation_heatmap.pdf",
  plot = urbanisation_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

# Correlation plot
out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/urbanisation_correlation.pdf"
pdf(out_file, width = 10, height = 10)
plot(urbanisation_meta %>% dplyr::select(-trap_name))
dev.off()

######### Improved grasslands ######### 

# Loop
meta_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i] 
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")
  meta_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      Improved_grassland = sum(habitat_type == "Improved grassland", na.rm = TRUE),
      Other = sum(!habitat_type %in% c("Improved grassland"), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      total_pixels = Improved_grassland + Other,
      !!paste0("improved_grassland_ratio_", base_buf_name) :=
        (Improved_grassland) / total_pixels
    ) %>% 
    dplyr::select(trap_name,
                  starts_with("improved_grassland_ratio_"))
}
# Join
improved_grassland_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                                  meta_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_improved_grassland_%s.csv",
  today_stamp
)
write.csv(improved_grassland_meta, file_out, row.names = FALSE)
# Heatmap
grassland_long <- improved_grassland_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    buffer = sub("^improved_grassland_ratio_", "", buffer),
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
# Order traps by their mean ratio
grassland_long <- grassland_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
grassland_heatmap <- ggplot(grassland_long, aes(x = buffer,
                                                y = fct_rev(trap_name),  
                                                fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Improved grasslands\nratio",
    option = "C",    
    direction = -1,
    limits = c(0, 1) 
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),   # long labels fit
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/improved_grassland_heatmap.pdf",
  plot = grassland_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

# Correlation plot
out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/improved_grassland_correlation.pdf"
pdf(out_file, width = 10, height = 10)
plot(improved_grassland_meta %>% dplyr::select(-trap_name))
dev.off()

######### Coastal ######### 

# Loop
meta_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i] 
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")
  meta_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      Coastal = sum(habitat_type %in% 
                      c("Saltwater", "Supralittoral sediment", "Supralittoral rock", "Littoral rock", "Littoral sediment", "Saltmarsh"), na.rm = TRUE),
      Other = sum(!(habitat_type %in% 
                         c("Saltwater", "Supralittoral sediment", "Supralittoral rock", "Littoral rock", "Littoral sediment", "Saltmarsh")), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      total_pixels = Coastal + Other,
      !!paste0("coastal_ratio_", base_buf_name) :=
        (Coastal) / total_pixels
    ) %>% 
    dplyr::select(trap_name,
                  starts_with("coastal_ratio_"))
}
# Join
coastal_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                       meta_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_coastal_meta_%s.csv",
  today_stamp
)
write.csv(coastal_meta, file_out, row.names = FALSE)

# Heatmap
costal_long <- coastal_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    buffer = sub("^coastal_ratio_", "", buffer),
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
# Order traps by their mean ratio
costal_long <- costal_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
coastal_heatmap <- ggplot(costal_long, aes(x = buffer,
                                           y = fct_rev(trap_name),  
                                           fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Coastal landscape\nratio",
    option = "C",    
    direction = -1,
    limits = c(0, 1) 
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),   # long labels fit
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/coastal_heatmap.pdf",
  plot = coastal_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

# Correlation plot
out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/coastal_correlation.pdf"
pdf(out_file, width = 10, height = 10)
plot(coastal_meta %>% dplyr::select(-trap_name))
dev.off()

######### Bog & heather ######### 

# Loop
meta_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i] 
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")
  meta_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      Heather_mountain_bog = sum(habitat_type %in% 
                                   c("Heather", "Heather grassland", "Bog", "Inland rock"), na.rm = TRUE),
      Other = sum(!(habitat_type %in% 
                         c("Heather", "Heather grassland", "Bog", "Inland rock")), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      total_pixels = Heather_mountain_bog + Other,
      !!paste0("heather_mountain_bog_ratio_", base_buf_name) :=
        (Heather_mountain_bog) / total_pixels
    ) %>% 
    dplyr::select(trap_name,
                  starts_with("heather_mountain_bog_ratio_"))
}
# Join
heather_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                       meta_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_heather_mountain_bog_meta_%s.csv",
  today_stamp
)
write.csv(heather_meta, file_out, row.names = FALSE)

# Heatmap
heather_long <- heather_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    buffer = sub("^heather_mountain_bog_ratio_", "", buffer),
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
# Order traps by their mean ratio
heather_long <- heather_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
heather_heatmap <- ggplot(heather_long, aes(x = buffer,
                                            y = fct_rev(trap_name),
                                            fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Heather, mountain, bog\nratio",
    option = "C",    
    direction = -1,
    limits = c(0, 1) 
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/heather_mountain_bog_heatmap.pdf",
  plot = heather_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

# Correlation plot
out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/heather_mountain_bog_correlation.pdf"
pdf(out_file, width = 10, height = 10)
plot(heather_meta %>% dplyr::select(-trap_name))
dev.off()

######### Freshwater ######### 

# Loop
meta_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i] 
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")
  meta_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      Freshwater = sum(habitat_type %in% c("Freshwater"), na.rm = TRUE),
      Other = sum(!(habitat_type %in% c("Freshwater")), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      total_pixels = Freshwater + Other,
      !!paste0("freshwater_ratio_", base_buf_name) :=
        (Freshwater) / total_pixels
    ) %>% 
    dplyr::select(trap_name,
                  starts_with("freshwater_ratio_"))
}
# Join
freshwater_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                          meta_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_freshwater_%s.csv",
  today_stamp
)
write.csv(freshwater_meta, file_out, row.names = FALSE)
# Heatmap
freshwater_long <- freshwater_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    buffer = sub("^agriculture_ratio_", "", buffer),
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
# Order traps by their mean ratio
freshwater_long <- freshwater_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
fresh_heatmap <- ggplot(freshwater_long, aes(x = buffer,
                                             y = fct_rev(trap_name), 
                                             fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Freshwater \nratio",
    option = "C",    
    direction = -1,
    limits = c(0, 1) 
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),   # long labels fit
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/freshwater_heatmap.pdf",
  plot = fresh_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

# Correlation plot
out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/freshwater_correlation.pdf"
pdf(out_file, width = 10, height = 10)
plot(freshwater_meta %>% dplyr::select(-trap_name))
dev.off() 

######### Unique habitat type ######### 
# Count unique land cover types per working_trap
land_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i]
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")   
  land_list[[i]] <- df %>% 
    group_by(trap_name) %>% 
    summarise(
      !!paste0("unique_land_types_", base_buf_name) :=
        n_distinct(habitat_type),
      .groups = "drop"
    )
}
# Combine unique-type counts from every buffer
landcover_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                         land_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_unique_land_types_%s.csv",
  today_stamp
)
write.csv(landcover_meta, file_out, row.names = FALSE)

# Heatmap
landcover_long <- landcover_meta %>% 
  pivot_longer(-trap_name,
               names_to = "buffer",
               values_to = "ratio") %>% 
  mutate(
    buffer = sub("^urbanisation_ratio_", "", buffer),
    buffer = fct_reorder(buffer,
                         as.numeric(stringr::str_extract(buffer, "\\d+")))
  )
# Order traps by their mean ratio
landcover_long <- landcover_long %>% 
  group_by(trap_name) %>% 
  mutate(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(trap_name = fct_reorder(trap_name, mean_ratio)) %>% 
  dplyr::select(-mean_ratio)
# Plot
landscape_heatmap <- ggplot(landcover_long, aes(x = buffer,
                                                y = fct_rev(trap_name), 
                                                fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "No. Unique\nhabitat types",
    option = "C",    
    direction = -1
  ) +
  labs(
    x = "Buffer size",
    y = "Trap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),   # long labels fit
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/unique_landscape_heatmap.pdf",
  plot = landscape_heatmap,
  device = "pdf",
  width = 11, height = 20, units = "cm"
)

######### Dominant habitat type ######### 
# Calculate dominant category per working_trap based on the sum of fractions
dom_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i]
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")  
  dom_list[[i]] <- df %>% 
    group_by(trap_name, habitat_type) %>% 
    summarise(total_pixels = sum(fraction, na.rm = TRUE), .groups = "drop") %>% 
    group_by(trap_name) %>% 
    mutate(percentage = 100 * total_pixels / sum(total_pixels)) %>% 
    slice_max(percentage, n = 1, with_ties = FALSE) %>% 
    ungroup() %>% 
    transmute(
      trap_name,
      !!paste0("dominant_habitat_", base_buf_name) := habitat_type,
      !!paste0("dominant_pct_",     base_buf_name) := percentage
    )
}
# Combine the dominant-habitat info from all buffers
dominant_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                        dom_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_dominant_land_type_%s.csv",
  today_stamp
)
write.csv(dominant_meta, file_out, row.names = FALSE)

######### Shannon Diversity Index ######### 
# Weighted richness + evenness; Sensitive to rare types
# Range: 0 to ln(S), where S = number of categories (habitat types)
# For example, if a trap has 5 habitat types: max H′ = ln(5) ≈ 1.61
# This is called Shannon’s Equitability Index or Pielou’s Evenness (normalised):
# shannon_index_normalised = H′ / ln(S)

######### Simpsons Diversity Index ######### 
# Dominance of most abundant types; Sensitive to dominant types
# Formula (1 - sum(p^2)) gives the Gini-Simpson Index
# Range: 0 to (1 - 1/S), which approaches 1 as richness and evenness increase
# 0 = only one habitat type
# 1 = infinite number of perfectly even types (theoretical)

######### Shannon & Simpson habitat diversity per buffer ######### 
div_list <- vector("list", length(buffer_objs))
for (i in seq_along(buffer_objs)) {
  obj_name <- buffer_objs[i]
  df <- get(obj_name)
  bits <- strsplit(obj_name, "_")[[1]]
  base_buf_name <- paste(bits[2:3], collapse = "_")   
  div_list[[i]] <- df %>% 
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
}
# Combine all diversity indices
diversity_meta <- Reduce(function(x, y) full_join(x, y, by = "trap_name"),
                         div_list)
# Save
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/02_LandCover_diversity_indices_%s.csv",
  today_stamp
)
write.csv(diversity_meta, file_out, row.names = FALSE)

# Plot

radii <- names(diversity_meta) %>%   
  stringr::str_extract("(?<=_buffer_)\\d+$") %>%
  na.omit() %>%
  unique() %>%
  as.numeric() %>%                                 
  sort()
out_dir <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/"
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (r in radii) {
  sh_col <- paste0("shannon_buffer_",  r)
  sim_col <- paste0("simpson_buffer_",  r)
  max_val <- max(
    max(diversity_meta[[sh_col]], na.rm = TRUE),
    max(diversity_meta[[sim_col]], na.rm = TRUE)
  )
  min_val <- min(
    min(diversity_meta[[sh_col]], na.rm = TRUE),
    min(diversity_meta[[sim_col]], na.rm = TRUE)
  )
  p <- ggplot(diversity_meta,
              aes_string(x = sh_col, y = sim_col)) +  
    geom_point(color = "#d69ee8", size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", color = "#5d1075",
                se = TRUE, linetype = "dashed") +
    labs(
      x = "Shannon Diversity Index",
      y = "Simpson Diversity Index",
      title = paste("Shannon vs Simpson:", r, "m buffer")
    ) +
    theme_classic() +
    coord_equal(xlim = c(min_val, max_val), ylim = c(min_val, max_val))
  ggsave(
    filename = file.path(out_dir,
                         paste0("indice_correlation_X", r, "_buffer.pdf")),
    plot = p,
    device = "pdf",
    width = 10, height = 10, units = "cm"
  )
}

# Correlate per index across radi - Shannon
# Build every unique combination of two different radii 
radius_pairs <- combn(radii, 2, simplify = FALSE)   
for (pair in radius_pairs) {
  r1 <- pair[1]
  r2 <- pair[2]
  col_x <- paste0("shannon_buffer_", r1)
  col_y <- paste0("shannon_buffer_", r2)
  max_val <- max(
    max(diversity_meta[[col_x]], na.rm = TRUE),
    max(diversity_meta[[col_y]], na.rm = TRUE)
  )
  min_val <- min(
    min(diversity_meta[[col_x]], na.rm = TRUE),
    min(diversity_meta[[col_y]], na.rm = TRUE)
  )
  p <- ggplot(diversity_meta, 
              aes_string(x = col_x, y = col_y)) +
    geom_point(color = "#d69ee8", size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", color = "#5d1075",
                se = TRUE, linetype = "dashed") +
    labs(
      x = paste0("Shannon Diversity Index (", r1, "-m buffer)"),
      y = paste0("Shannon Diversity Index (", r2, "-m buffer)"),
      title = paste("Shannon correlation:", r1, "m vs", r2, "m buffers")
    ) + coord_fixed(ratio = 1) +
    theme_classic() +
    coord_equal(xlim = c(min_val, max_val), ylim = c(min_val, max_val))
  ggsave(
    filename = file.path(out_dir,
                         paste0("shannon_correlation_X", r1,
                                "_vs_X", r2, "_buffer.pdf")),
    plot = p,
    device = "pdf",
    width = 10, height = 10, units = "cm"
  )
}

# Shannon per trap/radius 
div_long <- diversity_meta %>%
  pivot_longer(
    cols = starts_with(c("shannon", "simpson")),
    names_to = c("index", "buffer"),
    names_pattern = "(shannon|simpson)_buffer_(\\d+)",
    values_to = "value"
  )

# Rescale Shannon within each trap
div_long_scaled <- div_long %>%
  filter(index == "shannon") %>%
  group_by(trap_name) %>%
  mutate(
    value_scaled = (value - min(value, na.rm = TRUE)) /
      (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
  ) %>%
  ungroup()

# Order traps by variability across large buffers
var_order <- div_long_scaled %>%
  filter(buffer %in% radii) %>%
  group_by(trap_name) %>%
  summarise(
    rng = diff(range(value_scaled, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(rng)) %>%
  pull(trap_name)

# Plots (these plots really make sense only when very large sizes like 1000, 2500, 5000, 10000m are included)
shannon_heatmap <- ggplot(
  div_long_scaled %>%
    mutate(
      trap_name = factor(trap_name, levels = var_order),
      buffer = factor(buffer, levels = c("25", "50","100","500","1000"))
    ),
  aes(x = trap_name, y = buffer, fill = value_scaled)
) +
  geom_tile() +
  scale_fill_viridis_c(name = "Shannon [scaled by trap]", option = "magma") +
  labs(
    x = "Trap",
    y = "Buffer size [m]"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )
# ggsave(
#   "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/habitat_complexity/output/processing_plots/shannon_variability_heatmap_2024.pdf",
#   plot = shannon_heatmap,
#   device = "pdf",
#   width = 25, height = 15, units = "cm"
# )

shannon_lines <- ggplot(div_long_scaled, 
                        aes(x = factor(trap_name, levels = var_order), y = value, 
                            colour = factor(buffer, levels = c("25", "50","100","500","1000")), 
                            group = factor(buffer, levels = c("25", "50","100","500","1000")))) +
  geom_point(size = 2, alpha = 0.7, position = position_dodge(width = 0.5)) +
  geom_line(aes(group = interaction(index, buffer)),
            alpha = 0.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(
    "25" = "lightpink",
    "50" = "#ce6bed",
    "100" = "#933dad",
    #"500" = "#5f2570",
    # "2500" = "#0c010f", 
    "500" = "#076378",
    "1000" = "#46a3b8" 
    #"10000" = "#84cfe0"
  )) +
  labs(
    x = "Trap",
    y = "Shannon index",
    colour = "Buffer size [m]"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "bottom"
  )
# ggsave(
#   "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/habitat_complexity/output/processing_plots/shannon_variability_lineplot_2024.pdf",
#   plot = shannon_lines,
#   device = "pdf",
#   width = 25, height = 15, units = "cm"
# )

# The most and least diverse traps 

summary_stats <- div_long %>%
  group_by(trap_name, index) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    .groups = "drop"
  )

summary_stats_shannon <- summary_stats %>% filter(index == "shannon")
summary_stats_simpson <- summary_stats %>% filter(index == "simpson")

order_shannon <- summary_stats %>%
  filter(index == "shannon") %>%
  arrange(median_value) %>%
  pull(trap_name)

mean_index <- ggplot(summary_stats, 
                     aes(x = factor(trap_name, levels = order_shannon), y = mean_value, colour = index))+
  geom_point(size = 2, alpha = 0.7, position = position_dodge(width = 0.5)) +
  geom_line(aes(group = interaction(index)),
            alpha = 0.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(
    "shannon" = "#ce6bed",
    "simpson" = "#84cfe0"
  )) +
  labs(
    x = "Trap",
    y = "Mean index values",
    colour = "Diversity indices:"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "bottom"
  )

median_index <- ggplot(summary_stats, 
                       aes(x = factor(trap_name, levels = order_shannon), y = median_value, colour = index))+
  geom_point(size = 2, alpha = 0.7, position = position_dodge(width = 0.5)) +
  geom_line(aes(group = interaction(index)),
            alpha = 0.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(
    "shannon" = "#ce6bed",
    "simpson" = "#84cfe0"
  )) +
  labs(
    x = "Trap",
    y = "Median index values",
    colour = "Diversity indices:"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "bottom"
  )

indice_ordered_traps <- median_index / mean_index

ggsave(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/habitat_indice_ordered_traps_2024.pdf",
  plot = indice_ordered_traps,
  device = "pdf",
  width = 25, height = 15, units = "cm"
)
