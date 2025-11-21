
# This script calculates biodiversity metrics per catch lot
# The output is a series of pots and a csv file:
# biodiversity_catch_lot.csv
# With Simpson and Shannon diversity per catch lot 
# It's calculated based on all samples that have BINs 

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

cat(paste("Number of samples:", nrow(working_set),
          "\nNumber of sequenced samples:", nrow(working_set %>% filter(sequenced == TRUE)),
          "\nNumber of partners:", length(unique(working_set$partner)),
          "\nNumber of traps:", length(unique(working_set$trap_name)),
          "\nNumber of catch lots:", length(unique(working_set$sts_CATCH_LOT)),
          "\nNumber of (partially)sequenced catch lots:", length(working_set %>% 
                                                                  filter(sequenced_lot_level == TRUE) %>%
                                                                  pull(sts_CATCH_LOT) %>% unique())
))

# Filter only [fully/partially] sequenced catch lots 
working_set <- working_set %>% 
  filter(sequenced_lot_level == TRUE)

# Remove samples that have not been sequenced yet
working_set <- working_set %>% filter(bioscan_qc_sanger_qc_result != "None")

# Add information about the succes rate per catch lot
seq_lot_success <- working_set  %>%
  group_by(sts_CATCH_LOT) %>% 
  summarise(
    n_true  = sum(bioscan_qc_sanger_qc_result == "PASS", na.rm = TRUE),
    n_false = sum(bioscan_qc_sanger_qc_result != "PASS", na.rm = TRUE),
    .groups = "drop"
  ) %>% ungroup() %>% 
  mutate(
    seq_success_percentage = n_true*100/(n_true + n_false)
  ) %>% dplyr::select(sts_CATCH_LOT, seq_success_percentage)

working_set <- merge(working_set, seq_lot_success, by = "sts_CATCH_LOT")

ggplot(working_set, aes(x = seq_success_percentage)) +
  geom_histogram(binwidth = 1, color = "black", fill = "steelblue") +
  labs(
    x = "Sequencing success per catch lot (%)",
    y = "Count"
  ) +
  theme_minimal()

# Subset samples with BINs only [needed for abundance]
cat(paste("No. samples without BINs:", 
          working_set %>% filter(bioscan_qc_sanger_qc_result == "PASS") %>% filter(bold_bin_uri == "None") %>% nrow(),
          "\nPercentage of all samples:", 
          (working_set %>% filter(bioscan_qc_sanger_qc_result == "PASS") %>% filter(bold_bin_uri == "None") %>% nrow())*100/nrow(working_set), "%"
))
cat("Selecting only samples with BINs and that passed QC")
working_set_bins <- working_set %>% filter(bioscan_qc_sanger_qc_result == "PASS") %>% filter(bold_bin_uri != "None")

# length(unique(working_set$sts_CATCH_LOT))
# length(unique(working_set_bins$sts_CATCH_LOT))

# Shannon & Simpson [using vegan package for biodiveristy]
diversity_results <- working_set_bins %>%
  group_by(sts_CATCH_LOT) %>%
  summarise(
    Shannon = vegan::diversity(as.vector(table(bold_bin_uri)), index = "shannon"),
    Simpson = vegan::diversity(as.vector(table(bold_bin_uri)), index = "simpson"),
    n_unique_bins = n_distinct(bold_bin_uri),
    .groups = "drop"
  ) %>% ungroup()

# Add success rates per catch lot 
nrow(diversity_results)
diversity_results <- merge(diversity_results, seq_lot_success, by = "sts_CATCH_LOT")
nrow(diversity_results)

# Add functional diversity 
# Number of genera / number of families = taxonomic redundancy index
# Because:
# - many genera clustered in the same family = functionally redundant, low phylogenetic breadth
# - fewer genera per family (i.e., many families represented) = functionally broader, higher ecosystem integrity

# Should we filter out non-arthropods?????

taxo_table <- working_set_bins %>%
  group_by(sts_CATCH_LOT) %>%
  summarise(
    n_families = n_distinct(bold_family),
    n_genera   = n_distinct(bold_genus),
    n_bins     = n_distinct(bold_bin_uri),
    .groups = "drop"
  )

taxo_table <- taxo_table %>%
  mutate(
    genera_per_family = n_genera / n_families
  ) %>%
  mutate(
    bin_per_family = n_bins / n_families
  ) %>%
  mutate(
    bin_per_genus = n_bins / n_genera
  )

# Merge 

nrow(diversity_results)
diversity_results <- merge(diversity_results, taxo_table, by = "sts_CATCH_LOT")
nrow(diversity_results)

# Save file
file_out <- sprintf(
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output/04_biodiversity_catch_lot_%s.csv",
  today_stamp
)
write.csv(diversity_results, file_out, row.names = FALSE)

# Catch lot success results
all_CL <- ggplot(working_set %>% dplyr::select(sts_CATCH_LOT, seq_success_percentage) %>% unique(),
                 aes(x = seq_success_percentage)) +
  geom_histogram(binwidth = 1, color = "black", fill = "steelblue") +
  labs(
    x = "Sequencing success per catch lot (%)",
    y = "Count"
  ) +
  theme_minimal()

passedQC <- ggplot(diversity_results, aes(x = seq_success_percentage)) +
  geom_histogram(binwidth = 1, color = "black", fill = "pink") +
  labs(
    x = "Sequencing success per catch lot (%)",
    y = "Count"
  ) +
  theme_minimal()

all_CL + passedQC

# Diversity plots 
diversity_results$n_unique_bins <- factor(
  diversity_results$n_unique_bins,
  levels = min(diversity_results$n_unique_bins):max(diversity_results$n_unique_bins),
  ordered = TRUE
)

shannon_distribution <- ggplot(diversity_results,
                               aes(x = Shannon)) +
  geom_histogram(binwidth = 0.1, boundary = 0,
                 fill = "#5dc3e8", colour = "white") +
  scale_x_continuous(
    name   = "Shannon Index",
    breaks = seq(0, max(diversity_results$Shannon), 1)) +
  #expand = expansion(mult = c(0, 0.02))) +
  labs(y = "Frequency [catch lot]") +
  coord_cartesian(clip = "off") +   
  theme_classic(base_size = 10)

simpson_distribution <- ggplot(diversity_results,
                               aes(x = Simpson)) +
  geom_histogram(binwidth = 0.01, boundary = 0,
                 fill = "#5dc3e8", colour = "white") +
  scale_x_continuous(
    name   = "Simpson Index",
    breaks = seq(0, 1, 0.1)) +
  #expand = expansion(mult = c(0, 0.02))) +
  labs(y = "Frequency [catch lot]") +
  coord_cartesian(clip = "off") +   
  theme_classic(base_size = 10)

diversity_hist <- shannon_distribution/simpson_distribution

ggsave(      
  "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/shannon_simpson_lot_histogram.pdf",
  plot   = diversity_hist,
  device = "pdf",
  width  = 25, height = 24, units = "cm"
)

out_file <- "/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots/diversity_results_correlations.pdf"
pdf(out_file, width = 10, height = 10)
plot(diversity_results[-1])
dev.off()  
