## Scripts to download and process BIOSCAN samples 
🐛 🪰 🐝 🪲 

0. Download data
1. Subset samples
2. Calculate habitat metrics per trap
3. Fetch weather information per catch lot
4. Visualise which catch lots have not been sequenced yet
5. Calculate Shannon and Simpson indices per catch lot

### All required flat files are in: <br>
/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/required_files <br> 
• collection_time_codes.csv (time codes translated to numeric values) <br>
• mozz_to_partner.csv (older mozz plates translated to partner codes) <br>
• trap_to_partner.csv (catch lot location translated to trap name)  <br><br>
& /lustre/scratch126/tol/teams/lawniczak/projects/bioscan/bioscan_qc/mbrave_batch_data <br>
This directory stores all mBRAVE files - information of all samples that have been sequenced. <br><br>
& /lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/maps <br>
This directory contains all CEH (land cover) and MET office (weather) maps as well as general files required to visualise UK maps. <br><br>
These files must me manually updated when new MET or CEH data appears or when new partners join BIOSCAN. 

### All code should be run within:
/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/code

### All outputs are present in: 
/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/output <br>
& <br>
/lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/plots

### 00_manifest_fetch.sh 
To run execute: 
```bash
bsub < 00_manifest_fetch.sh
```
The code will create sts_manifests_[date].tsv file which contains all samples that are present in manifests submitted by the partners (sequenced and not sequenced)

### 01_subset_data.R
To run interactively connect to R studio and load the file: <br>
```bash
module load HGI/softpack/users/aw43/BOLDconnectR_bioscan/2
module load HGI/softpack/users/aw43/aw43_bioscan_habitat_complexity-2/3
module load rstudio
rstudio start --cpus 10 -M 100000 --queue normal --home /lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/code --pwd /lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/code -g team222
```
This can be done with all the R scripts in this repository. Remember to load in the specified softpack environments. <br>
Or submit as a job:
```bash
bsub < 01_subset_data.sh
```
This script takes the output of 00_manifest_fetch.sh and subsets the BIOSCAN data for further processing. It removes control samples and samples from all non-BIOSCAN partners. <br>
<b>NOTE:</b> if any new non-BIOSCAN samples are being processed and go through the standard BIOSCAN QC, these partner codes must be added manually within the 01_subset_data.R script. <br>
Currently included: "BGEP", "SNST", "JARO", "BGKU", "BSN", "BGEG", "OXHP", "POMS", "AYDI", "BGPT", "-BGE", "-TOL" <br>
The script retaines samples caught using Malaise trap only and from years 2021 onwards. <br>
Additional columns added to the data include: trap name, partner name, region, plate, day, month, year, 24h sampling selection, if any samples from a given plate / catch lot have been sequenced [TRUE/FALSE]. <br>
The script also corrects previously recognised mistakes in the coordinates of some traps. If these have been already corrected in the Portal, the script won't introduce any changes. <br> 
The script assesses completnes of the data using the mBRAVE input files to check if any BIOSCAN samples are missing the QC output information. <br>
<b>NOTE:</b> BOLDconnectR library requires an API key within the bold.apikey() function that you can obtain from your BOLD account. <br><br>
The output file: BIOSCAN_100k_samples_corrected[date].csv

### 02A_land_cover_maps.R

Do not run this script interactively because it would be very very slow, instead submit as a job:
```bash
bsub < 02A_land_cover_maps.sh
```
This script generates habitat-type summaries around BIOSCAN trap locations using 10m resolution CEH 2024 land cover maps. <br>
NI and GB are processed separately due to the use of different CEH maps. <br>
For each trap, the script: <br>
• Loads trap coordinates from a flat file trap_to_partner.csv <br>
• Loads land cover raster layers for Great Britain (GB) and Northern Ireland (NI) <br>
• Creates circular buffers of multiple radii (25, 50, 100, 500, 1000 m) - these sizes are currently specified in the script, please edit if you require bigger sizes <br>
• Extracts all raster pixels intersecting each buffer <br>
• Converts raster class codes into habitat labels <br>
• Calculates the number of pixels per trap per buffer <br>
• Saves one output csv per buffer size, separately for GB and NI <br><br>
The output files: [buffer size]_ buffer_[GBL or NIL]_ 2024_[date].csv

### 02B_land_cover_maps.R

To run from farm submit:
```bash
bsub < 02B_land_cover_maps.sh
```

This script loads the most recent csv files produced by 02A_land_cover_maps.R, extracts metadata and stores all datasets in a single rds file preparing them for downstream analyses. <br><br>
The output files: 02B_working_sets_radius.rds

### 02C_land_cover_maps.R

To run from farm submit:
```bash
bsub < 02C_land_cover_maps.sh
```
Or run via RStudio <br><br>

This script calculates habitat diversity indices and ratios of each habitat type for each BIOSCAN trap. The script uses the 02B_working_sets_radius.rds file as input. <br><br>
Habitat ratios are calculated as ratio = habitat_pixels / total_pixels_within_buffer for: <br>
• Arable & horticulture <br>
• Semi-natural grasslands <br>
• Forest (broadleaf + coniferous) <br>
• Urban & suburban & gardens <br>
• Improved grasslands <br>
• Coastal habitats <br>
• Heather / mountain / bog complexes <br>
• Freshwater <br>
<b>The script also returns: </b> <br>
• Number of unique habitat types per trap <br>
• Dominant habitat type and its percentage <br>
• Shannon diversity index (richness + evenness) <br>
• Simpson diversity index (evenness / dominance) <br>
Each metric is exported to the output directory as a standalone CSV summarising traps & buffers.<br>
<b>Generated plots include: </b> <br>
• Heatmaps showing spatial patterns of each habitat ratio across all traps <br>
• Correlation plots among buffers for each metric <br>
• Shannon–Simpson scatterplots per buffer radius <br>
• Cross-radius correlations of Shannon diversity <br>
• Summary plots ranking traps by mean and median diversity indices <br>

### 02D_land_cover_maps.R

### 03_weather_data_fetch.sh
(is this python??)

### 03_weather_data.R

### 04_visualise_present_catch_lots.R

### 05_biodiversity_per_catch_lot.R
Add functional diversity here 
Success per catch lot is calculated here and returned in the output data frame 

### 07_biodiversity_vs_habitat_models.R
for all models remember to subset only 24h and catch lots with good performance - depending on the analysis
also subset the months based on the 24h sampling distribution (winter excluded)
### 08_bins_vs_habitat_models.R

### 09_temporal_turnover.R
how to deal with the winter 24h+ samples??
### 10_taxonomy_plots.R

### 11_nbn_vs_bioscan.R
higher-level comparisons 















