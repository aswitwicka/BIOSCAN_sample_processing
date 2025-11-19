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
This directory contains all CEH (land cover) and MET office (weather) maps as well as general files required to visualise UK maps.

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














