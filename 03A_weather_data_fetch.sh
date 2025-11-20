# Get the weather data directly from the MET office

# Go to https://accounts.ceda.ac.uk/realms/ceda/account/
# Log in or create an account
# Navigate to Archive access -> Access Token -> Create new access token 
# Copy the new token and paste below: 

# Run as an interactive job: 

bsub -Is -n 4 -R "select[mem>20000] rusage[mem=20000] span[hosts=1]" -M 20000 -G team222 bash

export CEDA_TOKEN=" "

# Navigate to the weather directory
cd lustre/scratch126/tol/teams/lawniczak/projects/bioscan/100k_paper/met_weather


# Get rain data

curl -s -H "Authorization: Bearer $CEDA_TOKEN" \
  https://dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.1.ceda/5km/rainfall/day/v20250415/ \
  | grep -oP 'href="\K[^"]+' \
  | grep -E '2021|2022|2023|2024|2025' \
  | sed 's|^|https://dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.1.ceda/5km/rainfall/day/v20250415/|' \
  > METdata_rain_download_list.txt

wget --header "Authorization: Bearer $CEDA_TOKEN" -i METdata_rain_download_list.txt

rm rainfall_hadukgrid_uk_5km_day_1*
rm rainfall_hadukgrid_uk_5km_day_200*
rm rainfall_hadukgrid_uk_5km_day_2012*

mv *.nc ./rain_data/


# Get max temperature data

curl -s -H "Authorization: Bearer $CEDA_TOKEN" \
  https://dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.1.ceda/5km/tasmax/day/v20250415/ \
  | grep -oP 'href="\K[^"]+' \
  | grep -E '2021|2022|2023|2024|2025' \
  | sed 's|^|https://dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.1.ceda/5km/tasmax/day/v20250415/|' \
  > METdata_tasmax_download_list.txt

wget --header "Authorization: Bearer $CEDA_TOKEN" -i METdata_tasmax_download_list.txt

rm tasmax_hadukgrid_uk_5km_day_1*
rm tasmax_hadukgrid_uk_5km_day_200*
rm tasmax_hadukgrid_uk_5km_day_2012*

mv *.nc ./tasmax_data/


# Get min temperature data 

curl -s -H "Authorization: Bearer $CEDA_TOKEN" \
  https://dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.1.ceda/5km/tasmin/day/v20250415/ \
  | grep -oP 'href="\K[^"]+' \
  | grep -E '2021|2022|2023|2024|2025' \
  | sed 's|^|https://dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.1.ceda/5km/tasmin/day/v20250415/|' \
  > METdata_tasmin_download_list.txt

wget --header "Authorization: Bearer $CEDA_TOKEN" -i METdata_tasmin_download_list.txt

rm tasmin_hadukgrid_uk_5km_day_1*
rm tasmin_hadukgrid_uk_5km_day_200*
rm tasmin_hadukgrid_uk_5km_day_2012*

mv *.nc ./tasmin_data/

