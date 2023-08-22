# CV-MMI-WNV-Forecasting

This repository contains the processed data, code, and output figures/spreadsheets used for AGU submission 2023GH000855R, A spatially resolved and environmentally informed forecast model of West Nile virus in Coachella Valley, California., which has been submitted for possible publication in GeoHealth.

Within the ```data``` folder are the estimated annual WNV infection rates from 2006-2022 by NLDAS grid cell in Coachella Valley and the monthly temperature and evapotranspiration estimates from 2005-2022 by NLDAS grid cell. 

Within the ```model_outputs``` folder are the model parameters for models trained on data from 2006-2018 and 2006-2021. These files were generated using ```cv_nldas_all_years.R``` and ```cv_nldas_lootcv.R``` found in the root directory. 

The ```shapefiles``` folder contains geospatial data for the NLDAS cells comprising Coachella Valley and boundaries for the overall region.

The primary post-processing scripts are ```coachella_valley_forecasts_plots.R``` and ```LOYO_plotting.R```.

  1. ```coachella_valley_forecasts_plots.R``` contains all post-processing and plotting for the complete ensemble models trained on data from 2006-2018 and 2006-2021.
  2. ```LOYO_plotting.R``` contains all post-processing and plotting for the Leave-One-Year-Out ensemble models trained on data from 2006-2021.

Helper functions for the above scripts can be found in ```glmm_functions.R```.

Please email aman.patel1@mssm.edu for any questions about the repository.
