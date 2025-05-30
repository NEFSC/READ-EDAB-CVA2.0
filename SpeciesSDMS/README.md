# SpatialSDMS - README

This folder contains all of the codes used to generate the ensemble species distribution models (SDMs). As of May 2025, there is a single script for each species, but this may be modified to accomodate running multiple species on a container or in parallel.

The ensemble SDMs currently includes the following models:

1.  MAXENT
2.  GAM
3.  sdmTMB (GLMM)
4.  Random Forest (RF)

MAXENT and GAM models, as well as the final model ensemble, are generated using the [EFHSDM package.](https://github.com/alaska-groundfish-efh/EFHSDM) Random Forest models are generated using Random Forest Spatial Interpolation (RFSI) from the [mateo](https://github.com/AleksandarSekulic/Rmeteo) package. The GLMM is built using [sdmTMB](https://pbs-assess.github.io/sdmTMB/).

These scripts are currently written to run locally. The MAXENT, GAM, and RF models run quickly (within an hour total). sdmTMB, with the current mesh resolution, takes about 3 hours to run locally. Decreasing the resolution of the mesh decreases run times, with small changes in root mean squared error. Ensemble model creation may be moved to a container in the future, which may result in small changes in workflows.

This project is under active development.

# Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
