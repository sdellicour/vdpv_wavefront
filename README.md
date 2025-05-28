This repo gathers the R scripts related to the wavefront analyses performed in the study entitled "**Unravelling spatiotemporal heterogeneities of wild and vaccine-derived poliovirus spread: past and present**" (Candido *et al*., *submitted*), which are all gathered within the file `R_script_VDPV.r`.

Abstract: Outbreaks of vaccine-derived poliovirus type 2 (cVDPV2) have become a major threat to polio eradication. We analyse cVDPV2 cases and Wild-type poliovirus 1 (WPV1) sequences to uncover the spatiotemporal patterns and drivers of poliovirus spread. Between May 1, 2016 and September 29, 2023, 3120 cVDPV2 cases were reported across 76 outbreaks and 39 countries globally. Outbreaks have mostly been small (median = 5 cases, range 1-578), have spread to a median maximum distance of 231 km (0-4442) and for median duration of 202 days (0-1905). Wavefront velocity analysis of large outbreaks reveals a median velocity of spread of 2.3 km/day (1.0-4.4). International borders are associated with a slower velocity of spread (p < 0.001), when in the presence of high population immunity. Finally, phylogeographic analysis of 1572 global sequences, including 38 newly generated (1953-2015), reveals that historic WPV1 spread resembles recent cVDPV2 patterns and that international spread is largely sustained by unidirectional movement between neighbouring countries. Our findings offer insights for enhancing the geographical scope of poliovirus surveillance and vaccination response, crucial steps in the final phases of poliovirus eradication 

## Simulated data and other files

1. Due to data sharing restrictions, the data made available here is a simulation of an outbreak with characteristics similar to that of a poliovirus outbreak. The simulated dataset was generated with the "simulatorRRW2" available in the R package "[seraphim](https://github.com/sdellicour/seraphim)".

2. The shapefiles available in this repository were retrieved from [Natural Earth](https://www.naturalearthdata.com/).

3. Some raster files were too large to be deposited in this GitHub repo. This is the case for the rasters used for testing the impact of covariates on the wavefront velocity, but these files can be found as follows:

- the "friction" raster (friction to human movement) can be retrieved from [Weiss et al.](https://www.nature.com/articles/nature25181) (2018, *Nature*), moved to the main folder, and renamed as "Friction_2015.tif" before running the script.

- the population density raster can be downloaded from the [WorldPop](https://hub.worldpop.org/geodata/summary?id=24777) website, moved to the main folder, and renamed as "Human_popD.tif" before running the script.
