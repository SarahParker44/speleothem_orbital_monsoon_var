# speleothem_orbital_monsoon_var

## Code for "A data-model approach to interpreting speleothem oxygen isotope records from monsoon regions on orbital timescales"

### Authors:
Sarah Parker, Sandy P. Harrison, Laia Comas-Bru, Nikita Kaushal, Allegra LeGrande, Martin Werner

### Link to paper:

### Abstract:
Reconstruction of past changes in monsoon climate from speleothem oxygen isotope (δ18O) records is complex because δ18O signals can be influenced by multiple factors including changes in precipitation, precipitation recycling over land, temperature at the moisture source and changes in the moisture source region and transport pathway. Here, we analyse >150 speleothem records from version 2 of the Speleothem Isotopes Synthesis and Analysis (SISAL) database to produce composite regional trends in δ18O in monsoon regions; compositing minimises the influence of site-specific karst and cave processes that can influence individual site records. We compare speleothem δ18O observations with isotope-enabled climate model simulations to investigate the climatic drivers of these regional trends. We focus on differences in δ18O signals between interglacial (mid-Holocene and Last Interglacial) and glacial (Last Glacial Maximum) states, and on δ18O evolution through the Holocene. Differences in speleothem δ18O between the mid-Holocene and Last Interglacial in the East Asian and Indian monsoons are small, despite the larger summer insolation values during the Last Interglacial. Last Glacial Maximum δ18O values are significantly less negative than interglacial values. Comparison with simulated glacial-interglacial δ18O shows that changes are principally driven by global shifts in temperature and precipitation. Holocene speleothem δ18O records show distinct and coherent regional trends. Trends are similar to summer insolation in India, China and southwestern South America, but different in the Indonesian-Australian region. Redundancy analysis shows that 37% of Holocene variability can be accounted for by latitude and longitude, supporting the differentiation of records into individual monsoon regions. Regression analysis of simulated precipitation δ18O and climate variables show that global monsoon Holocene δ18O trends are driven by changes in precipitation, atmospheric circulation and (to a lesser extent) temperature, whilst precipitation recycling is non-significant. However, there are differences in drivers at a regional scale; clear relationships between precipitation and δ18O are seen for India, southwestern South America and the Indonesian-Australian regions, but not the East Asian monsoon. Changes in atmospheric circulation contributes to δ18O trends in the East Asian, Indian and Indonesian-Australian monsoons, and a weak temperature effect is observed over southern and central America and Asia. Recycling is influential in southwestern South America and southern Africa.  Overall, our analyses show that it is possible to differentiate the impacts of specific climatic drivers of precipitation δ18O and use this analysis to interpret orbital-scale changes in speleothem δ18O.

### About this repository:
This repository contains the codes used for the statistical analyses and production of figures for “A data-model approach to interpreting speleothem oxygen isotope records from monsoon regions on orbital timescales”. 

### Data availability: 
#### SISAL 
http://dx.doi.org/10.17864/1947.242
#### ECHAM
https://doi.org/10.1594/PANGAEA.902347
https://doi.pangaea.de/10.1594/PANGAEA.879229
#### GISS
NetCDF files included in "Data" file of this repository. 

## Analyses

### Generation of Figure 1, “sites_map.R”
Generates a world map of sites used in this paper, overlying WOKAM karst data.
Requires world land .shp file from https://www.naturalearthdata.com/downloads/110m-physical-vectors/ and WOKAM karst .shp file from https://produktcenter.bgr.de/terraCatalog/OpenSearch.do?search=473d851c-4694-4050-a37f-ee421170eca8&type=/Query/OpenSearch.do.
Csv files with site latitude and longitude are created within the codes for each analysis. 

### PCoA/RDA analysis (section 2.3), generation of Figure 2, “PCoA_RDA.R”
Principle Coordinate Analysis and Redundancy analysis (using ‘vegan’ R package). Generates a multiplot of PCoA biplot and RDA triplot.
Requires connection to SISALv2 database in MySQL. Requires a function sourced from “d18O_record_length.R”, which calculates the temporal coverage of each site, accounting for gaps, hiatuses and overlaps/gaps between entities. 

### Boxplot analysis (section 2.4) and generation of Figure 3
Calculation of speleothem d18O anomalies from modern: “boxplots_site_mean.R”. Requires connection to SISALv2 database in MySQL.
Preprocessing of ECHAM precipitation, temperature and d18O-precip variables:
Extract precipitation, temperature and d18O-precip, and time, lon and lat dimensions from netCDF files: “Extract_ECHAM_MH_LGM_LIG.m”
Calculate mean annual precipitation and temperature, annual precipitation-weighted mean d18O-precip: “Step1_calcd18O_MH_LGM_LIG.m”
Calculate mean across time periods (30 years long): “Step2_ECHAM_ave_MH_LGM.m”
Calculate anomalies to control run and resample LIG to same grid resolution as MH/LGM: “Step3_calc_anom.m”
Extraction of ECHAM d18Oprecip, temperature and precipitation around each speleothem site used in boxplot analysis, using distance weighted means +/-3 degrees: “Extract_ppt_temp_d18O_boxplots.R”. Uses “ECHAM_dat.mat” file from preprocessing steps. 
Plot results in boxplots, generating figure 3 in “mult_axes_boxplots”. Uses output from “Extract_ppt_temp_d18O_boxplots.R”. 

### Holocene evolution (section 2.5) and generation of Figure 4
Regional speleothem d18O composites for the Holocene are produced, by calculating Z-scores, binning, bootstrap resampling by site (to obtain confidence intervals) and lowess smoothing. Smoothed fit and confidence intervals for each monsoon region are generated and saved to a .csv file, using “region_spel_composites.R”. Requires connection to SISALv2 database in MySQL.
Regional evolution of GISS d18Oprecip extracted using +/-4 degrees around each site (used in speleothem regional composite analysis) distance-weighted means calculated and bootstrap resampling between simulated site d18O trends: “GISS_Hol_trends.R”. 
Speleothem regional composite data (calculated in “region_spel_composites.R”) and simulated GISS d18Oprecip trends (calculated in “GISS_Hol_trends.R”) plotted and exported as pdf in “Holocene_trend_plots.R”. 

### Multiple regression analysis (section 2.6) and generation of Figure 6
Extraction of regional simulated d18Oprecip and climate variables for each Holocene time slice:
	Extract regional d18Oprecip: “extract_Hol_region_d18O.R”
	Extract regional summer precipitation: “extract_Hol_summer_ppt.R”
Extract summer surface air temperature over moisture source areas: “extract_Hol_summer_temp.R”
Extract summer surface winds: “extract_Hol_summer_wdir.R”
Calculate precipitation recycling ratio over regions: “extract_Hol_ppt_recycling.R”
Multiple linear regression analysis carried out in “mult_linear_reg.R”, using variables, saved into csv. Files, from the “Extract_Hol_...R” files. Results from analysis are plotted as partial residual plots and exported as a pdf file. 

### Simulated d18Oprecip versus observed d18Ospel (supplementary figures)
Plotting simulated d18Oprecip with d18Ospel, converted to drip water equivalent, % agreement between model and observations calculated, using "Supp_figs.R", requires connection to SISALv2 database in MySQL.
