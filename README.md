# Multivariate Bayesian models with flexible shared interactions for analyzing spatio-temporal patterns of rare cancers

This repository contains R code for fitting spatio-temporal models with independent interactions and time-varying shared interactions using INLA. It also includes code to replicate and reproduce the results described in "Multivariate Bayesian models with flexible shared interactions for
analyzing spatio-temporal patterns of rare cancers".

## Table of contents

1.  [Note](#Note)
2.  [R code](#Rcode)
3.  [Acknowledgements](#Acknowledgements)
4.  [References](#Ref)

# Note <a name="Note"/>

The data for both pancreatic cancer and leukaemia incidence and mortality among males in Great Britain during the biennial periods of 2002-2003 and 2018-2019 were analyzed to demonstrate how spatio-temporal models with time-varying shared interactions offer advantages in estimating rare cancers over other spatio-temporal modelling techniques. Precisely, the region under study corresponds to three of the four nations that make up the United Kingdom: England, Wales and Scotland, which comprise the entire island of Great Britain, including small adjacent islands. Although the three nations belong to the same country, the national health system of each operates independently, thus the data have been collected separately and merged into a single database. The URLs from which the cancer and population data for each nation have been obtained are as follows:
1. [NHS England](https://www.cancerdata.nhs.uk/incidence_and_mortality) and [Population England](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/adhocs/12429clinicalcommissioninggroupsccgspopulationestimatesbysingleyearofageandsexenglandmid2001tomid2019)
2. [NHS Wales](https://phw.nhs.wales/services-and-teams/welsh-cancer-intelligence-and-surveillance-unit-wcisu/) and [Population Wales](https://statswales.gov.wales/Catalogue/Population-and-Migration/Population/Estimates/Local-Authority/populationestimates-by-localauthority-year)
3. [NHS Scotland](https://www.opendata.nhs.scot/dataset) and [Population Scotland](https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/population/population-estimates/mid-year-population-estimates/population-estimates-time-series-data)

# R code <a name="Rcode"/>
This folder contains the R code to replicate and reproduce the spatio-temporal multivariate models and to replicate and reproduce the results described in the paper.

The code to replicate and reproduce the spatio-temporal multivariate models proposed in this paper can be found [here](https://github.com/spatialstatisticsupna/Shared_interactions/tree/main/R/All_Models).

The code to replicate and reproduce the results obtained for pancreatic cancer and leukaemia cancer in Great Britain can be found [here](https://github.com/spatialstatisticsupna/Shared_interactions/tree/main/R/Results_Pancreatic) and [here](https://github.com/spatialstatisticsupna/Shared_interactions/tree/main/R/Results_Leukaemia) respectively.



# Acknowledgements <a name="Acknowledgements"/>
This work has been supported by Project PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033 and Ayudas Predoctorales Santander UPNA 2021-2022.
![plot](https://github.com/spatialstatisticsupna/Estimating_LOCP_cancer_mortality_rates/blob/main/micin-aei.jpg)

# References <a name="Ref"/>

Retegui, R., Etxeberria, J.  and Ugarte, M.D (2024). Multivariate Bayesian models with flexible shared interactions for analyzing spatio-temporal
patterns of rare cancers. _Submitted_. 
