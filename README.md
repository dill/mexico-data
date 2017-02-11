---
title: `mexico-data`: line transect survey of pan-tropical spotted dolphins in the Gulf of Mexico
author: David L Miller

---

This repository contains data from NOAA surveys with observations of pan-tropical spotted dolphins (as extracted from OBIS-SEAMAP with additional perpendicular distances to observations provided by NOAA), code to take this data and put it into a format suitable for the R packages [`Distance`](https://github.com/distancedevelopment/Distance) and [`dsm`](https://github.com/distancedevelopment/dsm), and an analysis of these data.

A complete example analysis (and description of the data) is provided here and at [http://distancesampling.org/R/vignettes/mexico-analysis.html](http://distancesampling.org/R/vignettes/mexico-analysis.html).

Note that the compiled HTML here has been run using the latest versions of the `Distance`, `mrds` and `dsm` packages available on github. Ensure you have the latest versions before trying out this analysis.

# Data overview

Data from a combination of several NOAA shipboard surveys conducted on pan-tropical spotted dolphins in the Gulf of Mexico. 47 observations of groups of dolphins. The group size was recorded, as well as the Beaufort sea state at the time of the observation. Coordinates for each observation and bathymetry data were also available as covariates for the analysis.

# Repository structure

The repository is organised as follows:

- `data/` contains the pre and post-processed data
  - `DProject/`: [Distance](http://distancesampling.org/Distance/) original project files that contained the data.
  - `csv/`: the comma-separated value files exported from the Distance project.
  - `Study_ar.*`: shapefiles with the polygon that defines the study area
  - `dolphins.RData`: the data processed and ready for analysis.
- `flatfile/`: data processed into the "flatfile" format required for [`Distance2`](https://github.com/DistanceDevelopment/Distance2) and a script to get it into the right format
- `paperplots/`: code and figures from Miller et al (2013).
- `support/`: more shapefiles and code to get the data into the right format for R
  - `makedata.R`: script to take the raw data and make the `dolphins.RData` file
  - `prediction-grid.*`: prediction grid shapefile(s)
  - `survey-area.*`: survey area shapefile(s)
- `mexico-analysis.Rmd`: example analysis of the data in RMarkdown format
- `build.sh`: bash script to render the HTML report
- `mexico-analysis.html`: example analysis of the data rendered to HTML


# References

Halpin, P.N., A.J. Read, E. Fujioka, B.D. Best, B. Donnelly, L.J. Hazen, C. Kot, K. Urian, E. LaBrecque, A. Dimatteo, J. Cleary, C. Good, L.B. Crowder, and K.D. Hyrenbach. 2009. OBIS-SEAMAP: The world data center for marine mammal, sea bird, and sea turtle distributions. *Oceanography* **22**(2):104-115

Miller, D.L., Burt, M.L., Rexstad, E.A., & Thomas, L. (2013). Spatial models for distance sampling data: recent developments and future directions. Methods in Ecology and Evolution, 4(11), 1001–1010. [http://doi.org/10.1111/2041-210X.12105](http://doi.org/10.1111/2041-210X.12105).

NOAA Southeast Fisheries Science Center. 1996. Report of a Cetacean Survey of Oceanic and Selected Continental Shelf Waters of the Northern Gulf of Mexico aboard NOAA Ship Oregon II (Cruise 220)


