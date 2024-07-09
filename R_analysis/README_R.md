# Lara Calvo--BMSC Fall Program 2020
## Directed Studies in Marine Science - R-script
## Advisor: Jacob Usinowicz

## R Analyses
This folder contains the main R code for downloading data, importing it into R, fitting SDMs with either MaxEnt or GAMs, then applying competition models between different populations of Daphnia.

The main output is in 3 forms: 
1. Images of maps or other summary plots, found in the output folder.
2. Rasters (.tif) that are saved in the QGIS folder to produce prettier maps.
3. An sqlite database (daphniabc.sqlite) that can be written to with the rasters and then read from (either by downloading or accessing remotely)

This folder contains the files: 
**Final_2021_edits.R**: Laura's final(ish) working code incorporating all of the analyses and generating plots for final papers and presentations. 

**presence_absence_SDM_popdyn.R**: Jacob's final(ish) working code incorporating all of the analyses, generating plots, generating .TIFs, and writing to the SQLite database. 

**old_code**: This folder contains a mix of code drafts, educational code documents, and works-in-progress to continue to build on this work. 

**cmip5, wc2-t**: database files which have been downloaded from WorldClim