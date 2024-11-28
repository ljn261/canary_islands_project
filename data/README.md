This folder contains the (raw) data used in our study

## Contents
- `./DEM/`: folder containing the Digital Elevation Model (DEM) of the Canary Islands from the [CanaryClim (topographical) dataset](https://figshare.com/articles/dataset/CanaryClim_-_Topographic_variables/22060433)
- `./environmental_data/`: folder containing the environmental raster layers used in this study to extract environmental variables, including the IUCN habitat classification map, FAO WAPOR-3 evapotranspiration data and CanaryClim data.
- `./gadm_land/`: folder containing shape file of Canary Islands defining the coordinates of its land masses, as extracted from the [Database of Global Administrative Areas (GADM)](https://gadm.org/)
- `./output_data/`: folder where environmental dataset (.csv) will be written
- `./centroid_ref_custom.csv`: .csv file containing a custom list of island centroids used in filtering species occurrences
- `./growth_form_list.csv`: .csv file containing the list of species and associated growth forms
- `./institutions_ref_custom.csv`: .csv file containing a custom list of (botanical) institutions used in filtering species occurrences
- `./island_occurrence_data_custom.csv`: .csv file containing a (revised) Flora of Canary Islands, based on the publication of [Beierkuhnlein et al. (2021)](https://www.mdpi.com/1424-2818/13/10/480) indicating the presence (1) or absence (0) of a species on specific islands
- `./species_list.csv`: .csv file containing the list of species analysed in the research and whether they are missing from the GBIF repository
