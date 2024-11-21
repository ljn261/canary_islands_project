## ------------------------------
## Title: Canary Islands Plant Trait Data Extraction Project - Naturalis Biodiversity Center
# This file contains the script that collects environmental data for a list of plant species
# Author: Lucas Jansen (l.s.jansen98@gmail.com)
# Date: 04/04/2024
# Written using R version 4.3.3

# Preliminaries:----
# Libraries:
library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(stringr)
library(rgbif)
library(readxl)
library(sf)
library(raster)
library(terra)

# Define coordinate uncertainty and minimum number of occurrences
COORDINATE_UNCERTAINTY = 1000
MIN_NUM_OCCURRENCES = 25

# Define Coordinate Reference System (CRS):
COORD_REF_SYSTEM = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Source functions used in this script from '00_aux_functions.R'
source('scripts/00_aux_functions.R')

# LOAD IN FILES: ----
# Load in shapefile
shape_file <- read_sf("data/gadm_land_3/gadm_land_3.shp")

# Load in centroid/institutions reference dataset:
centroid_ref <- read.csv('data/centroid_ref_custom.csv')
institutions_ref <- read.csv('data/institutions_ref_custom.csv')

# Load in island occurrence dataset (Beierkuhnlein et al. 2021)
species_island_data <- read.csv('data/island_occurrence_data_custom.csv')
species_island_data[, names(species_island_data) != "species_name"] <- 
  apply(species_island_data[, names(species_island_data) != "species_name"], 2, function(x){x == 1}) # change binary to logical values

# Load in DEM
DEM_rast <- raster("data/DEM/canaryclim_dem.tif")
DEM_rast <- terra::rast(DEM_rast)

# Calculate Terrain Ruggedness Index (TRI) raster from DEM
TRI_rast <- spatialEco::tri(DEM_rast)

# Convert DEM and TRI_rast to RasterLayer
DEM_rast <- raster(DEM_rast)
TRI_rast <- raster(TRI_rast)

# Load in habitat map:
IUCN_hab_map <- raster('data/environmental_data/IUCN/iucn_habitatclassification_composite_lvl2_ver003.tif')
IUCN_hab_map <- projectRaster(IUCN_hab_map, crs = COORD_REF_SYSTEM) # Change CRS, takes a long time!

# Load in BioClim variable 1 (mean annual temperature):
bioclim1_rast <- raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_01_1979_2013.tif')
bioclim1_rast <- projectRaster(bioclim1_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 2 (mean diurnal temperature range):
bioclim2_rast <- raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_02_1979_2013.tif')
bioclim2_rast <- projectRaster(bioclim2_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 4 (temperature seasonality):
bioclim4_rast <- raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_04_1979_2013.tif')
bioclim4_rast <- projectRaster(bioclim4_rast, crs = COORD_REF_SYSTEM)
#bioclim4_rast <- bioclim4_rast/100

# Load in BioClim variable 9 (mean daily temperature of driest quarter):
bioclim9_rast <- raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_09_1979_2013.tif')
bioclim9_rast <- projectRaster(bioclim9_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 11 (mean daily temperature of coldest quarter):
bioclim11_rast <- raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_11_1979_2013.tif')
bioclim11_rast <- projectRaster(bioclim11_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 12 (mean annual precipitation):
bioclim12_rast <- raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_12_1979_2013.tif')
bioclim12_rast <- projectRaster(bioclim12_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 15 (precipitation seasonality):
bioclim15_rast <- raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_15_1979_2013.tif')
bioclim15_rast <- projectRaster(bioclim15_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 17 (mean monthly precipitation of driest quarter):
bioclim17_rast <- raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_17_1979_2013.tif')
bioclim17_rast <- projectRaster(bioclim17_rast, crs = COORD_REF_SYSTEM)

# Load in evapotranspiration data:
evapotranspiration_rast <- raster('data/environmental_data/WaPOR-3/AETI/WaPORV3_average_annualAETI.tif')
evapotranspiration_rast <- projectRaster(evapotranspiration_rast, crs = COORD_REF_SYSTEM)

# Load in evaporation data:
evaporation_rast <- raster('data/environmental_data/WaPOR-3/E/WaPORV3_average_annualE.tif')
evaporation_rast <- projectRaster(evaporation_rast, crs = COORD_REF_SYSTEM)

# Load in transpiration data:
transpiration_rast <- raster('data/environmental_data/WaPOR-3/T/WaPORV3_average_annualT.tif')
transpiration_rast <- projectRaster(transpiration_rast, crs = COORD_REF_SYSTEM)

# Load in interception data:
interception_rast <- raster('data/environmental_data/WaPOR-3/I/WaPORV3_average_annualI.tif')
interception_rast <- projectRaster(interception_rast, crs = COORD_REF_SYSTEM)

# Load in species list:
species_list <- read.csv('data/species_list.csv')

# Load in growth form list:
growth_form_list <- read.csv('data/growth_form_list.csv') 

# Remove missing species (n = 44) from species_list
species_list <- species_list[!species_list$missing, ] # exclude species from list (no occurrence data)

# OCCURRENCE DATA EXTRACTION/FILTERING + ENVIRONMENTAL VARIABLES EXTRACTION:----
environmental_data <- compile_dataset(species_list = species_list, 
                                      shape_file = shape_file, 
                                      centroid_ref = centroid_ref,
                                      institutions_ref = institutions_ref,
                                      species_island_data = species_island_data,
                                      uncertainty_limit = COORDINATE_UNCERTAINTY,
                                      DEM_rast = DEM_rast,
                                      TRI_rast = TRI_rast,
                                      habitat_rast = IUCN_hab_map,
                                      bioclim1_rast = bioclim1_rast,
                                      bioclim2_rast = bioclim2_rast,
                                      bioclim4_rast = bioclim4_rast,
                                      bioclim9_rast = bioclim9_rast,
                                      bioclim11_rast = bioclim11_rast,
                                      bioclim12_rast = bioclim12_rast,
                                      bioclim15_rast = bioclim15_rast,
                                      bioclim17_rast = bioclim17_rast,
                                      evapotranspiration_rast = evapotranspiration_rast,
                                      evaporation_rast = evaporation_rast,
                                      transpiration_rast = transpiration_rast,
                                      interception_rast = interception_rast)


 # DATASET COMPILATION (AGGREGATION TO SPECIES-LEVEL):----   
environmental_data_summary_500_13 <- environmental_data_500m |> 
  group_by(phylum, order, family, genus, species) |> 
  summarise(num_occurrences = length(species),
            elevation_median = median(elevation, na.rm = T), # ELEVATION---
            elevation_mad = mad(elevation, na.rm = T),
            TRI_median = median(TRI, na.rm = T), # TRI---
            annual_temp = median(annual_temp, na.rm = T), # TEMPERATURE---
            diurnal_temp_range_median = median(diurnal_temp_range, na.rm = T),
            temp_seasonality = median(temp_seasonality, na.rm = T),
            temp_driest_quarter = median(temp_driest_quarter, na.rm = T),
            temp_coldest_quarter = median(temp_coldest_quarter, na.rm = T),
            annual_precip = median(annual_precip, na.rm = T), # PRECIPITATION---
            precip_driest_quarter = median(precip_driest_quarter, na.rm = T),
            precip_seasonality = median(precip_seasonality, na.rm = T),
            habitat = toString(unique(na.omit(habitat_information))), # turn habitat information to string of unique habitat codes
            habitat_num = length(unique(na.omit(habitat_information))),
            evapotranspiration = median(evapotranspiration, na.rm = T), # EVAPOTRANSPIRATION---
            evaporation = median(evaporation, na.rm = T),
            transpiration = median(transpiration, na.rm = T),
            interception = median(interception, na.rm = T),
            aridity_index = annual_precip/evapotranspiration) |>
  filter(num_occurrences >= MIN_NUM_OCCURRENCES) |> 
  mutate(range_size = geometry |>
           st_combine() |>
           st_concave_hull(0.1) |>
           st_area() |>
           as.numeric() / 1000000 |>
           round(2)) |>
  # mutate(#habitat = lapply(strsplit(habitat, ",", TRUE), as.numeric),
  #        range_size = round( # Calculate range size: ---
  #          as.numeric(st_area( # calculate area of convex hull
  #              st_convex_hull( # create convex hull around points
  #                st_combine(geometry)
  #              )
  #            )
  #          )/1000000, # convert range size from m^2 to km^2
  #          2)) |> # round to 2 decimals
  left_join(growth_form_list[, c('species', 'growth_form')], by = 'species') |> # Add growth form column 
  st_drop_geometry() |> # remove geometry attribute
  dplyr::ungroup() 
  
write.csv(environmental_data_1000m, file = 'data/output_data/environmental_dataset_full_1000m.csv')
