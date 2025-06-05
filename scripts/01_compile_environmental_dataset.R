## Title: Canary Islands Plant Trait Data Extraction Project - Naturalis Biodiversity Center----
# This file contains the script that collects environmental data for a list of plant species
# Author: Lucas Jansen (l.s.jansen98@gmail.com)
# Date: 13/02/2025
# Written using R version 4.4.1

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
COORDINATE_UNCERTAINTY = 100
MIN_NUM_OCCURRENCES = 25

# Define Coordinate Reference System (CRS):
COORD_REF_SYSTEM = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Source functions used in this script from '00_aux_functions.R'
source('scripts/00_aux_functions.R')

# LOAD IN FILES: ----
# Load in shapefile
shape_file <- sf::read_sf("data/gadm_land_3/gadm_land_3.shp")

# Load in centroid/institutions reference dataset:
centroid_ref <- utils::read.csv('data/centroid_ref_custom.csv')
institutions_ref <- utils::read.csv('data/institutions_ref_custom.csv')

# Load in island occurrence dataset (Beierkuhnlein et al. 2021)
species_island_data <- utils::read.csv('data/island_occurrence_data_custom.csv')
species_island_data[, names(species_island_data) != "species_name"] <- 
  base::apply(species_island_data[, names(species_island_data) != "species_name"], 2, function(x){x == 1}) # change binary to logical values

# Load in DEM
DEM_rast <- raster::raster("data/DEM/canaryclim_dem.tif")
DEM_rast <- raster::projectRaster(DEM_rast, crs = COORD_REF_SYSTEM)
DEM_rast <- terra::rast(DEM_rast)

# Calculate Terrain Ruggedness Index (TRI) raster from DEM
TRI_rast <- spatialEco::tri(DEM_rast)

# Convert DEM and TRI_rast to RasterLayer
DEM_rast <- raster::raster(DEM_rast)
TRI_rast <- raster::raster(TRI_rast)

# Load in BioClim variable 1 (mean annual temperature):
bioclim1_rast <- raster::raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_01_1979_2013.tif')
bioclim1_rast <- raster::projectRaster(bioclim1_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 2 (mean diurnal temperature range):
bioclim2_rast <- raster::raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_02_1979_2013.tif')
bioclim2_rast <- raster::projectRaster(bioclim2_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 4 (temperature seasonality):
bioclim4_rast <- raster::raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_04_1979_2013.tif')
bioclim4_rast <- raster::projectRaster(bioclim4_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 9 (mean daily temperature of driest quarter):
bioclim9_rast <- raster::raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_09_1979_2013.tif')
bioclim9_rast <- raster::projectRaster(bioclim9_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 11 (mean daily temperature of coldest quarter):
bioclim11_rast <- raster::raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_11_1979_2013.tif')
bioclim11_rast <- raster::projectRaster(bioclim11_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 12 (mean annual precipitation):
bioclim12_rast <- raster::raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_12_1979_2013.tif')
bioclim12_rast <- raster::projectRaster(bioclim12_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 15 (precipitation seasonality):
bioclim15_rast <- raster::raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_15_1979_2013.tif')
bioclim15_rast <- raster::projectRaster(bioclim15_rast, crs = COORD_REF_SYSTEM)

# Load in BioClim variable 17 (mean monthly precipitation of driest quarter):
bioclim17_rast <- raster::raster('data/environmental_data/bioclim/CHELSA_CanaryIslands_bio_17_1979_2013.tif')
bioclim17_rast <- raster::projectRaster(bioclim17_rast, crs = COORD_REF_SYSTEM)

# Load in evapotranspiration data:
evapotranspiration_rast <- raster::raster('data/environmental_data/WaPOR-3/AETI/WaPORV3_average_annualAETI.tif')
evapotranspiration_rast <- raster::projectRaster(evapotranspiration_rast, crs = COORD_REF_SYSTEM)

# Load in evaporation data:
evaporation_rast <- raster::raster('data/environmental_data/WaPOR-3/E/WaPORV3_average_annualE.tif')
evaporation_rast <- raster::projectRaster(evaporation_rast, crs = COORD_REF_SYSTEM)

# Load in transpiration data:
transpiration_rast <- raster::raster('data/environmental_data/WaPOR-3/T/WaPORV3_average_annualT.tif')
transpiration_rast <- raster::projectRaster(transpiration_rast, crs = COORD_REF_SYSTEM)

# Load in interception data:
interception_rast <- raster::raster('data/environmental_data/WaPOR-3/I/WaPORV3_average_annualI.tif')
interception_rast <- raster::projectRaster(interception_rast, crs = COORD_REF_SYSTEM)

# Load in species list:
species_list_extended <- read.csv('data/species_list_extended.csv')

# Load in growth form list:
growth_form_list_extended <- read.csv('data/growth_form_list_extended.csv') 

# Remove missing species (n = 44) from species_list
species_list <- species_list[!species_list$missing, ] # exclude species from list (no occurrence data)

# OCCURRENCE DATA EXTRACTION/FILTERING + ENVIRONMENTAL VARIABLES EXTRACTION:----
## The following step requires a (stable) internet connection and can have a long run-time depending on the number of species
environmental_data <- compile_dataset(species_list = species_list_extended, 
                                      shape_file = shape_file, 
                                      centroid_ref = centroid_ref,
                                      institutions_ref = institutions_ref,
                                      species_island_data = species_island_data,
                                      uncertainty_limit = COORDINATE_UNCERTAINTY,
                                      DEM_rast = DEM_rast,
                                      TRI_rast = TRI_rast,
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
environmental_data_aggregate_100_25 <- environmental_data_100m_extended |> 
  group_by(phylum, order, family, genus, species) |> 
  summarise(num_occurrences = length(species), # Count number of occurrences per species
            elevation_median = median(elevation, na.rm = T), # Elevation median
            elevation_mad = mad(elevation, na.rm = T), # Elevation (MAD)
            TRI_median = median(TRI, na.rm = T), # Terrain Ruggedness Index
            annual_temp = median(annual_temp, na.rm = T), # Annual temperature
            diurnal_temp_range_median = median(diurnal_temp_range, na.rm = T), # Diurnal temperature range
            temp_seasonality = median(temp_seasonality, na.rm = T), # Temperature seasonality
            temp_driest_quarter = median(temp_driest_quarter, na.rm = T), # Temperature driest quarter
            temp_coldest_quarter = median(temp_coldest_quarter, na.rm = T), # Temperature coldest quarter
            annual_precip = median(annual_precip, na.rm = T), # Annual precipitation
            precip_driest_quarter = median(precip_driest_quarter, na.rm = T), # Precipitation driest quarter
            precip_seasonality = median(precip_seasonality, na.rm = T), # Precipitation seasonality
            evapotranspiration = median(evapotranspiration, na.rm = T), # Evapotranspiration
            evaporation = median(evaporation, na.rm = T), # Evaporation
            transpiration = median(transpiration, na.rm = T), # Transpiration
            interception = median(interception, na.rm = T), # Interception
            aridity_index = annual_precip/evapotranspiration) |> # Aridity Index
  filter(num_occurrences >= MIN_NUM_OCCURRENCES) |> # Apply sampling approach
  left_join(growth_form_list_extended[, c('species', 'growth_form')], by = 'species') |> # Add growth form column 
  st_drop_geometry() |> # remove geometry attribute
  dplyr::ungroup() 
  
write.csv(environmental_data_100m_extended, file = 'data/output_data/environmental_dataset_extended_100m.csv')
