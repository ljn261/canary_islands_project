##-------------------------------
## Title: Canary Islands Plant Trait Data Extraction Project - Naturalis Biodiversity Center
# This file contains a collection of functions used to extract plant trait data
# Author: Lucas Jansen (l.s.jansen98@gmail.com)
# Date: 04/04/2024
# Written using R version 4.3.3


extract_locality_data <- function(binomial_name){
  # Load data from GBIF:
  message(sprintf('\nExtracting GBIF locality data for %s...', binomial_name))
  occurrence_data <- rgbif::occ_search(scientificName = binomial_name,                         # requires binomial name w/o author
                                       limit = 5000,                                           # limit = max. no. of occurrences
                                       hasCoordinate = T,                                      # hasCoordinate = T only downloads occurrences with coordinates
                                       decimalLongitude = "-18.341578, -13.265994",            # longitudinal boundary (roughly drawing a box around the CI)
                                       decimalLatitude = "27.508589, 29.518117")                          
  
  
  occurrence_data_relevant <- occurrence_data$data # select 'data' element from list
  
  if(!("coordinateUncertaintyInMeters" %in% colnames(occurrence_data_relevant))){ # check if coordinate uncertainty column is missing
    
    coordinateUncertaintyInMeters <- rep(0, nrow(occurrence_data_relevant)) # add vector of 0s as coordinateUncertainty
    
    occurrence_data_relevant <- cbind(occurrence_data_relevant, coordinateUncertaintyInMeters) 
  }
  
  occurrence_data_relevant <- occurrence_data_relevant %>% # select relevant columns
    dplyr::select(gbifID,
                  datasetKey,
                  phylum,
                  order,
                  family,
                  genus,
                  species,
                  taxonRank,
                  countryCode,
                  coordinateUncertaintyInMeters,
                  #locality,
                  #year, 
                  basisOfRecord, 
                  #institutionCode, 
                  #datasetName, 
                  issues,
                  identifier,
                  decimalLongitude, 
                  decimalLatitude)
  
  # Convert country code from ISO2c to ISO3c, the latter is required for CoordinateCleaner package
  occurrence_data_relevant$countryCode <- countrycode::countrycode(occurrence_data_relevant$countryCode, 
                                                                   origin =  'iso2c', 
                                                                   destination = 'iso3c') 
  
  # Convert dataframe to data.frame object 
  occurrence_data_relevant <- data.frame(occurrence_data_relevant)
  
  return(occurrence_data_relevant)
}

flag_locality_data <- function(locality_data, shape_file, centroid_ref, institutions_ref, species_island_data, uncertainty_limit = 5000){
  flags <- CoordinateCleaner::clean_coordinates(x = locality_data,
                                                verbose = T,
                                                lon = "decimalLongitude",
                                                lat = "decimalLatitude",
                                                countries = "countryCode",
                                                seas_ref = shape_file,
                                                country_ref = shape_file,     
                                                centroids_ref = centroid_ref,            # centroid_ref file contains coordinates of centroids of islands
                                                centroids_rad = 500,                      # define radius around centroid in which occurrences will be flagged
                                                inst_ref = institutions_ref,
                                                inst_rad = 500,
                                                outliers_method = 'quantile',
                                                outliers_mtp = 15,
                                                outliers_td = 1000,
                                                species = "species",
                                                tests = c("centroids", 
                                                          "equal",
                                                          "gbif", 
                                                          "institutions",
                                                          "duplicates",
                                                          "zeros", 
                                                          "countries", 
                                                          "seas", 
                                                          "capitals"))
  
  if(dim(species_island_data)[1] == 1){
    island_flag <- flag_island_occurrence(locality_data, shape_file, species_island_data) # for each entry, check if occurs on correct island (according to Beierkuhnlein)
    flags <- flags %>% mutate('island' = island_flag) # add column containing island flags to flags dataframe
  } else {
    stop('Incorrect dimension of species island occurrence data. \n')
  }
  
  cdrep_flag <- flag_locality_issue(locality_data) # for each entry, check whether 'cdrep' issue is present
  flags <- flags %>% mutate('cdrep_issue' = cdrep_flag)
  
  uncertainty_flag <- flag_coordinate_uncertainty(locality_data, uncertainty_limit)
  flags <- flags %>% mutate('coord_uncertainty' = uncertainty_flag)
  
  return(flags)
}

flag_island_occurrence <- function(locality_data, shape_file, species_island_data){
  
  message('Testing island occurrence----')
  
  shape_vector <- terra::vect(shape_file) # convert shapefile to SpatVector
  
  locality_coordinates <- locality_data[,c("decimalLongitude", "decimalLatitude")] # subset coordinate columns
  
  locality_vector <- terra::vect(locality_data, geom = c("decimalLongitude", "decimalLatitude"), crs = shape_vector) # convert coordinates to SpatVector 
  
  coordinate_island_data = terra::extract(shape_vector, locality_coordinates)[['island']] # extracting island data based on coordinates
  
  occurs_on_island = sapply(coordinate_island_data, check_island_occurrence, species_island_data = species_island_data) # return a logical vector containing, for each occurrence, whether it is in line with the Beierkuhnlein data
  
  message(sprintf("Flagged %s records for island occurrence", 
                  sum(!occurs_on_island, na.rm = TRUE)))
  
  return(unname(unlist(occurs_on_island)))
}

check_island_occurrence <- function(island_name, species_island_data){
  if (is.na(island_name)){
    return(FALSE)
  }
  if (island_name == 'Tenerife'){
    return(species_island_data[['T']])
  } 
  if (island_name == 'Gran Canaria'){
    return(species_island_data[['C']])
  }
  if (island_name == 'Fuerteventura'){
    return(species_island_data[['F']])
  } 
  if (island_name == 'Lanzarote'){
    return(species_island_data[['L']])
  } 
  if (island_name == 'La Gomera'){
    return(species_island_data[['G']])
  } 
  if (island_name == 'La Palma'){
    return(species_island_data[['P']])
  } 
  if (island_name == 'El Hierro'){
    return(species_island_data[['H']])
  }
}

flag_locality_issue <- function(locality_data){
  message('Testing cdrep issue----')
  
  issue_data <- locality_data$issues # subset issue data
  
  cdrep_issue <- sapply(issue_data, check_locality_issue) # make logical vector whether cdrep issue is found for instance
  
  message(sprintf("Flagged %s records for cdrep issue", 
                  sum(!cdrep_issue, na.rm = TRUE)))
  
  return(cdrep_issue)
}

check_locality_issue <- function(issue_data){
  issues <- unlist(strsplit(issue_data, ','))
  
  cdrep_issue <- !('cdrep' %in% issues)
  
  return(cdrep_issue) # return TRUE if 'cdrep' is *not* in issues
}

flag_coordinate_uncertainty <- function(locality_data, uncertainty_limit){
  
  message('Testing coordinate uncertainty----')
  coord_uncertainty <- locality_data$coordinateUncertaintyInMeters
  
  coord_uncertainty <- coord_uncertainty <= uncertainty_limit # check which coordinates are below uncertainty_limit
  
  coord_uncertainty[is.na(coord_uncertainty)] = TRUE # replace NA values with FALSE
  
  message(sprintf("Flagged %s records for coordinate uncertainty", 
                  sum(!coord_uncertainty, na.rm = TRUE)))
  
  return(coord_uncertainty)
}

compile_dataset <- function(species_list, 
                            shape_file, 
                            centroid_ref, 
                            institutions_ref, 
                            species_island_data, 
                            uncertainty_limit, 
                            DEM_rast, 
                            TRI_rast,
                            habitat_rast, 
                            bioclim1_rast,
                            bioclim2_rast,
                            bioclim4_rast,
                            bioclim5_rast,
                            bioclim6_rast,
                            bioclim9_rast,
                            bioclim10_rast,
                            bioclim11_rast,
                            bioclim12_rast,
                            bioclim14_rast,
                            bioclim15_rast,
                            bioclim17_rast,
                            evapotranspiration_rast,
                            evaporation_rast,
                            transpiration_rast,
                            interception_rast){
  
  start_time = Sys.time() # track time of compilation
  pb <- txtProgressBar(min = 0, max = length(species_list$species_name), style = 3)
  
  report_columns = c('family', 
                     'genus', 
                     'species_original',
                     'species_gbif',
                     'clean_none', 
                     'clean_summary', 
                     'clean_island',
                     'clean_cdrep',
                     'clean_uncertainty') # columns for the report dataframe
  
  cleaning_report = data.frame(matrix(nrow = 0, ncol = length(report_columns))) # assigning report_df as empty dataframe (will be filled row-wise)
  
  environmental_data = data.frame(matrix(nrow = 0, ncol = 0)) # assign empty dataframe to fill with environmental data
  
  
  for(i in seq_along(species_list$species_name)){ # loop through all species
    species_name <- species_list$species_name[i]
    
    setTxtProgressBar(pb, i) # adding to progress bar
    
    base::tryCatch(
      {
        # Locality data extraction:
        locality_data <- extract_locality_data(species_name) # extract locality information (coordinates)
        
        # Locality data flagging:
        species_island_data_subset <- species_island_data[species_island_data$species_name == species_name, ] # subset species occurrence data on island
        
        flags <- flag_locality_data(locality_data, 
                                    shape_file, 
                                    centroid_ref, 
                                    institutions_ref, 
                                    species_island_data_subset, 
                                    uncertainty_limit = uncertainty_limit) # flag locality data
        
        locality_data <- locality_data %>% 
          mutate(.summary = flags$.summary) %>%
          mutate(island = flags$island) %>%
          mutate(cdrep = flags$cdrep_issue) %>%
          mutate(coord_uncertainty = flags$coord_uncertainty)
        
        # Locality data cleaning:
        locality_data_clean_1 <- locality_data[locality_data$.summary, ] # cleaning step 1
        locality_data_clean_2 <- locality_data_clean_1[locality_data_clean_1$island, ] # cleaning step 2
        locality_data_clean_3 <- locality_data_clean_2[locality_data_clean_2$cdrep, ] # cleaning step 3
        locality_data_clean_4 <- locality_data_clean_3[locality_data_clean_3$coord_uncertainty, ] # cleaning step 4
        
        
        family_name = unique(locality_data$family) # extract family name
        genus_name = unique(locality_data$genus) # extract genus name
        species_name_gbif = unique(locality_data$species)
        
        clean_none = dim(locality_data)[1]
        clean_summary = dim(locality_data_clean_1)[1]
        clean_island = dim(locality_data_clean_2)[1]
        clean_cdrep = dim(locality_data_clean_3)[1]
        clean_uncertainty = dim(locality_data_clean_4)[1]
        
        report_vector <- c(family_name, 
                           genus_name, 
                           species_name, 
                           species_name_gbif,
                           clean_none, 
                           clean_summary, 
                           clean_island,
                           clean_cdrep,
                           clean_uncertainty)
        
        cleaning_report = rbind(cleaning_report, report_vector)
        
        
        if(!(dim(locality_data_clean_4)[1] == 0)){ # extract variables ONLY if there are occurrences after cleaning!
          cleaned_coordinates <- st_as_sf(x = locality_data_clean_4, 
                                          coords = c("decimalLongitude", "decimalLatitude"),
                                          crs = COORD_REF_SYSTEM) # convert cleaned coordinates to spatial points
          
          # Environmental variable extraction:
          elevation <- raster::extract(DEM_rast, cleaned_coordinates) # extract ELEVATION
          cleaned_coordinates <- cleaned_coordinates %>% mutate(elevation = elevation) # add elevation data to dataframe
          
          TRI <- raster::extract(TRI_rast, cleaned_coordinates) # extract TERRAIN RUGGEDNESS INDEX
          cleaned_coordinates <- cleaned_coordinates %>% mutate(TRI = TRI) # add TRI data to dataframe
          
          habitat_information <- raster::extract(habitat_rast, cleaned_coordinates) # extract HABITAT
          cleaned_coordinates <- cleaned_coordinates %>% mutate(habitat_information = habitat_information) # add habitat information to dataframe
          
          bioclim1 <- raster::extract(bioclim1_rast, cleaned_coordinates) # exctract annual temperature
          cleaned_coordinates <- cleaned_coordinates %>% mutate(annual_temp = bioclim1)
          
          bioclim2 <- raster::extract(bioclim2_rast, cleaned_coordinates) # extract diurnal temperature range
          cleaned_coordinates <- cleaned_coordinates %>% mutate(diurnal_temp_range = bioclim2)
          
          bioclim4 <- raster::extract(bioclim4_rast, cleaned_coordinates) # extract temperature seasonality
          cleaned_coordinates <- cleaned_coordinates %>% mutate(temp_seasonality = bioclim4)
          
          bioclim9 <- raster::extract(bioclim9_rast, cleaned_coordinates) # extract temperature of driest quarter
          cleaned_coordinates <- cleaned_coordinates %>% mutate(temp_driest_quarter = bioclim9)
          
          bioclim11 <- raster::extract(bioclim11_rast, cleaned_coordinates) # extract temperature of coldest quarter
          cleaned_coordinates <- cleaned_coordinates %>% mutate(temp_coldest_quarter = bioclim11)
          
          bioclim12 <- raster::extract(bioclim12_rast, cleaned_coordinates) # extract annual precipitation
          cleaned_coordinates <- cleaned_coordinates %>% mutate(annual_precip = bioclim12)
          
          bioclim15 <- raster::extract(bioclim15_rast, cleaned_coordinates) # extract precipitation seasonality
          cleaned_coordinates <- cleaned_coordinates %>% mutate(precip_seasonality = bioclim15)
          
          bioclim17 <- raster::extract(bioclim17_rast, cleaned_coordinates) # extract precipitation of driest quarter
          cleaned_coordinates <- cleaned_coordinates %>% mutate(precip_driest_quarter = bioclim17)
          
          evapotranspiration <- raster::extract(evapotranspiration_rast, cleaned_coordinates) # extract evapotranspiration
          cleaned_coordinates <- cleaned_coordinates %>% mutate(evapotranspiration = evapotranspiration)
          
          evaporation <- raster::extract(evaporation_rast, cleaned_coordinates) # extract evaporation
          cleaned_coordinates <- cleaned_coordinates %>% mutate(evaporation = evaporation)
          
          transpiration <- raster::extract(transpiration_rast, cleaned_coordinates) # extract transpiration
          cleaned_coordinates <- cleaned_coordinates %>% mutate(transpiration = transpiration)
          
          interception <- raster::extract(interception_rast, cleaned_coordinates) # extract interception
          cleaned_coordinates <- cleaned_coordinates %>% mutate(interception = interception)
          
          environmental_data <- rbind(environmental_data, cleaned_coordinates) # add information to final dataframe
        }
      },
      error = function(e) {
        cat(paste("Error encountered for species:", species_name, " \nError message:", e$message, " \n"))
        return(NULL)
      }
    ) # end of TryCatch
  }
  
  close(pb) # closing connection to progress bar
  
  colnames(cleaning_report) = report_columns # change column names of report_df
  class(cleaning_report$clean_none) <- 'integer'
  class(cleaning_report$clean_summary) <- 'integer'
  class(cleaning_report$clean_island) <- 'integer'
  class(cleaning_report$clean_cdrep) <- 'integer'
  class(cleaning_report$clean_uncertainty) <- 'integer'
  
  cleaning_report <<- cleaning_report
  
  end_time = Sys.time()
  
  diff_time = difftime(end_time, start_time, units = 'mins')
  
  message(sprintf("DATASET COMPILATION took %s minute(s)-----", round(diff_time, 2)))
  
  return(environmental_data)
}


