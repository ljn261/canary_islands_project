## ------------------------------
## Title: Canary Islands Plant Trait Data Extraction Project - Naturalis Biodiversity Center
# This file contains the script that performs analyses on the compiled environmental dataset
# Author: Lucas Jansen (l.s.jansen98@gmail.com)
# Date: 04/04/2024
# Written using R version 4.3.3

# Libraries:
library(ggcorrplot)
library(MASS)
library(dplyr)
library(tidyr)
library(MVN)
library(ggord)
library(vegan)
library(purrr)
library(FSA)


# CORRELATION ANALYSIS:----
# Calculate Kendall correlation between environmental variables
environmental_data_summary |> 
  dplyr::select(-order, -phylum, -family, -genus, -species, -habitat, -growth_form) |>
  cor(method = 'kendall', use="pairwise.complete.obs") |>
  ggcorrplot()

# A decision was made to keep the following (10) variables: 
  # - elevation (MAD), 
  # - temperature seasonality, 
  # - mean daily air temperature of the coldest quarter, 
  # - annual precipitation, 
  # - precipitation of the driest quarter, 
  # - precipitation seasonality, 
  # - evapotranspiration, 
  # - TRI, 
  # - aridity index, and 
  # - range size. 

# LINEAR DISCRIMINANT ANALYSIS:----
# Check for normality assumption:
environmental_data_summary |>
  dplyr::select(elevation_mad,
                TRI_median,
                temp_seasonality,
                temp_coldest_quarter,
                annual_precip,
                precip_driest_quarter,
                precip_seasonality,
                evapotranspiration,
                aridity_index,
                range_size) |>
  mvn(mvnTest = 'hz', multivariatePlot = 'qq')

# LDA data first needs to be scaled (in order to ensure equal contribution of environmental variables):
lda_data_scaled <- environmental_data_summary |> 
  dplyr::select(elevation_mad, temp_seasonality, temp_coldest_quarter, annual_precip, precip_driest_quarter, 
         precip_seasonality, evapotranspiration, TRI_median, aridity_index, range_size) |>
  scale() |>
  as.data.frame() |>
  mutate(growth_form = environmental_data_summary$growth_form)

# Perform LDA:
lda_model <- lda(growth_form ~ ., data = lda_data_scaled)
lda_model

# Plot 2D LDA model (biplot)
ggord(lda_model, lda_data_scaled$growth_form, repel = FALSE, arrow = 0.2)

# PERMANOVA:----
permanova_data <- environmental_data_summary |>
  dplyr::select(elevation_mad,
                TRI_median,
                temp_seasonality,
                temp_coldest_quarter,
                annual_precip,
                precip_driest_quarter,
                precip_seasonality,
                evapotranspiration,
                aridity_index,
                range_size)

dist_matrix <- vegdist(permanova_data, method = 'euclidian')
permanova_data <- cbind(permanova_data, 'growth_form' = environmental_data_summary$growth_form)
adonis2(dist_matrix ~ growth_form, permanova_data, method = 'euclidian')

# Apply Kruskal-Wallis tests:
environmental_variables <- c('elevation_mad', 
                             'TRI_median', 
                             'temp_seasonality', 
                             'annual_precip', 
                             'temp_coldest_quarter',
                             'precip_driest_quarter',
                             'precip_seasonality',
                             'evapotranspiration', 
                             'aridity_index',
                             'range_size')

lapply(environmental_variables, 
       function(x)kruskal.test(as.formula(paste0(x ," ~ growth_form")), data = permanova_data))

# Apply Dunn's tests:
lapply(environmental_variables, 
       function(x)dunnTest(as.formula(paste0(x ," ~ growth_form")), data = permanova_data, method = 'bonferroni'))


# GROUPED ANALYSIS:----
# Obtain grouped data (Herbaceous versus Woody):
lda_data_scaled_grouped <- lda_data_scaled |> 
  mutate(growth_form = case_when(
    grepl('Herbaceous', growth_form) ~ 'Herbaceous',
    !grepl('Herbaceous', growth_form) ~ 'Woody'))

# Perform grouped LDA:
lda_model_grouped <- lda(growth_form ~ ., data = lda_data_scaled_grouped)
lda_model_grouped

plot(lda_model_grouped)

# PERMANOVA:
permanova_data <- environmental_data_summary |>
  dplyr::select(elevation_mad,
                TRI_median,
                temp_seasonality,
                temp_coldest_quarter,
                annual_precip,
                precip_driest_quarter,
                precip_seasonality,
                evapotranspiration,
                aridity_index,
                range_size)

dist_matrix <- vegdist(permanova_data, method = 'euclidian')

permanova_data_grouped <- permanova_data |>
  mutate(growth_form = environmental_data_summary$growth_form) |> 
  mutate(growth_form = case_when(grepl('Herbaceous', growth_form) ~ 'Herbaceous',
                                 !grepl('Herbaceous', growth_form) ~ 'Woody'))

adonis2(dist_matrix ~ growth_form, permanova_data_grouped, method = 'euclidian')

# Apply Kruskal-Wallis tests:
lapply(environmental_variables, 
       function(x)kruskal.test(as.formula(paste0(x ," ~ growth_form")), data = permanova_data_grouped))

# Apply Dunn's tests:
lapply(environmental_variables, 
       function(x)dunnTest(as.formula(paste0(x ," ~ growth_form")), data = permanova_data_grouped, method = 'bonferroni'))
