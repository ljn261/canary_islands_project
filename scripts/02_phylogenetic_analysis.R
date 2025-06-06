## Title: Canary Islands Plant Trait Data Extraction Project - Naturalis Biodiversity Center----
# This file contains the script that collects information from the Phylogenetic ANOVA
# Author: Lucas Jansen (l.s.jansen98@gmail.com)
# Date: 05/06/2025
# Written using R version 4.4.1

# Source functions used in this script from '00_aux_functions.R'
source('../scripts/00_aux_functions.R')

# Environmental dataset of balanced sampling approach:
environmental_dataset <- read.csv("./data/output_data/environmental_dataset_balanced_aggregate.csv")

# 1. Correlation analysis----
# Compute and visualise Kendall's correlation per variable
environmental_dataset |> 
  dplyr::select(
    elevation_mad,
    elevation_median,
    TRI_median,
    annual_temp,
    diurnal_temp_range_median,
    temp_seasonality,
    temp_driest_quarter,
    temp_coldest_quarter,
    annual_precip,
    precip_driest_quarter,
    precip_seasonality,
    evapotranspiration,
    evaporation,
    transpiration,
    interception,
    aridity_index) |>
  stats::cor(method = 'kendall', use="pairwise.complete.obs") |>
  ggcorrplot::ggcorrplot(lab = TRUE)

# Some variables are very strongly correlated ($\tau > 0.7$): 
#   - elevation (median) and annual temperature/temperature of the coldest quarter 
# - annual temperature and temperature of the coldest quarter 
# - diurnal temperature range and temperature seasonality 
# - precipitation of the driest quarter and evapotranspiration 
# - precipitation of the driest quarter and transpiration 
# - evapotranspiration and transpiration/interception 
# - transpiration and interception
# 
# A selection of the following variables has been made: 
#   - elevation (MAD), 
# - temperature seasonality, 
# - mean daily air temperature of the coldest quarter, 
# - annual precipitation, 
# - precipitation seasonality, 
# - evaporation, 
# - interception, 
# - TRI, and 
# - aridity index

environmental_dataset |> 
  dplyr::select(
    elevation_mad,
    TRI_median,
    temp_seasonality,
    temp_coldest_quarter,
    annual_precip,
    precip_seasonality,
    evaporation,
    interception,
    aridity_index) |>
  stats::cor(method = 'kendall', use="pairwise.complete.obs") |>
  ggcorrplot::ggcorrplot(lab = TRUE)

# 2. Create phylogeny using V.PhyloMaker2----
phylo_tree_phylomaker_input <- environmental_dataset |> 
  dplyr::select(species, genus, family) |> 
  cbind(data.frame(species.relative = "", genus.relative = ""))

phylo_tree_phylomaker_relatives <- V.PhyloMaker2::bind.relative(sp.list = phylo_tree_phylomaker_input, 
                                                                tree = V.PhyloMaker2::GBOTB.extended.TPL, 
                                                                nodes = V.PhyloMaker2::nodes.info.1.TPL)

phylo_tree_phylomaker <- V.PhyloMaker2::phylo.maker(sp.list = phylo_tree_phylomaker_input, 
                                                    tree = phylo_tree_phylomaker_relatives$phylo, 
                                                    nodes = phylo_tree_phylomaker_relatives$nodes.info, 
                                                    scenarios = 'S3')

phylo_tree_phylomaker$scenario.3$tip.label <- gsub("_", " ", phylo_tree_phylomaker$scenario.3$tip.label)

# 3. Run Phylogenetic PCA----
phylo_pca <- run_phylo_pca(environmental_data = environmental_dataset, 
                                  phylo_tree = phylo_tree_phylomaker$scenario.3)

# Plot in 3D space using plotly
phylo_pca_plot <- plot_3d_phylo_pca(environmental_data = environmental_dataset, 
                                    phylo_tree = phylo_tree_phylomaker$scenario.3)
phylo_pca_plot

# 4. Run Phylogenetic ANOVA----
phytools::phylANOVA(tree = phylo_tree_phylomaker$scenario.3, 
                    x = environmental_dataset$growth_form, 
                    y = environmental_dataset$TRI_median, # Adjust environmental variable
                    nsim = 1000,
                    posthoc = TRUE,
                    p.adj = "holm")
