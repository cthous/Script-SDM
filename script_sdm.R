
#### PACKAGES ####
# Core spatial and plotting packages
library(sf)
library(ggplot2)
library(raster)
library(ggcorrplot)
library(viridis)
library(dplyr)
library(tidyr)
library(usdm)
library(terra)
library(dismo)
library(biomod2)
library(wesanderson)
library(ggthemes)

# Dependencies for biomod2 modeling
# uncomment to install if necessary:
# install.packages(c('mda', 'gam', 'earth', 'maxnet', 'randomForest', 'xgboost', 'ecospat'))
library(ecospat)

#### ENVIRONMENT SETTING ####
# Set the random seed to ensure reproducibility of results
set.seed(2222) 

# Set your working directory below:
# setwd("path/to/your/repository")

#### LOAD OCCURRENCE DATA ####
occurrence_data <- read_sf('1_OCCURENCES/species_occ.shp')

#### LOAD STUDY AREA ####
study_area <- read_sf('1_OCCURENCES/study_area.shp') # Add your study area file path here
plot(study_area)

#### LOAD ENVIRONMENTAL VARIABLES ####
# Load all raster layers from the environmental variables folder
env_variables <- stack(list.files('ENV_VARIABLES', pattern = ".tif", full.names = TRUE))
env_variables <- mask(env_variables, study_area)

# Here you should define the names of the the variables you'll use for modeling. The ones bellow are the ones I used
names(env_variables) <- c(
  "land_cover", "slope", "dist_protected_areas", "dist_roads", "elevation", "forests",
  "ghm", "habitat_heterogeneity", "perc_agropastoral", "perc_forests", "perc_other_forests", "annual_prec",
  "prec_dry_month", "prec_seasonality", "prec_wet_month", "max_temp_warm_month", "mean_temp",
  "min_temp_cold_month", "temp_seasonality", "topographic_diversity", "water"
)

#### ENVIRONMENTAL VARIABLES CORRELATION ####
# STEP 1: CORRELATION - Calculate correlation matrix using Pearson's coefficient
corMatrix <- as.data.frame(layerStats(env_variables, 'pearson', na.rm = TRUE)$`pearson correlation coefficient`)

# STEP 2: Plot correlogram
corr_plot <- ggcorrplot(
  corMatrix, 
  hc.order = TRUE, 
  lab = TRUE, 
  ggtheme = ggplot2::theme_light,
  show.legend = TRUE, 
  legend.title = "Pearson's r", 
  show.diag = FALSE,
  colors = inferno(3),
  outline.color = "white",
  lab_size = 2, 
  pch = .1, 
  pch.col = "black", 
  pch.cex = .1,
  tl.cex = 5, 
  digits = 2
) +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, hjust = 0.95, vjust = 0.5), 
    axis.text.y = element_text(size = 7, hjust = 0.95, vjust = 0.5)
  )

# Display the plot
print(corr_plot)

# STEP 3: SAVE OUTPUTS
# Ensure the '4_RESULTS' directory exists before saving
ggsave("4_RESULTS/correlationMatrix.png", plot = corr_plot, scale = 1.25, dpi = 600)
write.csv(corMatrix, "4_RESULTS/correlationMatrix.csv", row.names = TRUE)

#### IDENTIFY HIGHLY CORRELATED VARIABLES ####
# STEP 1: Filter variable pairs with high correlation (|r| > 0.7)
corVariables <- layerStats(env_variables, 'pearson', na.rm = TRUE)$`pearson correlation coefficient` %>%
  as.data.frame() %>%
  mutate(var1 = rownames(.)) %>%
  gather(var2, value, -var1) %>%
  arrange(desc(value)) %>%
  filter(abs(value) > 0.7 & abs(value) < 0.999999)

# STEP 2: Count how many times each variable is highly correlated
corVariablesNumber <- as.data.frame(
  table(corVariables$var1)[order(table(corVariables$var1), decreasing = TRUE)]
)

#### VARIABLE REMOVAL (CORRELATION) ####
# CRITICAL USER STEP: Check 'corVariables' and 'corVariablesNumber' above.
# Identify pairs with correlation |r| > 0.7. From each pair, choose one 
# variable to exclude (e.g., the one with less ecological relevance).
# Replace the empty vector below with the names of the variables you want to drop.
toRemove <- c() # Example: toRemove <- c("elevation", "mean_temp")

if (length(toRemove) > 0) {
  env_variables <- dropLayer(env_variables, toRemove)
}

# Visual check of remaining layers
plot(env_variables)

#### MULTICOLLINEARITY ANALYSIS (VIF) ####
# STEP 1: INITIAL VIF CALCULATION
vifCollinearity <- vif(as.data.frame(env_variables))
write.csv(vifCollinearity, '4_RESULTS/vifCollinearity_raw.csv', row.names = FALSE)

# STEP 2: MULTICOLLINEARITY REDUCTION (OPTIONAL ROUND 1)
# Check the 'vifCollinearity' output. IF any variable still has VIF > 4 
# (or high correlation), enter its name below to remove it. Otherwise, leave it empty.
toRemove2 <- c() # Example: toRemove2 <- c("min_temp_cold_month")

if (length(toRemove2) > 0) {
  env_variables <- dropLayer(env_variables, toRemove2)
  
  # Recalculate and save VIF only if variables were removed
  vifCollinearity <- vif(as.data.frame(env_variables))
  write.csv(vifCollinearity, '4_RESULTS/vifCollinearity_final16.csv', row.names = FALSE)
}

# STEP 3: MULTICOLLINEARITY REDUCTION (OPTIONAL ROUND 2)
# Repeat the process IF necessary. Otherwise, leave it empty.
toRemove3 <- c() # Example: toRemove3 <- c("annual_prec")

if (length(toRemove3) > 0) {
  env_variables <- dropLayer(env_variables, toRemove3)
  
  # Final VIF calculation and export only if variables were removed
  vifCollinearity <- vif(as.data.frame(env_variables))
  write.csv(vifCollinearity, '4_RESULTS/vifCollinearity_final15.csv', row.names = FALSE)
}

#### SAVE FINAL SELECTED RASTERS ####
# Export the final clean set of environmental variables
terra::writeRaster(rast(env_variables), "4_RESULTS/selectedVariables.tif", overwrite = TRUE)


#### MODELING ####

#### SPATIAL FILTERING ####
# Thin occurrence data to retain only one point per raster cell (pixel)
s <- gridSample(data.frame(x = occ_species$long, y = occ_species$lat), env_variables[[1]], n = 1) 
occ_species <- occ_species[row.names(s), ]

# Adjust presence dataframe structure
occ_species$presence <- 1
occ_species["specie"] <- NULL
occ_species$species <- 'Nasua_nasua' #Replace for the name of your species of interest

###############Acho que dá pra tirar né. Pelo que entendi é só ajustar o nome da coluna aqui###############
# Standardize column name formatting
colnames(occ_species)[colnames(occ_species) == "Ponto"] <- "ponto"

# Save filtered occurrences
write.csv(occ_species, '1_OCCURENCES/nasua_occurrences_filtered.csv', row.names = FALSE)

#### BACKGROUND POINTS GENERATION ####
# Sample 10,000 random background points within the study area
back_data <- st_sample(
  study_area,
  10000,
  type = "random",
  exact = TRUE,
  warn_if_not_integer = TRUE,
  by_polygon = FALSE,
  progress = TRUE
) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  st_as_sf(coords = c("X", "Y"), remove = FALSE)

# Merge presence and background data frames
back_data <- data.frame(
  "specie" = unique(occ_species$species), 
  "occurence" = NA, 
  "lat" = back_data$Y, 
  "long" = back_data$X, 
  "ref" = NA, 
  "method" = NA, 
  "geometry" = back_data$geometry, 
  "presence" = NA
)

nasua <- rbind(as.data.frame(occ_species), back_data)

#### BIOMOD2 DATA FORMATTING ####
nasua_df <- as.data.frame(nasua)
nasua_vect <- vect(nasua_df, geom = c('long', 'lat'))

# Format data for biomod2 packaging

myBiomodData <- BIOMOD_FormatingData(
  resp.var = nasua_df$presence,
  expl.var = stack(env_variables),
  resp.xy = nasua_df[, c("long", "lat")],
  resp.name = unique(nasua_df$specie),
  PA.nb.rep = 0,
  PA.nb.absences = 10000,
  PA.dist.min = 1000, 
  PA.strategy = 'random'
)

#### EXPORT MODELING DATA TABLE ####
# Extract environmental variables at presence and pseudo-absence locations
dataTable <- as.data.frame(myBiomodData@data.env.var)
dataTable <- cbind(myBiomodData@coord, dataTable)
dataTable$Presence <- myBiomodData@data.species
dataTable$Presence <- ifelse(is.na(dataTable$Presence), "Pseudo-absence", "Presence")
dataTable$Species <- "Nasua nasua"

# Save final dataset
write.csv(dataTable, "4_RESULTS/myBiomodData_final.csv", row.names = FALSE)

#### VISUALIZATION: DENSITY PLOTS ####

# STEP 1: PREPARE DATA FOR FACETING
# Automatically select environmental variable columns by excluding metadata
metadata_cols <- c("x", "y", "Presence", "Species")
var_names <- setdiff(names(dataTable), metadata_cols)

dataTable2 <- dataTable %>%
  dplyr::select(all_of(var_names)) %>%
  gather(key = "text", value = "value") %>%
  mutate(
    text = gsub("_", " ", text),
    value = as.numeric(value)
  )

# Add the presence/pseudo-absence label back for grouping
dataTable2$presence <- rep(dataTable$Presence, length(var_names))

# STEP 2: GENERATE DENSITY PLOT
gg1 <- ggplot(dataTable2, aes(x = value, fill = presence)) +
  geom_density(colour = "black", alpha = .7, linewidth = .1) +
  scale_fill_manual(
    values = wes_palette(n = 3, name = 'Darjeeling1'), 
    breaks = c("Pseudo-absence", "Presence")
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text = element_text(size = 5),
    strip.text = element_text(size = 7) 
  ) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~ text, scales = 'free')

# Display the plot
print(gg1)

# STEP 3: SAVE PLOT
ggsave("4_RESULTS/density_distribution_presence_absence.png", plot = gg1, dpi = 600)

#### ENSEMBLE OF SMALL MODELS (ESM) ####

# STEP 1: CALIBRATION OF SIMPLE BIVARIATE MODELS
set.seed(2222)
my.ESM <- ecospat.ESM.Modeling(
  data = myBiomodData,
  models = c('GLM', "RF", "ANN", "MARS", "FDA", "CTA", "SRE", "XGBOOST"),
  NbRunEval = 1,
  DataSplit = 80,
  weighting.score = c("TSS"),
  parallel = TRUE,
  modeling.id = 'AllModels',
  Prevalence = NULL
)

# STEP 2: ENSEMBLE ECOSPAT MODELS
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM, weighting.score = c("SomersD"), threshold = 0)

# STEP 3: CALCULATE THRESHOLDS TO PRODUCE BINARY MAPS
my.ESM_thresholds <- ecospat.ESM.threshold(my.ESM_EF)


#### STANDARD SINGLE MODELS CALIBRATION ####
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  modeling.id = 'AllModels',
  models = c("GLM", "RF", "ANN", "MARS", "FDA", "CTA", "SRE", "XGBOOST"),
  CV.nb.rep = 100,
  data.split.perc = 80,
  var.import = 0,
  metric.eval = c('TSS'),
  CV.do.full.models = FALSE,
  prevalence = NULL,
  seed.val = 2222,
  nb.cpu = 18
)


#### SINGLE MODELS EVALUATION PLOTS ####

# Plot 2: Evaluation by algorithm
gg2 <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
gg2_1 <- gg2$plot + 
  theme_light() +
  scale_fill_brewer(palette = 'Accent') +
  geom_hline(data = data.frame(yint = 0.8, Eval.metric = "TSS"), aes(yintercept = yint), linetype = "dotted") +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(gg2_1)

# Plot 3: Evaluation by run
gg3 <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('run', 'run'))
gg3 <- gg3$plot + 
  theme_light() +
  scale_color_viridis(discrete = TRUE, option = 'plasma') +
  geom_hline(data = data.frame(yint = 0.8, Eval.metric = "TSS"), aes(yintercept = yint), linetype = "dotted") +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(gg3)

# Save Evaluation Plots
ggsave("4_RESULTS/gg2_1_evaluation_algorithms.png", plot = gg2_1, dpi = 600, width = 4199, height = 3199, units = "px")
ggsave("4_RESULTS/gg3_evaluation_runs.png", plot = gg3, dpi = 600, width = 19000, height = 8000, units = "px", limitsize = FALSE)


#### BIOMOD2 ENSEMBLE MODELING ####

# Load previous single models output
model_out_path <- "Nasua.nasua/Nasua.nasua.AllModels.models.out"
myBiomodModelOut2 <- assign("myBiomodModelOut2", get(load(model_out_path)))

myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut2,
  models.chosen = 'all',
  em.by = 'all',
  metric.select = c('TSS'),
  metric.select.thresh = c(0.8),
  var.import = 3,
  metric.eval = c('TSS'),
  em.algo = 'EMwmean',
  EMwmean.decay = 'proportional',
  seed.val = 2222
)


#### ENSEMBLE EVALUATION & VARIABLE IMPORTANCE METRICS ####

# Export evaluation scores
EMevaluations <- get_evaluations(myBiomodEM)
write.csv(EMevaluations, "4_RESULTS/EMevaluations_final.csv", row.names = FALSE)

EMevaluations <- EMevaluations[, c(1, 3:10)]
EMevaluations$algo <- "EM"

# Export variable importance (detailed and summarized metrics)
EMvarImportance <- get_variables_importance(myBiomodEM)
write.csv(EMvarImportance, "4_RESULTS/EMvarImportance_run_final.csv", row.names = FALSE)

EMvarImportance_mean <- EMvarImportance %>%
  group_by(expl.var) %>%
  summarize(m = mean(var.imp)) %>%
  arrange(desc(m), expl.var)

write.csv(EMvarImportance_mean, "4_RESULTS/EMvarImportance_mean_final.csv", row.names = FALSE)


#### ENSEMBLE VISUALIZATION PLOTS ####

# Plot 4: Combined Model + Ensemble Evaluation
gg2_2 <- gg2
EMevaluations3 <- data.frame(
  'full.name' = EMevaluations$full.name, 
  'PA' = 'allData',
  'run' = 'RUN1',
  'algo' = EMevaluations$algo,
  'metric.eval' = EMevaluations$metric.eval,
  'cutoff' = EMevaluations$cutoff,
  'sensitivity' = as.numeric(EMevaluations$sensitivity),
  'specificity' = as.numeric(EMevaluations$specificity),
  'calibration' = as.numeric(EMevaluations$sensitivity) / 100,
  'validation' = NA,
  'evaluation' = NA
)

gg2_2$plot$data <- rbind(gg2_2$plot$data, EMevaluations3)
gg2_2 <- gg2_2$plot + 
  theme_light() + 
  geom_hline(data = data.frame(yint = 0.8, Eval.metric = "TSS"), aes(yintercept = yint), linetype = "dotted") +
  scale_fill_brewer(palette = 'Accent') +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(gg2_2)

# Plot 5: Ensemble Variable Importance
gg4 <- ggplot() +
  geom_bar(data = EMvarImportance_mean, aes(x = reorder(expl.var, m), y = m, fill = reorder(expl.var, m)), stat = "identity") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.7, direction = 1, option = 'inferno') +
  coord_flip() +
  theme_light() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11),
    plot.margin = unit(c(0.1, 0.3, 0, -0.3), "cm")
  ) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

print(gg4)

# Plot 6: Response Curves
gg5 <- bm_PlotResponseCurves(
  bm.out = myBiomodEM,
  models.chosen = 'all',
  fixed.var = 'median'
)$plot

gg5$data <- subset(gg5$data, pred.name == "Nasua.nasua_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo")

gg5 <- gg5 + 
  theme_light() + 
  geom_line(aes(color = expl.name)) +
  scale_color_viridis(discrete = TRUE, alpha = 1, direction = -1, option = 'inferno') +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "none"
  ) +
  facet_wrap(~expl.name, scales = 'free', ncol = 3) +
  labs(x = NULL, y = 'Probability of occurrence') + 
  ggtitle("")

print(gg5)

# Save Ensemble Plots & Data
ggsave("4_RESULTS/gg2_2_combined_evaluation.png", plot = gg2_2, dpi = 600, width = 4199, height = 3199, units = "px")
ggsave("4_RESULTS/gg4_variable_importance.png", plot = gg4, dpi = 600, width = 4199, height = 4199, units = "px")
ggsave("4_RESULTS/gg5_response_curves.png", plot = gg5, dpi = 600, width = 4199, height = 4199, units = "px")

write.csv(gg5$data, "4_RESULTS/EMresponseCurves.csv", row.names = FALSE)


#### ENSEMBLE FORECASTING (PROJECTION) ####

# Project onto current environmental conditions
myBiomodProj <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  proj.name = 'CurrentEM',
  new.env = env_variables,
  models.chosen = 'all',
  metric.binary = c('TSS'),
  metric.filter = c('TSS'),
  build.clamping.mask = TRUE,
  output.format = '.img',
  seed.val = 2222
)
