###############################################################################
# Field evaluation of 3-(N-acetyl-n-butyl) aminopropionic acid ethyl ester
# (IR3535) as a spatial repellent to control malaria:
# A Randomised, Before-After-Control-Intervention trial
#
# R script: geospatial analysis of entomological data
#
# Purpose:
#   1. Load and prepare mosquito collection and household/residency data.
#   2. Aggregate mosquito counts by household/location, collection round,
#      treatment arm, indoor/outdoor setting and species.
#   3. Fit a negative-binomial Bayesian geospatial model using INLA/SPDE.
#   4. Add a spillover-gradient covariate for control households near treated
#      areas/locations and refit the model.
#   5. Produce fitted-value maps, post-pre change maps, spillover-gradient
#      diagnostic plots and fixed-effect forest plots.
#
# Inputs expected:
#   - BDEnto.xlsx: entomological records.
#   - AF_Resid.xlsx: household/residency data with AFnovo and NumRes.
#   - studyAreaControl.shp: control-area polygon.
#   - studyAreaTreat.shp: repellent/treatment-area polygon.
###############################################################################

#install.packages("BiocManager")
#BiocManager::install(c("graph", "Rgraphviz"), ask = FALSE, update = TRUE)
#options(pkgType = "binary")
#install.packages(
#  "INLA",
#  repos = c("https://cloud.r-project.org",
#            INLA = "https://inla.r-inla-download.org/R/stable"),
#  dependencies = c("Depends", "Imports", "LinkingTo")
#)
# install.packages('dplyr')
# install.packages('readxl')
# install.packages('sf')
# install.packages('spdep')
# install.packages('sp')

library(INLA)
#inla.version()
library(dplyr)
library(readxl)
library(sf)
library(sp)
library(fmesher)
library(tmap)
library(ggplot2)
library(raster)
library(terra)

# Load and merge data
ento_file_path <- "BDEnto.xlsx"
entodf <- read_excel(ento_file_path)

resid_data <- read_excel("AF_Resid.xlsx")
entodf <- left_join(entodf, resid_data, by = "AFnovo")

# Filter + aggregate
agg_df <- entodf %>%
  filter(
    !is.na(lat), !is.na(lon), !is.na(Ncol), !is.na(tratPulv),
    !is.na(inOut), !is.na(genus), !is.na(NumRes),
    sp != 7
  ) %>%
  mutate(
    post = ifelse(Ncol > 1, 1, 0),
    treated_time = tratPulv * post,
    sp = factor(sp, levels = c(1, 2, 3, 4),
                labels = c("An_funestus", "An_gambiae", "Culex", "An_tenebrosus")),
    inOut = factor(inOut, levels = c(1, 2), labels = c("Indoor", "Outdoor")),
    genus = factor(genus, levels = c(1, 3), labels = c("Anopheles", "Culex"))
  ) %>%
  group_by(lat, lon, Ncol, tratPulv, post, treated_time, sp, inOut, genus) %>%
  summarise(
    count = n(),
    mean_NumRes = mean(NumRes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    norm_count = count / mean_NumRes,
    intercept = 1,
    y = count  # Response variable
  )

#Check dataframe
agg_df

# Create spatial mesh
coords_sf <- st_as_sf(agg_df, coords = c("lon", "lat"), crs = 4326)
coords_xy <- st_coordinates(coords_sf)

mesh <- fmesher::fm_mesh_2d_inla(
  loc = coords_xy,
  max.edge = c(0.02, 0.2),
  cutoff = 0.01
)

# SPDE and A matrix
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(0.02, 0.01),
                            prior.sigma = c(1, 0.01))
s_index <- inla.spde.make.index("space", spde$n.spde)

coords_sf <- sf::st_as_sf(agg_df, coords = c("lon", "lat"), crs = 4326)
coords_xy <- sf::st_coordinates(coords_sf)[, 1:2, drop = FALSE]  # X/Y numeric matrix

mesh <- fmesher::fm_mesh_2d_inla(
  loc = coords_xy,
  max.edge = c(0.02, 0.2),
  cutoff = 0.01
)

A <- INLA::inla.spde.make.A(mesh = mesh, loc = coords_xy)

# Create design matrix (auto dummy encoding)
X <- model.matrix(~ 0 + intercept + tratPulv + post + treated_time +
                    inOut + sp +
                    treated_time:sp +
                    inOut:sp +
                    treated_time:inOut,
                  data = agg_df)

# Create INLA Stack
stack <- inla.stack(
  data = list(y = agg_df$y, E = agg_df$mean_NumRes),
  A = list(A, 1),
  effects = list(
    space = s_index,
    cbind(as.data.frame(X), time = agg_df$Ncol)
  ),
  tag = "est"
)

# Define model formula
formula <- y ~ -1 + intercept + tratPulv + post + treated_time +
  inOutOutdoor +
  spAn_gambiae + spCulex + spAn_tenebrosus +
  treated_time:spAn_gambiae + treated_time:spCulex + treated_time:spAn_tenebrosus +
  inOutOutdoor:spAn_gambiae + inOutOutdoor:spCulex + inOutOutdoor:spAn_tenebrosus +
  treated_time:inOutOutdoor +
  f(space, model = spde) +
  f(time, model = "rw1")

# Run INLA model
result_combined <- inla(
  formula,
  family = "nbinomial",
  data = inla.stack.data(stack),
  E = inla.stack.data(stack)$E,
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  control.inla = list(strategy = "simplified.laplace")
)

# Review model
summary(result_combined)
plot(result_combined)
names(result_combined)



################################################################################
###### SPILOVER MODEL ##########
################################################################################

###############################################################################
# Add an explicit spillover/gradient covariate
#
# It will:
#   1) Build a proximity-to-opposite-arm gradient (in METERS)
#   2) Add spillover covariate(s) to the dataset
#   3) Refit the INLA model including the spillover term
#
# IMPORTANT:
# - This code creates a NEW dataset and refits a NEW model object
# - It reuses the existing mesh/SPDE/A matrix from the baseline run so the
#   modeling code stays as close as possible to the previous model.
###############################################################################


# Define spillover/gradient settings 
r_m <- 250

# Ensure lon/lat exist as numeric columns 
stopifnot(all(c("lon", "lat", "tratPulv", "post") %in% names(agg_df)))

# Build an sf object from the aggregated dataset 
agg_sf_ll <- st_as_sf(
  agg_df,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# The median lon/lat is uded to represent the study area.
lon_med <- stats::median(agg_df$lon, na.rm = TRUE)
lat_med <- stats::median(agg_df$lat, na.rm = TRUE)
utm_zone <- floor((lon_med + 180) / 6) + 1
epsg_utm <- ifelse(lat_med < 0, 32700 + utm_zone, 32600 + utm_zone)

# Transform points to the chosen UTM CRS (units become meters).
agg_sf_m <- st_transform(agg_sf_ll, crs = epsg_utm)

# Split into treated and control subsets (in metric CRS)
sf_treat_m <- agg_sf_m[agg_sf_m$tratPulv == 1, ]
sf_ctrl_m  <- agg_sf_m[agg_sf_m$tratPulv == 0, ]

# both groups must exist
if (nrow(sf_treat_m) == 0 || nrow(sf_ctrl_m) == 0) {
  stop("Spillover gradient requires BOTH treated (tratPulv==1) and control (tratPulv==0) rows.")
}

# For each treated row, find nearest control row index
idx_near_ctrl_for_treat <- st_nearest_feature(sf_treat_m, sf_ctrl_m)

# Compute treated-to-nearest-control distances (meters)
d_treat_to_ctrl_m <- as.numeric(
  st_distance(sf_treat_m, sf_ctrl_m[idx_near_ctrl_for_treat, ], by_element = TRUE)
)

# For each control row, find nearest treated row index
idx_near_treat_for_ctrl <- st_nearest_feature(sf_ctrl_m, sf_treat_m)

# Compute control-to-nearest-treated distances (meters)
d_ctrl_to_treat_m <- as.numeric(
  st_distance(sf_ctrl_m, sf_treat_m[idx_near_treat_for_ctrl, ], by_element = TRUE)
)

# Create a full-length distance vector aligned with agg_df row order
dist_to_opposite_m <- rep(NA_real_, nrow(agg_df))

# Fill treated rows with treated-to-control distances
dist_to_opposite_m[agg_df$tratPulv == 1] <- d_treat_to_ctrl_m

# Fill control rows with control-to-treated distances
dist_to_opposite_m[agg_df$tratPulv == 0] <- d_ctrl_to_treat_m

# Convert distance into a smooth spillover weight
# exponential decay
#   w = exp(-d / r_m)
# where:
#   d   = distance to nearest opposite-arm location (meters)
#   r_m = decay scale (meters)
# Properties:
# - d=0 => w=1
# - d=r_m => w=exp(-1) ~ 0.37
# - as d increases => w approaches 0
w_boundary <- exp(-dist_to_opposite_m / r_m)

# Define the actual spillover covariates 
# We encode "spillover exposure" for CONTROLS only:
# - controls close to treated locations may be partially exposed
# - treated locations get 0 for spill_ctrl (because they are "fully treated")
# This covariate is NOT the treatment effect - it is a contamination term.
spill_ctrl <- ifelse(agg_df$tratPulv == 0, w_boundary, 0)

# allow spillover to matter mainly after intervention starts.
spill_ctrl_post <- spill_ctrl * agg_df$post

# Create a new dataset for the spillover model
agg_df_spill <- agg_df

# Add distance and spillover columns for documentation + modeling.
agg_df_spill$dist_to_opposite_m <- dist_to_opposite_m
agg_df_spill$w_boundary         <- w_boundary
agg_df_spill$spill_ctrl         <- spill_ctrl
agg_df_spill$spill_ctrl_post    <- spill_ctrl_post

# Build a new design matrix INCLUDING spillover terms
X_spill <- model.matrix(
  ~ 0 + intercept + tratPulv + post + treated_time +
    spill_ctrl +
    spill_ctrl_post +
    inOut + sp +
    treated_time:sp +
    inOut:sp +
    treated_time:inOut,
  data = agg_df_spill
)

# Rebuild the INLA stack for the spillover model 
# Reuse:
#   - A (the projector matrix) from the baseline code
#   - s_index from the baseline SPDE setup
#   - the same time index (Ncol)
stack_spill <- inla.stack(
  data = list(y = agg_df_spill$y, E = agg_df_spill$mean_NumRes),
  A = list(A, 1),
  effects = list(
    space = s_index,
    cbind(as.data.frame(X_spill), time = agg_df_spill$Ncol)
  ),
  tag = "est_spill"
)

# Define the INLA formula for the spillover model
# start from the baseline fixed effects and add:
#   + spill_ctrl
#   + spill_ctrl_post

# Base fixed effects part (copied from the baseline model formula)
base_formula_text <- paste(
  "y ~ -1 + intercept + tratPulv + post + treated_time",
  "+ inOutOutdoor",
  "+ spAn_gambiae + spCulex + spAn_tenebrosus",
  "+ treated_time:spAn_gambiae + treated_time:spCulex + treated_time:spAn_tenebrosus",
  "+ inOutOutdoor:spAn_gambiae + inOutOutdoor:spCulex + inOutOutdoor:spAn_tenebrosus",
  "+ treated_time:inOutOutdoor",
  sep = " "
)

# Spillover terms
spill_terms <- c()
if ("spill_ctrl" %in% colnames(X_spill)) {
  spill_terms <- c(spill_terms, "spill_ctrl")
}
if ("spill_ctrl_post" %in% colnames(X_spill)) {
  spill_terms <- c(spill_terms, "spill_ctrl_post")
}

spill_formula_text <- if (length(spill_terms) > 0) {
  paste("+", paste(spill_terms, collapse = " + "))
} else {
  ""
}

# Random effects part (same as baseline)
re_formula_text <- paste(
  "+ f(space, model = spde)",
  "+ f(time, model = 'rw1')"
)

# Combine into final formula text and convert to an R formula
formula_spill <- as.formula(paste(base_formula_text, spill_formula_text, re_formula_text))

# Fit the spillover model 
# same settings as the baseline run for comparability.
result_spillover <- inla(
  formula_spill,
  family = "nbinomial",
  data = inla.stack.data(stack_spill),
  E = inla.stack.data(stack_spill)$E,
  control.predictor = list(A = inla.stack.A(stack_spill), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  control.inla = list(strategy = "simplified.laplace")
)

# Inspect results 
summary(result_spillover)

# Optional: directly compare key fit criteria vs baseline
cat("\n--- Fit criteria comparison ---\n")
cat("Baseline WAIC:", result_combined$waic$waic, "\n")
cat("Spillover WAIC:", result_spillover$waic$waic, "\n")
cat("Baseline DIC :", result_combined$dic$dic, "\n")
cat("Spillover DIC :", result_spillover$dic$dic, "\n")

# Optional: check CPO failures (should ideally be 0)
cat("\n--- CPO failures ---\n")
cat("Baseline failures :", sum(result_combined$cpo$failure, na.rm = TRUE), "\n")
cat("Spillover failures:", sum(result_spillover$cpo$failure, na.rm = TRUE), "\n")

# Optional: show spillover coefficient(s) quickly if present
if ("spill_ctrl" %in% rownames(result_spillover$summary.fixed)) {
  cat("\n--- Spillover effect (spill_ctrl) ---\n")
  print(result_spillover$summary.fixed["spill_ctrl", ])
}
if ("spill_ctrl_post" %in% rownames(result_spillover$summary.fixed)) {
  cat("\n--- Spillover post interaction (spill_ctrl_post) ---\n")
  print(result_spillover$summary.fixed["spill_ctrl_post", ])
}


################################################################################
# Load Study Area Shapefiles
################################################################################

control_area_sf <- st_read("studyAreaControl.shp")
treated_area_sf <- st_read("studyAreaTreat.shp")

control_area_sf_poly <- st_make_valid(control_area_sf)
treated_area_sf_poly <- st_cast(treated_area_sf, "POLYGON") %>% st_make_valid()

control_area_sf_poly$area_type <- "Control"
treated_area_sf_poly$area_type <- "Repellent"

study_areas <- rbind(
  control_area_sf_poly[, c("geometry", "area_type")],
  treated_area_sf_poly[, c("geometry", "area_type")]
)

# Extract predictions from the spillover model
idx_est_spill <- inla.stack.index(stack_spill, tag = "est_spill")$data

lambda_mean_spill <- result_spillover$summary.fitted.values$mean[idx_est_spill]
lambda_q025_spill <- result_spillover$summary.fitted.values$`0.025quant`[idx_est_spill]
lambda_q975_spill <- result_spillover$summary.fitted.values$`0.975quant`[idx_est_spill]

# Build ONE sf object that matches exactly the estimation stack rows
lambda_df_spill <- agg_df_spill %>%
  mutate(
    lambda_mean = lambda_mean_spill,
    lambda_q025 = lambda_q025_spill,
    lambda_q975 = lambda_q975_spill,
    # per-res version (protect against division issues)
    lambda_per_res = lambda_mean / ifelse(mean_NumRes > 0, mean_NumRes, NA_real_)
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# ----------------------------
# Plot Map (Predicted mosquito density)
# ----------------------------
tmap_mode("plot")

tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill_alpha = 0.3,
    fill.scale  = tm_scale(values = c("Control" = "skyblue", "Repellent" = "orange")),
    fill.legend = tm_legend(title = "Study Area")
  ) +
  tm_shape(lambda_df_spill) +
  tm_symbols(
    fill = "lambda_per_res",  # or "lambda_mean" if you prefer the raw count scale
    fill.scale  = tm_scale_intervals(style = "quantile", values = "viridis"),
    size = 0.6,
    fill.legend = tm_legend(title = expression(Predicted~mu))
  ) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_title("Predicted Mosquito Density (Spillover-adjusted model)") +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.text.size = 0.9,
    legend.title.size = 1.0,
    legend.width = 8
  )


################################################################################
## OTHER MAPS
################################################################################

###############################################################################
# POST − PRE CHANGE MAP (spillover-adjusted model)
###############################################################################

# Helper: build the SAME design matrix used in the spillover model 
# This must match the model.matrix() call used when fitting result_spillover.
build_X_spill <- function(df) {
  X <- model.matrix(
    ~ 0 + intercept + tratPulv + post + treated_time +
      inOut + sp +
      treated_time:sp +
      inOut:sp +
      treated_time:inOut +
      spill_ctrl + spill_ctrl_post,
    data = df
  )
  return(X)
}

# Create PRE and POST prediction datasets (same rows, different scenario)
pred_pre <- agg_df_spill %>%
  mutate(
    post = 0,
    treated_time = 0,
    # spill_ctrl is distance-based and does NOT depend on post; keep as is
    spill_ctrl_post = 0,
    # use a consistent time index for baseline (pre = Ncol 1)
    Ncol = 1,
    intercept = 1
  )

pred_post <- agg_df_spill %>%
  mutate(
    post = 1,
    treated_time = tratPulv * 1,
    # spillover interaction activates post
    spill_ctrl_post = spill_ctrl * 1,
    # use a consistent time index for first post (post = Ncol 2)
    Ncol = 2,
    intercept = 1
  )

# Build A matrices for prediction points (same coordinates as data) ----
# We create sf points to ensure lon/lat are numeric and properly used.
pred_coords_sf <- st_as_sf(agg_df_spill, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
pred_coords_xy <- st_coordinates(pred_coords_sf)

A_pred <- inla.spde.make.A(mesh = mesh, loc = pred_coords_xy)

# Build prediction stacks
# Build design matrices
X_pre  <- build_X_spill(pred_pre)
X_post <- build_X_spill(pred_post)

# Create stacks (no response; these are prediction-only)
stack_pre <- inla.stack(
  data = list(y = NA, E = pred_pre$mean_NumRes),
  A = list(A_pred, 1),
  effects = list(
    space = s_index,
    cbind(as.data.frame(X_pre), time = pred_pre$Ncol)
  ),
  tag = "pred_pre"
)

stack_post <- inla.stack(
  data = list(y = NA, E = pred_post$mean_NumRes),
  A = list(A_pred, 1),
  effects = list(
    space = s_index,
    cbind(as.data.frame(X_post), time = pred_post$Ncol)
  ),
  tag = "pred_post"
)

# Combine original estimation stack + prediction stacks
stack_all <- inla.stack(stack_spill, stack_pre, stack_post)

# Refit to obtain predictions for the new rows
# reuse the spillover model formula
result_spill_pp <- inla(
  formula_spill,  # must be the SAME formula used for result_spillover
  family = "nbinomial",
  data = inla.stack.data(stack_all),
  E = inla.stack.data(stack_all)$E,
  control.predictor = list(A = inla.stack.A(stack_all), compute = TRUE),
  control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
  control.mode = list(theta = result_spillover$mode$theta, restart = FALSE),
  control.inla = list(strategy = "simplified.laplace")
)

# Extract predicted means for PRE and POST rows 
idx_pre  <- inla.stack.index(stack_all, tag = "pred_pre")$data
idx_post <- inla.stack.index(stack_all, tag = "pred_post")$data

mu_pre  <- result_spill_pp$summary.fitted.values$mean[idx_pre]
mu_post <- result_spill_pp$summary.fitted.values$mean[idx_post]

# Convert to per-residency rate
E_i <- agg_df_spill$mean_NumRes
rate_pre  <- mu_pre  / E_i
rate_post <- mu_post / E_i

# Difference and ratio
delta_rate <- rate_post - rate_pre
ratio_rate <- rate_post / pmax(rate_pre, 1e-9)  # guard against divide-by-zero

# Build sf object for mapping
change_df <- agg_df_spill %>%
  mutate(
    rate_pre = rate_pre,
    rate_post = rate_post,
    delta_rate = delta_rate,
    ratio_rate = ratio_rate
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Plot maps
tmap_mode("plot")

# Difference map (post - pre)
tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill.scale  = tm_scale(values = c("Control" = "skyblue", "Repellent" = "orange")),
    fill.alpha = 0.25,
    fill.legend = tm_legend(title = "Study Area")
  ) +
  tm_shape(change_df) +
  tm_symbols(
    fill = "delta_rate",
    size = 0.6,
    fill.scale = tm_scale_intervals(style = "quantile"),
        fill.legend = tm_legend(title = "Δ predicted density\n(post − pre)\n(per residency)")
  ) +
  tm_title("Change in predicted mosquito density (spillover model)") +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "right",
            legend.text.size = 0.9,
            legend.title.size = 1.0,
            legend.width = 10)

# Ratio map (post / pre) — optional
tm_shape(study_areas) +
  tm_polygons(fill = "area_type", fill.alpha = 0.25, fill.legend = tm_legend(title = "Study Area")) +
  tm_shape(change_df) +
  tm_symbols(
    fill = "ratio_rate",
    size = 0.6,
    fill.scale = tm_scale_intervals(style = "quantile"),
    fill.legend = tm_legend(title = "Predicted ratio\n(post / pre)\n(per residency)")
  ) +
  tm_title("Ratio of predicted mosquito density (spillover model)") +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "right",
            legend.text.size = 0.9,
            legend.title.size = 1.0,
            legend.width = 10)


################################################################################
# MAPS AS GRADIENTS
################################################################################

# CONTINUOUS SURFACE MAPS from the INLA spillover entomology model


# CHECK REQUIRED OBJECTS
req <- c("study_areas","agg_df_spill","stack_spill","result_spillover","mesh","spde","s_index","formula_spill")
miss <- req[!sapply(req, exists)]
if (length(miss) > 0) stop("Missing objects in your environment: ", paste(miss, collapse=", "))


# SETTINGS
grid_cellsize_deg <- 0.001   # ~200–250m --> 0.001 for finer (slower)
r_m <- 250                   # MUST match the r_m used to compute spill_ctrl in your model
map_inOut <- "Indoor"        # choose: levels(agg_df_spill$inOut)
map_sp    <- "An_funestus"   # choose: levels(agg_df_spill$sp)

E_map <- 1                   # E=1 => map interpretable "per residency"
time_pre  <- 1               # baseline Ncol index used for "pre" scenario
time_post <- 2               # post Ncol index used for "post" scenario


# STUDY AREA GEOMETRY
study_areas <- st_make_valid(study_areas)
study_areas <- st_transform(study_areas, 4326)

study_union <- st_union(study_areas)

treated_geom <- study_areas %>% filter(area_type == "Repellent") %>% st_union()
control_geom <- study_areas %>% filter(area_type == "Control")   %>% st_union()

# compute distances in UTM zone 36S (EPSG:32736) - Mozambique
utm_epsg <- 32736

treated_geom_utm <- st_transform(treated_geom, utm_epsg)

# MAKE A GRID OF POINTS (CELL CENTERS) USING sf
grid_pts <- st_make_grid(
  study_union,
  cellsize = grid_cellsize_deg,
  what = "centers"
)

grid_sf <- st_sf(geometry = grid_pts) %>%
  st_intersection(study_union)

# Extract lon/lat
xy <- st_coordinates(grid_sf)
grid_sf$lon <- xy[,1]
grid_sf$lat <- xy[,2]

# Determine treated vs control membership
is_in_treated <- suppressWarnings(st_within(grid_sf, treated_geom, sparse = FALSE))
if (inherits(is_in_treated, "try-error") || ncol(is_in_treated) == 0) {
  # fallback: treat as "treated" if distance to treated geom is ~0
  d0 <- as.numeric(st_distance(st_transform(grid_sf, utm_epsg), treated_geom_utm))
  grid_sf$tratPulv <- ifelse(d0 < 1, 1, 0)
} else {
  grid_sf$tratPulv <- ifelse(is_in_treated[,1], 1, 0)
}

# Distance from each grid point to treated geometry (meters)
grid_sf_utm <- st_transform(grid_sf, utm_epsg)
grid_sf$dist_to_treated_m <- as.numeric(st_distance(grid_sf_utm, treated_geom_utm))

# spill_ctrl = exp(-dist/r_m) for CONTROLS only; treated gets 0
grid_sf$spill_ctrl <- ifelse(grid_sf$tratPulv == 0,
                             exp(-grid_sf$dist_to_treated_m / r_m),
                             0)


# BUILD PREDICTION DATA (PRE and POST)

# Ensure factor levels match fitted data
stopifnot(map_inOut %in% levels(agg_df_spill$inOut))
stopifnot(map_sp    %in% levels(agg_df_spill$sp))

grid_base <- grid_sf %>%
  st_drop_geometry() %>%
  mutate(
    inOut = factor(map_inOut, levels = levels(agg_df_spill$inOut)),
    sp    = factor(map_sp,    levels = levels(agg_df_spill$sp)),
    intercept = 1,
    mean_NumRes = E_map
  )

pred_pre <- grid_base %>%
  mutate(
    post = 0,
    treated_time = 0,
    spill_ctrl_post = 0,
    Ncol = time_pre
  )

pred_post <- grid_base %>%
  mutate(
    post = 1,
    treated_time = tratPulv,
    spill_ctrl_post = spill_ctrl,
    Ncol = time_post
  )

# match the model.matrix used in the spill model
X_formula <- ~ 0 + intercept + tratPulv + post + treated_time +
  spill_ctrl + spill_ctrl_post +
  inOut + sp +
  treated_time:sp +
  inOut:sp +
  treated_time:inOut

X_pre  <- model.matrix(X_formula, data = pred_pre)
X_post <- model.matrix(X_formula, data = pred_post)

# Align columns to the fitted model matrix
fit_cols <- colnames(result_spillover$model.matrix)

# Add missing cols (if any) with zeros
add_missing_cols <- function(X, cols_needed) {
  missing <- setdiff(cols_needed, colnames(X))
  if (length(missing) > 0) {
    for (m in missing) X <- cbind(X, setNames(rep(0, nrow(X)), m))
  }
  X[, cols_needed, drop = FALSE]
}

X_pre  <- add_missing_cols(X_pre,  fit_cols)
X_post <- add_missing_cols(X_post, fit_cols)


# BUILD A MATRIX FOR GRID
A_grid <- inla.spde.make.A(
  mesh = mesh,
  loc  = as.matrix(grid_base[, c("lon","lat")])
)


# BUILD STACKS + REFIT FOR PREDICTIONS
stack_pre <- inla.stack(
  data = list(y = NA, E = pred_pre$mean_NumRes),
  A = list(A_grid, 1),
  effects = list(
    space = s_index,
    cbind(as.data.frame(X_pre), time = pred_pre$Ncol)
  ),
  tag = "pred_pre"
)

stack_post <- inla.stack(
  data = list(y = NA, E = pred_post$mean_NumRes),
  A = list(A_grid, 1),
  effects = list(
    space = s_index,
    cbind(as.data.frame(X_post), time = pred_post$Ncol)
  ),
  tag = "pred_post"
)

stack_all <- inla.stack(stack_spill, stack_pre, stack_post)

result_grid <- inla(
  formula_spill,
  family = "nbinomial",
  data = inla.stack.data(stack_all),
  E = inla.stack.data(stack_all)$E,
  control.predictor = list(A = inla.stack.A(stack_all), compute = TRUE),
  control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE),
  control.mode      = list(theta = result_spillover$mode$theta, restart = FALSE),
  control.inla      = list(strategy = "simplified.laplace")
)

# Indices
idx_pre  <- inla.stack.index(stack_all, tag = "pred_pre")$data
idx_post <- inla.stack.index(stack_all, tag = "pred_post")$data

# Posterior means for mu (data scale). With E_map=1, mu is "per residency".
mu_pre  <- result_grid$summary.fitted.values$mean[idx_pre]
mu_post <- result_grid$summary.fitted.values$mean[idx_post]

delta_rate <- mu_post - mu_pre

# Attach back to grid_sf
grid_sf$rate_post  <- mu_post
grid_sf$rate_pre   <- mu_pre
grid_sf$delta_rate <- delta_rate

# RASTERIZE WITH terra
# Create raster template
bb <- st_bbox(study_union)
r_template <- rast(
  xmin = as.numeric(bb["xmin"]),
  xmax = as.numeric(bb["xmax"]),
  ymin = as.numeric(bb["ymin"]),
  ymax = as.numeric(bb["ymax"]),
  resolution = grid_cellsize_deg,
  crs = "EPSG:4326"
)

# Convert points to terra vector
grid_v <- vect(grid_sf)

# Rasterize (mean, since one point per cell, it’s effectively direct assignment)
r_post  <- rasterize(grid_v, r_template, field = "rate_post", fun = "mean", background = NA)
r_delta <- rasterize(grid_v, r_template, field = "delta_rate", fun = "mean", background = NA)

# Mask outside study area polygon
study_v <- vect(study_union)
r_post  <- mask(r_post,  study_v)
r_delta <- mask(r_delta, study_v)

# HOUSEHOLD POINTS: use the sf already built for mapping 
# This should be an sf object with POINT geometries and CRS = 4326 (lon/lat).
hh_sf <- lambda_df_spill

# make sure it is sf and has a CRS
stopifnot(inherits(hh_sf, "sf"))
stopifnot(!is.na(sf::st_crs(hh_sf)))

# Make sure it matches the map CRS
hh_sf <- sf::st_transform(hh_sf, sf::st_crs(study_areas))


################################################################################
# PLOT WITH tmap
################################################################################

# Map A: predicted POST density per residency
map_post <- tm_shape(study_areas) +
  tm_polygons(fill = "area_type", fill_alpha = 0.25) +
  tm_shape(r_post) +
  tm_raster(
    title = "Predicted density\n(post, per residency)",
    style = "quantile"
  ) +
  tm_shape(study_areas) +
  tm_borders(col = "black", lwd = 1) +
  tm_title(paste0("Predicted mosquito density surface (post)\n",
                  map_sp, ", ", map_inOut, ", spillover-adjusted")) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "right",
            legend.text.size = 0.9,
            legend.title.size = 1.0,
            legend.width = 10)

# Map B: change surface (post - pre)
map_delta <- tm_shape(study_areas) +
  tm_polygons(fill = "area_type", fill_alpha = 0.25) +
  tm_shape(r_delta) +
  tm_raster(
    title = "Change in density\n(post − pre, per residency)",
    style = "quantile",
    midpoint = NA   # prevents tmap from forcing midpoint=0 behaviour
  ) +
  tm_shape(study_areas) +
  tm_borders(col = "black", lwd = 1) +
  tm_title(paste0("Change in predicted density surface (post − pre)\n",
                  map_sp, ", ", map_inOut, ", spillover-adjusted")) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "right",
            legend.text.size = 0.9,
            legend.title.size = 1.0,
            legend.width = 12)

# Print maps
map_post
map_delta


tmap_mode("plot")

map_post <- tm_shape(study_areas) +
  tm_polygons(fill = "area_type", fill_alpha = 0.25) +
  
  tm_shape(r_post) +
  tm_raster(
    title = "Predicted density\n(post, per residency)",
    style = "quantile"
  ) +
  
  tm_shape(hh_sf) +
  tm_dots(col = "black", fill = "white", size = 0.4, alpha = 0.8) +
  
  # Borders LAST so they show clearly
  tm_shape(study_areas) +
  tm_borders(col = "black", lwd = 1.5) +
  
  # --- ADD NATIVE TMAP GRID FOR AXES ---
  tm_grid(
    labels.show = TRUE,
    labels.inside.frame = FALSE, # This forces the coordinates outside the box
    lines = FALSE,               # Keeps the map clean without internal grid lines
    x = lon_ticks,               # Optional: Reuses your custom tick spacing
    y = lat_ticks,               # Optional: Reuses your custom tick spacing
    labels.size = 0.7,
    labels.format = list(digits = 3) # Forces 3 decimal places
  ) +
  
  tm_title(paste0(
    "Predicted mosquito density surface (post) + household locations\n",
    map_sp, ", ", map_inOut, ", spillover-adjusted"
  )) +
  
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.width = 10,
    legend.text.size = 0.9,
    legend.title.size = 1.0,
    # You may need to increase your bottom/left margins slightly 
    # if the outside labels get cut off:
    outer.margins = c(0.05, 0.05, 0.02, 0.02), 
    clip = TRUE
  )

map_post

tmap_mode("plot")

map_delta <- tm_shape(study_areas) +
  tm_polygons(fill = "area_type", fill_alpha = 0.25) +
  
  tm_shape(r_delta) +
  tm_raster(
    title = "Change in density\n(post \u2212 pre, per residency)",
    style = "quantile",
    midpoint = NA
  ) +
  
  # OPTIONAL (recommended): show household locations on top (same as post map)
  tm_shape(hh_sf) +
  tm_dots(col = "black", fill = "white", size = 0.4, alpha = 0.8) +
  
  # Borders LAST so they show clearly
  tm_shape(study_areas) +
  tm_borders(col = "black", lwd = 1.5) +
  
  # --- NATIVE TMAP GRID FOR AXES (same as your post map) ---
  tm_grid(
    labels.show = TRUE,
    labels.inside.frame = FALSE,  # puts coordinates outside the box
    lines = FALSE,                # no internal grid lines
    x = lon_ticks,                # reuse your custom ticks
    y = lat_ticks,                # reuse your custom ticks
    labels.size = 0.7,
    labels.format = list(digits = 3)
  ) +
  
  tm_title(paste0(
    "Change in predicted mosquito density surface (post \u2212 pre) + household locations\n",
    map_sp, ", ", map_inOut, ", spillover-adjusted"
  )) +
  
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.width = 12,
    legend.text.size = 0.9,
    legend.title.size = 1.0,
    outer.margins = c(0.05, 0.05, 0.02, 0.02),
    clip = TRUE
  )

map_delta


# SAVE 
tmap_save(map_post,  "Fig_surface_post_per_res.png",  width = 8, height = 6, units = "in", dpi = 300)
tmap_save(map_delta, "Fig_surface_change_per_res.png",width = 8, height = 6, units = "in", dpi = 300)


################################################################################
# BLACK & WHITE STUDY MAP
# - Control vs Repellent areas as outlines + light fill
# - Household points: treated = solid black; control = open circle
################################################################################

# Ensure the household sf exists
stopifnot(exists("hh_sf"))
stopifnot(inherits(hh_sf, "sf"))

# Ensure we have a treatment indicator for points
# If your points use tratPulv (0/1), keep it.
# If your points use trat (0/1), adapt here.
if (!("tratPulv" %in% names(hh_sf))) {
  stop("hh_sf must contain a treatment indicator column named 'tratPulv' (0/1).")
}

# Create a point group factor for clean legend control
hh_sf$arm <- factor(
  ifelse(hh_sf$tratPulv == 1, "Repellent households", "Control households"),
  levels = c("Control households", "Repellent households")
)

# Make a B/W palette for polygons
# Control: white-ish; Repellent: light grey-ish (still B/W)
bw_fill_values <- c("Control" = "white", "Repellent" = "grey85")

tmap_mode("plot")

bw_map <- tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill_alpha = 1,
    fill.scale = tm_scale(values = bw_fill_values),
    fill.legend = tm_legend(title = "Study area")
  ) +
  # Area borders on top
  tm_shape(study_areas) +
  tm_borders(col = "black", lwd = 2) +
  # Household points: control = open circle; repellent = filled black
  tm_shape(hh_sf) +
  tm_symbols(
    # Border color for all points
    col = "black",
    # Fill differs by arm
    fill = "arm",
    # Make "Control households" white fill, "Repellent households" black fill
    fill.scale = tm_scale(values = c("Control households" = "white",
                                     "Repellent households" = "black")),
    # Point size
    size = 0.4,
    shape = 21,
    # Legend title
    fill.legend = tm_legend(title = "Households")
  ) +
  tm_compass(position = c("right", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_title("Study areas and household locations") +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.width = 10
  )

bw_map


bw_map <- tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill_alpha = 1,
    fill.scale = tm_scale(values = bw_fill_values),
    fill.legend = tm_legend(title = "Study area")
  ) +
  
  # Area borders on top
  tm_shape(study_areas) +
  tm_borders(col = "black", lwd = 2) +
  
  # Household points: ALL black
  tm_shape(hh_sf) +
  tm_symbols(
    col = "black",      # border
    fill = "black",     # interior fill
    size = 0.4,         # adjust as needed (0.25–0.5 typical)
    shape = 21          # circle
    # No legend needed for households now
  ) +
  
  tm_compass(position = c("right", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_title("Study areas and household locations") +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.width = 10
  )

bw_map

# Save the last map (repeat after each map if you want all three saved)
# tmap_save(last_map(), "Fig_post_surface.png", width = 8, height = 6, units = "in", dpi = 300)
# tmap_save(last_map(), "Fig_delta_surface.png", width = 8, height = 6, units = "in", dpi = 300)
# tmap_save(last_map(), "Fig_unc_surface.png", width = 8, height = 6, units = "in", dpi = 300)

# QC: check for repellent-labelled households inside the control polygon
# study_areas_v <- sf::st_make_valid(study_areas)
# hh_v <- sf::st_make_valid(hh_sf)
# 
# study_areas_v <- sf::st_transform(study_areas_v, sf::st_crs(hh_v))
# 
# control_poly <- study_areas_v %>%
#   dplyr::filter(area_type == "Control") %>%
#   sf::st_union()
# 
# bad_repellent_in_control <- hh_v %>%
#   dplyr::filter(arm == "Repellent households") %>%
#   dplyr::mutate(
#     in_control = as.logical(sf::st_within(., control_poly, sparse = FALSE))
#   ) %>%
#   dplyr::filter(in_control)
# 
# if (nrow(bad_repellent_in_control) > 0) {
#   warning(
#     nrow(bad_repellent_in_control),
#     " repellent-labelled household(s) fall inside the control polygon."
#   )
#   print(bad_repellent_in_control, width = Inf)
# }


################################################################################
## GRAPHS
################################################################################

# from the entomology spillover model code):
#   agg_df_spill      : aggregated data used for the spillover model (with dist_to_opposite_m, mean_NumRes, etc.)
#   result_spillover  : INLA model fit (nbinomial)
#   stack_spill       : INLA stack used to fit result_spillover

# This code creates 3 gradient plots:
#   1) Post-only (controls only)
#   2) Indoor-only (controls only)
#   3) An. funestus-only (controls only)

# Each plot is shown for a max distance window (e.g., 500 m or 1000 m).

# Extract fitted means for the estimation rows
idx_est <- inla.stack.index(stack_spill, tag = "est_spill")$data
mu_hat  <- result_spillover$summary.fitted.values$mean[idx_est]

# Attach predictions to the exact data used
grad_base <- agg_df_spill %>%
  mutate(
    mu_hat    = mu_hat,
    rate_hat  = mu_hat / mean_NumRes  # predicted density per residency
  ) %>%
  # We interpret spillover as contamination into CONTROLS
  filter(tratPulv == 0) %>%
  # Keep only rows with distance defined
  filter(!is.na(dist_to_opposite_m))

# Choose distance window 
max_dist_m <- 500
# max_dist_m <- 1000

grad_base <- grad_base %>% filter(dist_to_opposite_m <= max_dist_m)

# Helper function to plot a gradient
plot_gradient <- function(df, title_text) {
  ggplot(df, aes(x = dist_to_opposite_m, y = rate_hat)) +
    geom_point(alpha = 0.6) +
    geom_smooth(se = TRUE, method = "loess") +
    labs(
      x = "Distance to nearest treated location (m)",
      y = "Predicted mosquito density per residency",
      title = title_text
    ) +
    theme_minimal()
}

# Post-only gradient (controls, post=1)
grad_post_only <- grad_base %>%
  filter(post == 1)

p_post <- plot_gradient(
  grad_post_only,
  paste0("Spillover gradient (controls, post-only, ≤ ", max_dist_m, " m)")
)
print(p_post)

# Indoor-only gradient (controls, indoor only)
grad_indoor_only <- grad_base %>%
  filter(inOut == "Indoor")

p_indoor <- plot_gradient(
  grad_indoor_only,
  paste0("Spillover gradient (controls, indoor-only, ≤ ", max_dist_m, " m)")
)
print(p_indoor)

# An. funestus-only gradient (controls, species = An_funestus)
grad_funestus_only <- grad_base %>%
  filter(sp == "An_funestus")

p_funestus <- plot_gradient(
  grad_funestus_only,
  paste0("Spillover gradient (controls, An. funestus-only, ≤ ", max_dist_m, " m)")
)
print(p_funestus)


####################################
# Extract fixed effects summary 
####################################

# Fixed effects table from INLA (newer INLA stores it in summary.fixed)
fix_df <- as.data.frame(result_spillover$summary.fixed)

# Add term names as a column
fix_df$term <- rownames(result_spillover$summary.fixed)

# Reorder columns for readability
fix_df <- fix_df[, c("term", "mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode", "kld")]


###############################################################################
# Figure - Forest plot (IRR) of key fixed effects from the NB model
###############################################################################

# terms to show
terms_keep <- c(
  "treated_time",      # primary intervention effect (DiD)
  "tratPulv",          # baseline treated vs control
  "post",              # overall time effect
  "inOutOutdoor",      # outdoor vs indoor
  "spill_ctrl",        # spillover in controls
  "spill_ctrl_post"    # spillover change post
)

# Keep only terms that actually exist in the model
plot_df <- fix_df %>%
  filter(term %in% terms_keep) %>%
  mutate(
    IRR = exp(mean),
    IRR_low = exp(`0.025quant`),
    IRR_high = exp(`0.975quant`),
    term = factor(term, levels = rev(terms_keep))
  )

# Plot
ggplot(plot_df, aes(x = IRR, y = term)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_point() +
  geom_errorbarh(aes(xmin = IRR_low, xmax = IRR_high), height = 0.2) +
  scale_x_log10() +
  labs(
    x = "Incidence Rate Ratio (log scale)",
    y = NULL,
    title = "Key fixed effects (spillover-adjusted mosquito density model)"
  ) +
  theme_minimal()


###############################################################################
# Figure - Gradient plot — predicted density vs distance to opposite arm
###############################################################################

# Extract fitted mean mu for the estimation rows
idx_est <- inla.stack.index(stack_spill, tag = "est_spill")$data
mu_hat <- result_spillover$summary.fitted.values$mean[idx_est]

# Build a plotting dataset
grad_df <- agg_df_spill %>%
  mutate(
    mu_hat = mu_hat,
    # Rate per residency
    rate_hat = mu_hat / mean_NumRes
  ) %>%
    filter(tratPulv == 0) %>%
  # Guard against missing distances
  filter(!is.na(dist_to_opposite_m))

# Plot predicted rate vs distance (with a smooth)
ggplot(grad_df, aes(x = dist_to_opposite_m, y = rate_hat)) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = TRUE, method = "loess") +
  labs(
    x = "Distance to nearest treated location (m)",
    y = "Predicted mosquito density per residency",
    title = "Spillover gradient in controls: predicted density vs distance"
  ) +
  theme_minimal()

grad_df2 <- grad_df %>% filter(dist_to_opposite_m <= 500)

ggplot(grad_df2, aes(x = dist_to_opposite_m, y = rate_hat)) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = TRUE, method = "loess") +
  labs(
    x = "Distance to nearest treated location (m)",
    y = "Predicted mosquito density per residency",
    title = "Spillover gradient near the interface (≤ 500 m)"
  ) +
  theme_minimal()


###############################################################################
# Figure - Observed vs predicted counts (or per-residency rates)
###############################################################################

# Extract fitted mean mu again
idx_est <- inla.stack.index(stack_spill, tag = "est_spill")$data
mu_hat <- result_spillover$summary.fitted.values$mean[idx_est]

ppc_df <- agg_df_spill %>%
  mutate(
    mu_hat = mu_hat,
    obs = y,
    obs_rate = y / mean_NumRes,
    pred_rate = mu_hat / mean_NumRes
  )

# Observed vs predicted RATE per residency
ggplot(ppc_df, aes(x = pred_rate, y = obs_rate)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(
    x = "Predicted rate (per residency)",
    y = "Observed rate (per residency)",
    title = "Observed vs predicted mosquito density (per residency)"
  ) +
  theme_minimal()