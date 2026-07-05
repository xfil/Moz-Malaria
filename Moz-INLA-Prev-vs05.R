###############################################################################
# Prevalence geospatial analysis script
#
# Study:
#   Field evaluation of 3-(N-acetyl-n-butyl) aminopropionic acid ethyl ester
#   (IR3535) as a spatial repellent to control malaria:
#   A Randomised, Before-After-Control-Intervention trial
#
# Dataset:
#   Malaria prevalence survey data collected at three timepoints:
#     P1 = baseline / pre-intervention
#     P2 = post-intervention
#     P3 = post-intervention
#
# Main objective:
#   Fit sequential Bayesian geospatial prevalence models using R-INLA to assess
#   malaria infection status in relation to treatment arm, post-intervention
#   period, treatment-by-time interaction, individual-level covariates, temporal
#   structure, spatial autocorrelation, and an optional treated/control boundary
#   spillover covariate.
#
# Main modelling structure:
#   - Outcome: Malaria infection status (Mal; 0/1)
#   - Likelihood: Binomial
#   - Fixed effects:
#       tratPulv, post, treated_time
#       plus sequential covariates: sex, idade, UltMal_grouped, CriancasAF
#   - Random effects:
#       temporal RW1 over survey timepoint
#       spatial SPDE field
#   - Spillover model:
#       distance-based exponential spillover term for control observations near
#       treated observations, including post-intervention interaction
###############################################################################

# Install packages
# install.packages('dplyr')
# install.packages('ggplot2')
# install.packages('readxl')
# install.packages('spData')
# install.packages('sn')
# install.packages('foreach')
# install.packages('terra')
# install.packages('sf')
# install.packages('spdep')
# install.packages('stringr')
# install.packages('tidyr')
# install.packages('reshape2')
# install.packages("ggmap")
# install.packages("leaflet")
# install.packages("raster")
# install.packages("viridis")
# install.packages("writexl")
# install.packages("fmesher")
# install.packages("INLA")
# install.packages("spatstat")
# install.packages("spatstat.core")
# install.packages("spatstat.geom")
# install.packages("spatstat.data")
# install.packages("spatstat.explore")
# install.packages("lubridate")
# install.packages("flextable")
# install.packages("officer")
# install.packages("rsatscan")
# install.packages("tmap")
# install.packages("ggmap")
# install.packages("htmlwidgets")

# Load libraries
library(INLA)
library(readxl)
library(stringr)
library(tidyr)
library(dplyr)
library(spdep)
library(sn)
library(ggplot2)
library(sp)
library(sf)
library(terra)
library(fmesher)
library(reshape2)
library(ggmap)
library(leaflet)
library(raster)
library(viridis)
library(writexl)
library(spatstat.geom)
library(spatstat.core)
library(spatstat.explore)
library(lubridate)
library(flextable)
library(officer)
library(gganimate)
library(rsatscan)
library(tmap)
library(ggmap)
library(htmlwidgets)
library(scales)

# packageVersion("INLA")


#########################################
# Import data
#########################################

# Load data  
prev_file_path <- "BDPrev.xlsx"

# Read the first sheet of the file
prevdf <- read_excel(prev_file_path)


#########################################
# Check data
#########################################

# check missing data
sum(is.na(prevdf))

# check missing data by column
colSums(is.na(prevdf))

# The complete.cases() function returns a logical vector indicating which rows have no NA values (TRUE) and which have at least one (FALSE).
complete.cases(prevdf)

#Percentage of NAs per column:
col_na_counts <- colSums(is.na(prevdf))
col_na_percentages <- (col_na_counts / nrow(prevdf)) * 100
print(col_na_percentages)

# Create a data frame with two columns: "Variable" and "NA_Percentage"
na_summary_df <- data.frame(
  Variable = names(col_na_percentages),  # Column names of original df
  NA_Percentage = as.numeric(col_na_percentages) # The calculated percentages
)
print(na_summary_df)

# Filter the data frame to show only rows where NA_Percentage is not zero
na_summary_filtered_df <- na_summary_df[na_summary_df$NA_Percentage != 0, ]

print("\nFiltered NA Summary (NA_Percentage != 0):")
print(na_summary_filtered_df)


#########################################
# Transform data
#########################################

prevdf_long <- prevdf %>%
  mutate(across(matches("^P[1-3]"), as.character)) %>%  # force all repeated vars to character
  pivot_longer(
    cols = matches("^P[1-3]"),
    names_to = c("timepoint", "variable"),
    names_pattern = "(P[1-3])(.*)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = "variable",
    values_from = "value"
  )

#Inspect non-numeric values
unique(prevdf_long$Idade[!grepl("^[0-9.]+$", prevdf_long$Idade)])

prevdflg <- prevdf_long %>%
  mutate(
    timepoint = factor(timepoint, levels = c("P1", "P2", "P3")),
    # Standardize and convert malaria outcome
    Mal = str_trim(Mal),
    Mal = na_if(Mal, ""),                  # Convert blanks to NA
    Mal = factor(Mal, levels = c("0", "1")),  # "0" = negative, "1" = positive
    
    Idade = str_trim(Idade),
    Idade = na_if(Idade, ""),
    Idade = as.numeric(Idade),
    CriancasAF = as.numeric(CriancasAF),
    sex = factor(sex),
    Afnovo <- factor(Afnovo),
    UltimaVezMalaria = as.character(UltimaVezMalaria)
  )


prevdflg <- prevdflg %>%
  # Create a unique ID based on 'ordem5abr25' and 'timepoint'
  mutate(
    unique_id = paste0(ordem5abr25, "_", timepoint),
    
    # Add year and month columns based on timepoint
    year = case_when(
      timepoint == "P1" ~ 2021,
      timepoint == "P2" ~ 2022,
      timepoint == "P3" ~ 2023
    ),
    
    month = case_when(
      timepoint == "P1" ~ 6,  # June
      timepoint == "P2" ~ 6,  # June
      timepoint == "P3" ~ 7   # July
    )
  )

#df <- df %>%
#  rename(AFnovo = afnovo)

prevdflg <- rename(prevdflg, AFnovo = Afnovo)
prevdflg <- rename(prevdflg, idade = Idade)

table(prevdflg$Mal, useNA = "always")


prevdfclean <- prevdflg %>% 
  filter(
    !is.na(Mal),
  )


#########################################
# Preprocess data
#########################################

#Percentage of NAs per column:
col_na_counts <- colSums(is.na(prevdfclean))
col_na_percentages <- (col_na_counts / nrow(prevdfclean)) * 100
print(col_na_percentages)

# Create a data frame with two columns: "Variable" and "NA_Percentage"
na_summary_df <- data.frame(
  Variable = names(col_na_percentages),  # Column names of original df
  NA_Percentage = as.numeric(col_na_percentages) # The calculated percentages
)
print(na_summary_df)

# Filter the data frame to show only rows where NA_Percentage is not zero
na_summary_filtered_df <- na_summary_df[na_summary_df$NA_Percentage != 0, ]

print("\nFiltered NA Summary (NA_Percentage != 0):")
print(na_summary_filtered_df)

# Create post indicator (0 = P1, 1 = P2 or P3)
prevdfclean$post <- ifelse(prevdfclean$timepoint == "P1", 0, 1)

# Create numeric time_index for rw1 (e.g., 1 = P1, 2 = P2, 3 = P3)
prevdfclean$time_index <- as.numeric(factor(prevdfclean$timepoint, levels = c("P1", "P2", "P3")))

# Ensure tratPulv is numeric/factor as needed
prevdfclean$tratPulv <- as.numeric(prevdfclean$tratPulv)

# Create interaction term manually 
prevdfclean$treated_time <- prevdfclean$tratPulv * prevdfclean$post

# Clean or impute CriancasAF (e.g., replace NA with median)
prevdfclean$CriancasAF[is.na(prevdfclean$CriancasAF)] <- median(prevdfclean$CriancasAF, na.rm = TRUE)

prevdfclean$Mal <- as.numeric(as.character(prevdfclean$Mal))

prevdfclean$sex[prevdfclean$timepoint == "P3" & prevdfclean$sex == 1] <- 0
prevdfclean$sex[prevdfclean$timepoint == "P3" & prevdfclean$sex == 2] <- 1

prevdfclean$sex <- as.numeric(as.character(prevdfclean$sex))

str(prevdfclean$idade)
summary(prevdfclean$idade)
sum(is.na(prevdfclean$idade))

# Drop data with missing values age
prevdfclean <- prevdfclean[!is.na(prevdfclean$idade), ]

prevdfclean$lon <- as.numeric(prevdfclean$lon)
prevdfclean$lat <- as.numeric(prevdfclean$lat)

summary(prevdfclean$UltimaVezMalaria)
table(prevdfclean$UltimaVezMalaria, useNA = "ifany")

# 2 entries with '99'
prevdfclean <- prevdfclean[prevdfclean$UltimaVezMalaria != 99, ]

# Create a new grouped variable
prevdfclean$UltMal_grouped <- NA

# Recode based on original values
prevdfclean$UltMal_grouped[prevdfclean$UltimaVezMalaria %in% c(1, 2)] <- 1  # recent
prevdfclean$UltMal_grouped[prevdfclean$UltimaVezMalaria == 3] <- 2         # past
prevdfclean$UltMal_grouped[prevdfclean$UltimaVezMalaria == 4] <- 3         # unknown
prevdfclean$UltMal_grouped[prevdfclean$UltimaVezMalaria == 5] <- 4         # never

# Convert to factor with labels
prevdfclean$UltMal_grouped <- factor(prevdfclean$UltMal_grouped,
                                     levels = c(1, 2, 3, 4),
                                     labels = c("recent", "past", "unknown", "never")
)

summary(prevdfclean$CriancasAF)
sum(is.na(prevdfclean$CriancasAF))

prevdfclean <- prevdfclean[prevdfclean$CriancasAF != 0, ]


#########################################
# Model 001
#########################################

# Create spatial coordinates matrix for SPDE
coords <- as.matrix(prevdfclean[, c("lon", "lat")])

# Build the SPDE mesh
mesh <- inla.mesh.2d(loc = coords, max.edge = c(0.02, 0.2), cutoff = 0.01)

# Define the SPDE model
spde <- inla.spde2.pcmatern(mesh = mesh,
                            alpha = 2,
                            prior.range = c(0.1, 0.01),
                            prior.sigma = c(1, 0.01))

# Link observations to the mesh
A <- inla.spde.make.A(mesh = mesh, loc = coords)

# Create spatial index
spatial_index <- inla.spde.make.index(name = "spatial_field", n.spde = spde$n.spde)


# Create INLA stack
stack <- inla.stack(
  data = list(Mal = prevdfclean$Mal),
  A = list(A, 1),
  effects = list(
    spatial_field = spatial_index,
    data.frame(
      intercept = 1,
      tratPulv = prevdfclean$tratPulv,
      post = prevdfclean$post,
      treated_time = prevdfclean$treated_time,
      time_index = prevdfclean$time_index
    )
  ),
  tag = "est"
)

# Model 1: Treatment + Time + Interaction + Spatiotemporal Random Effects
formula1 <- Mal ~ 1 + tratPulv + post + treated_time +
  f(time_index, model = "rw1") +
  f(spatial_field, model = spde)

# Fit the model
result1 <- inla(
  formula1,
  family = "binomial",
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose = FALSE
)

summary(result1)


#########################################
# Model 2 - add sex
#########################################

# check if sex is as factor
table(prevdfclean$sex)

# stack model 2
stack2 <- inla.stack(
  data = list(Mal = prevdfclean$Mal),
  A = list(A, 1),
  effects = list(
    spatial_field = spatial_index,
    data.frame(
      intercept = 1,
      tratPulv = prevdfclean$tratPulv,
      post = prevdfclean$post,
      treated_time = prevdfclean$treated_time,
      sex = prevdfclean$sex,
      time_index = prevdfclean$time_index
    )
  ),
  tag = "est"
)

# Model formula
formula2 <- Mal ~ 1 + tratPulv + post + treated_time + sex +
  f(time_index, model = "rw1") +
  f(spatial_field, model = spde)

# Fit model
result2 <- inla(
  formula2,
  family = "binomial",
  data = inla.stack.data(stack2),
  control.predictor = list(A = inla.stack.A(stack2), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose = FALSE
)

summary(result2)


#########################################
# Model 3 - add age 
#########################################

# Model 3 formula
formula3 <- Mal ~ 1 + tratPulv + post + treated_time + sex + idade +
  f(time_index, model = "rw1") +
  f(spatial_field, model = spde)

# Stack for model 3
stack3 <- inla.stack(
  data = list(Mal = prevdfclean$Mal),
  A = list(A, 1),
  effects = list(
    spatial_field = spatial_index,
    data.frame(
      intercept = 1,
      tratPulv = prevdfclean$tratPulv,
      post = prevdfclean$post,
      treated_time = prevdfclean$treated_time,
      sex = prevdfclean$sex,
      idade = prevdfclean$idade,
      time_index = prevdfclean$time_index
    )
  ),
  tag = "est"
)

# Fit model 3
result3 <- inla(
  formula3,
  family = "binomial",
  data = inla.stack.data(stack3),
  control.predictor = list(A = inla.stack.A(stack3), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose = FALSE
)

summary(result3)


#########################################
# Model 4 - add last time had malaria
#########################################

# Rebuild stack for Model 4 
stack4 <- inla.stack(
  data = list(Mal = prevdfclean$Mal),
  A = list(A, 1),
  effects = list(
    spatial_field = spatial_index,
    data.frame(
      intercept = 1,
      tratPulv = prevdfclean$tratPulv,
      post = prevdfclean$post,
      treated_time = prevdfclean$treated_time,
      sex = prevdfclean$sex,
      idade = prevdfclean$idade,
      UltMal_grouped = prevdfclean$UltMal_grouped,
      time_index = prevdfclean$time_index
    )
  ),
  tag = "est"
)

# Define formula for Model 4
formula4 <- Mal ~ 1 + tratPulv + post + treated_time + sex + idade + UltMal_grouped +
  f(time_index, model = "rw1") +
  f(spatial_field, model = spde)

# Fit Model 4
result4 <- inla(
  formula4,
  family = "binomial",
  data = inla.stack.data(stack4),
  control.predictor = list(A = inla.stack.A(stack4), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose = FALSE
)

summary(result4)


#########################################
# Model 5 - build from model 3 add number of children AF
#########################################

# Build Model 5 stack
stack5 <- inla.stack(
  data = list(Mal = prevdfclean$Mal),
  A = list(A, 1),
  effects = list(
    spatial_field = spatial_index,
    data.frame(
      intercept = 1,
      tratPulv = prevdfclean$tratPulv,
      post = prevdfclean$post,
      treated_time = prevdfclean$treated_time,
      sex = prevdfclean$sex,
      idade = prevdfclean$idade,
      CriancasAF = prevdfclean$CriancasAF,
      time_index = prevdfclean$time_index
    )
  ),
  tag = "est"
)

# Define model 5
formula5 <- Mal ~ 1 + tratPulv + post + treated_time + sex + idade + CriancasAF +
  f(time_index, model = "rw1") +
  f(spatial_field, model = spde)

# Fit model 5
result5 <- inla(
  formula5,
  family = "binomial",
  data = inla.stack.data(stack5),
  control.predictor = list(A = inla.stack.A(stack5), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose = FALSE
)

summary(result5)


################################################################################
# Model 6 — Model 5 + spillover/gradient covariate (treated/control interface)
#
# Goal:
#   Extend Model 5 (malaria prevalence) by adding a distance-based spillover term
#   that captures potential interference/contamination near the treated/control
#   boundary (no buffer area).
#
#   1) Compute, for each observation, the distance (meters) to the nearest
#      household in the *opposite* treatment arm.
#   2) Convert distance into a smooth weight w = exp(-d / r_m).
#   3) Define a spillover exposure variable for CONTROLS only:
#        spill_ctrl = w if tratPulv==0 else 0
#   4) Optionally allow spillover only in post period:
#        spill_ctrl_post = spill_ctrl * post
#
#   - Distances are computed in a projected CRS (UTM) so r_m is in meters.
################################################################################

# Choose the spillover decay scale (meters)
r_m <- 250

# If TRUE, include spill_ctrl_post (spillover active mainly post-intervention).
# If FALSE, include only spill_ctrl (spillover regardless of period).
include_spillover_post <- TRUE

# Ensure coordinates and treatment indicators are present
stopifnot(all(c("lon", "lat", "tratPulv", "post") %in% names(prevdfclean)))

# Build sf object (WGS84 lon/lat)
prev_sf_ll <- st_as_sf(
  prevdfclean,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE  # keep lon/lat columns
)

# Select an appropriate UTM CRS automatically (meters)
lon_med <- stats::median(prevdfclean$lon, na.rm = TRUE)
lat_med <- stats::median(prevdfclean$lat, na.rm = TRUE)
utm_zone <- floor((lon_med + 180) / 6) + 1
epsg_utm <- ifelse(lat_med < 0, 32700 + utm_zone, 32600 + utm_zone)

# Transform to UTM CRS so st_distance returns meters
prev_sf_m <- st_transform(prev_sf_ll, crs = epsg_utm)

# Split treated vs control and compute nearest-opposite distances
sf_treat_m <- prev_sf_m[prev_sf_m$tratPulv == 1, ]
sf_ctrl_m  <- prev_sf_m[prev_sf_m$tratPulv == 0, ]

# both groups must exist
if (nrow(sf_treat_m) == 0 || nrow(sf_ctrl_m) == 0) {
  stop("Spillover gradient requires BOTH treated (tratPulv==1) and control (tratPulv==0) observations.")
}

# For each treated obs, find nearest control obs index
idx_near_ctrl_for_treat <- st_nearest_feature(sf_treat_m, sf_ctrl_m)

# Distances treated -> nearest control (meters)
d_treat_to_ctrl_m <- as.numeric(
  st_distance(sf_treat_m, sf_ctrl_m[idx_near_ctrl_for_treat, ], by_element = TRUE)
)

# For each control obs, find nearest treated obs index
idx_near_treat_for_ctrl <- st_nearest_feature(sf_ctrl_m, sf_treat_m)

# Distances control -> nearest treated (meters)
d_ctrl_to_treat_m <- as.numeric(
  st_distance(sf_ctrl_m, sf_treat_m[idx_near_treat_for_ctrl, ], by_element = TRUE)
)

# Assemble a full-length vector aligned with prevdfclean row order
dist_to_opposite_m <- rep(NA_real_, nrow(prevdfclean))
dist_to_opposite_m[prevdfclean$tratPulv == 1] <- d_treat_to_ctrl_m
dist_to_opposite_m[prevdfclean$tratPulv == 0] <- d_ctrl_to_treat_m

# Convert distance to a smooth decay weight
# Exponential decay: w = exp(-d/r_m)
w_boundary <- exp(-dist_to_opposite_m / r_m)

# Define spillover covariates (controls only) 
# Controls near treated areas may be partially "exposed", so they get positive
# spillover values. Treated observations are set to 0 for this contamination term.
spill_ctrl <- ifelse(prevdfclean$tratPulv == 0, w_boundary, 0)

# spillover active only after intervention
spill_ctrl_post <- spill_ctrl * prevdfclean$post

# Add to prevdfclean for tracking/documentation
prevdfclean$dist_to_opposite_m <- dist_to_opposite_m
prevdfclean$w_boundary         <- w_boundary
prevdfclean$spill_ctrl         <- spill_ctrl
prevdfclean$spill_ctrl_post    <- spill_ctrl_post

# Rebuild coordinates matrix for A if needed 
coords_prev <- as.matrix(prevdfclean[, c("lon", "lat")])

# Rebuild A with the same mesh so models remain comparable
A_prev <- inla.spde.make.A(mesh = mesh, loc = coords_prev)

# Build Model 6 stack (Model 5 + spillover)
# the same structure as Model 5, adding spill_ctrl (+ spill_ctrl_post).
stack6 <- inla.stack(
  data = list(Mal = prevdfclean$Mal),
  A = list(A_prev, 1),
  effects = list(
    spatial_field = spatial_index,
    data.frame(
      intercept    = 1,
      tratPulv     = prevdfclean$tratPulv,
      post         = prevdfclean$post,
      treated_time = prevdfclean$treated_time,
      sex          = prevdfclean$sex,
      idade        = prevdfclean$idade,
      CriancasAF   = prevdfclean$CriancasAF,
      spill_ctrl   = prevdfclean$spill_ctrl,
      spill_ctrl_post = prevdfclean$spill_ctrl_post,
      time_index   = prevdfclean$time_index
    )
  ),
  tag = "est"
)

# Define Model 6 formula
if (include_spillover_post) {
  formula6 <- Mal ~ 1 + tratPulv + post + treated_time + sex + idade + CriancasAF +
    spill_ctrl + spill_ctrl_post +
    f(time_index, model = "rw1") +
    f(spatial_field, model = spde)
} else {
  formula6 <- Mal ~ 1 + tratPulv + post + treated_time + sex + idade + CriancasAF +
    spill_ctrl +
    f(time_index, model = "rw1") +
    f(spatial_field, model = spde)
}

# Fit Model 6
result6_spillover <- inla(
  formula6,
  family = "binomial",
  data = inla.stack.data(stack6),
  control.predictor = list(A = inla.stack.A(stack6), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose = FALSE
)

# Summarize Model 6
summary(result6_spillover)

# Compare Model 5 vs Model 6
cat("\n--- Model 5 vs Model 6 (spillover) comparison ---\n")
cat("Model 5 DIC :", result5$dic$dic, "\n")
cat("Model 6 DIC :", result6_spillover$dic$dic, "\n")
cat("Model 5 WAIC:", result5$waic$waic, "\n")
cat("Model 6 WAIC:", result6_spillover$waic$waic, "\n")

# CPO failures should ideally be 0
cat("\nCPO failures:\n")
cat("Model 5 failures:", sum(result5$cpo$failure, na.rm = TRUE), "\n")
cat("Model 6 failures:", sum(result6_spillover$cpo$failure, na.rm = TRUE), "\n")

# Optional: print the spillover coefficient(s) directly
if ("spill_ctrl" %in% rownames(result6_spillover$summary.fixed)) {
  cat("\n--- Spillover coefficient: spill_ctrl ---\n")
  print(result6_spillover$summary.fixed["spill_ctrl", ])
}
if ("spill_ctrl_post" %in% rownames(result6_spillover$summary.fixed)) {
  cat("\n--- Spillover coefficient: spill_ctrl_post ---\n")
  print(result6_spillover$summary.fixed["spill_ctrl_post", ])
}

################################################################################
# End of Model 6
################################################################################


################################################################################
# PREVALENCE MODEL 6 (SPILLOVER) - POSTERIOR, FOREST PLOT, MAPS
# Uses: result6_spillover, stack6, prevdfclean, study_areas, mesh, spde
################################################################################

# fail checks 
stopifnot(exists("result6_spillover"))      # fitted model exists
stopifnot(exists("stack6"))                 # estimation stack exists
stopifnot(exists("prevdfclean"))            # prevalence data exists
stopifnot(all(c("lon","lat") %in% names(prevdfclean)))   # coords exist
stopifnot(exists("study_areas"))            # polygons exist
stopifnot(exists("mesh"))                   # mesh exists
stopifnot(exists("spde"))                   # spde exists


# Posterior distribution of treated_time (and OR version)

# Extract marginal posterior for treated_time (log-odds scale)
m_treat <- result6_spillover$marginals.fixed[["treated_time"]]

# Stop with clear message if term name is different
if (is.null(m_treat)) {
  stop(
    "Could not find 'treated_time' in result6_spillover$marginals.fixed.\n",
    "Run: rownames(result6_spillover$summary.fixed) to see the exact term name."
  )
}

# Convert marginal to a data.frame for plotting 
df_post <- as.data.frame(m_treat)                         # columns: x, y
colnames(df_post) <- c("log_odds", "density")             # rename columns

#Compute posterior mean and 95% credible interval 
post_mean <- INLA::inla.emarginal(function(x) x, m_treat) # E[beta]
post_q025 <- INLA::inla.qmarginal(0.025, m_treat)         # 2.5%
post_q975 <- INLA::inla.qmarginal(0.975, m_treat)         # 97.5%

# Plot posterior density (log-odds)
p_post <- ggplot(df_post, aes(x = log_odds, y = density)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +       # no-effect line
  geom_vline(xintercept = post_mean, linetype = "solid") +
  labs(
    title = "Posterior distribution: treated_time (prevalence Model 6)",
    subtitle = paste0(
      "Mean = ", round(post_mean, 3),
      " | 95% CrI [", round(post_q025, 3), ", ", round(post_q975, 3), "]"
    ),
    x = "Effect on log-odds",
    y = "Posterior density"
  ) +
  theme_minimal()

print(p_post)

# Plot posterior of Odds Ratio (OR = exp(beta)) 
set.seed(123)                                              # reproducible sampling
beta_samp <- INLA::inla.rmarginal(20000, m_treat)           # draw samples from marginal
or_samp   <- exp(beta_samp)                                 # convert to odds ratios

df_or <- data.frame(or = or_samp)

p_or <- ggplot(df_or, aes(x = or)) +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed") +        # OR=1 no effect
  labs(
    title = "Posterior distribution: Odds Ratio for treated_time (Model 6)",
    x = "Odds Ratio (exp(treated_time))",
    y = "Posterior density"
  ) +
  theme_minimal()

print(p_or)


# Forest plot of fixed effects (Odds Ratios)

# Extract fixed effect summaries 
fx <- as.data.frame(result6_spillover$summary.fixed)        # table with mean and CrIs
fx$term <- rownames(fx)                                     # keep coefficient names

# Convert to Odds Ratios (logit link)
fx <- fx %>%
  mutate(
    OR      = exp(mean),                                    # mean OR
    OR_low  = exp(`0.025quant`),                             # lower CrI OR
    OR_high = exp(`0.975quant`)                              # upper CrI OR
  )

# Remove intercept
fx_plot <- fx %>%
  filter(term != "(Intercept)") %>%                         # common name
  filter(term != "intercept")                               # just in case

# Order terms by OR
fx_plot <- fx_plot %>%
  arrange(OR) %>%
  mutate(term = factor(term, levels = term))

# Forest plot
p_forest <- ggplot(fx_plot, aes(x = OR, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed") +        # no-effect OR line
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2) +
  geom_point(size = 2) +
  scale_x_continuous(trans = "log10", labels = label_number()) +
  labs(
    title = "Fixed effects (prevalence Model 6): Posterior Odds Ratios",
    x = "Odds Ratio (log scale)",
    y = NULL
  ) +
  theme_minimal()

print(p_forest)

################################################################################
# Predicted risk maps (points)
################################################################################

# Get estimation indices from the stack (tag = 'est')
idx_est <- INLA::inla.stack.index(stack6, tag = "est")$data

# Extract fitted values (posterior mean + CrI)
risk_mean <- result6_spillover$summary.fitted.values$mean[idx_est]
risk_q025 <- result6_spillover$summary.fitted.values$`0.025quant`[idx_est]
risk_q975 <- result6_spillover$summary.fitted.values$`0.975quant`[idx_est]

# Attach predictions back to the same data rows
prev_map_sf <- prevdfclean %>%
  mutate(
    risk_mean = risk_mean,                                  # predicted prevalence (mean)
    risk_q025 = risk_q025,                                  # 2.5% CrI
    risk_q975 = risk_q975,                                  # 97.5% CrI
    risk_ci_width = risk_q975 - risk_q025                   # uncertainty width
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Point map: posterior mean prevalence
tmap_mode("plot")

m_prev_mean <- tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill.scale = tm_scale(values = c("Control" = "skyblue", "Repellent" = "orange")),
    fill.legend = tm_legend(title = "Study area")
  ) +
  tm_shape(prev_map_sf) +
  tm_symbols(
    fill = "risk_mean",
    size = 0.55,
    fill.scale = tm_scale_intervals(style = "quantile", values = "viridis"),
    fill.legend = tm_legend(title = "Predicted prevalence\n(posterior mean)")
  ) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "right",
            legend.text.size = 0.9,
            legend.title.size = 1.0,
            legend.width = 10) +
  tm_title("Predicted malaria prevalence")

m_prev_mean

# Point map: uncertainty (CrI width)
m_prev_unc <- tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill.scale = tm_scale(values = c("Control" = "skyblue", "Repellent" = "orange")),
    fill.legend = tm_legend(title = "Study area")
  ) +
  tm_shape(prev_map_sf) +
  tm_symbols(
    fill = "risk_ci_width",
    size = 0.55,
    fill.scale = tm_scale_intervals(style = "quantile", values = "viridis"),
    fill.legend = tm_legend(title = "Uncertainty\n(95% CrI width)")
  ) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "right",
            legend.text.size = 0.9,
            legend.title.size = 1.0,
            legend.width = 10) +
  tm_title("Uncertainty in predicted prevalence")

m_prev_unc


tmap_mode("plot")

m_prev_mean <- tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill_alpha = 0.25,
    fill.scale  = tm_scale(values = c("Control" = "skyblue", "Repellent" = "orange")),
    fill.legend = tm_legend(title = "Study area")
  ) +

  tm_shape(prev_map_sf) +
  tm_symbols(
    fill = "risk_mean",
    size = 0.55,
    fill.scale  = tm_scale_intervals(style = "quantile", values = "viridis"),
    fill.legend = tm_legend(title = "Predicted prevalence\n(posterior mean)")
  ) +

  # show household/sampling locations on top 
  tm_shape(hh_sf) +
  tm_dots(col = "black", fill = "white", size = 0.3, alpha = 0.8) +

  # Borders LAST so they overlay everything
  tm_shape(study_areas) +
  tm_borders(col = "black", lwd = 1.5) +

  # Axes labels outside, no grid lines
  tm_grid(
    labels.show = TRUE,
    labels.inside.frame = FALSE,
    lines = FALSE,
    x = lon_ticks,  
    y = lat_ticks,
    labels.size = 0.7,
    labels.format = list(digits = 3)
  ) +

  tm_title("Predicted malaria prevalence (posterior mean)") +

  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.width = 10,
    legend.text.size = 0.9,
    legend.title.size = 1.0,
    outer.margins = c(0.05, 0.05, 0.02, 0.02),
    clip = TRUE
  )

m_prev_mean


tmap_mode("plot")

m_prev_mean <- tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill_alpha = 0.25,
    fill.scale  = tm_scale(values = c("Control" = "skyblue", "Repellent" = "orange")),
    fill.legend = tm_legend(title = "Study area")
  ) +
  
  tm_shape(prev_map_sf) +
  tm_symbols(
    fill = "risk_mean",
    size = 0.55,
    fill.scale  = tm_scale_intervals(style = "quantile", values = "viridis"),
    fill.legend = tm_legend(title = "Predicted prevalence\n(posterior mean)")
  ) +
  
  # show household/sampling locations on top
  tm_shape(hh_sf) +
  tm_dots(col = "black", fill = "white", size = 0.3, alpha = 0.8) +
  
  # Borders LAST so they overlay everything
  tm_shape(study_areas) +
  tm_borders(col = "black", lwd = 1.5) +
  
  # Axes labels outside, no grid lines 
  tm_grid(
    labels.show = TRUE,
    labels.inside.frame = FALSE,
    lines = FALSE,
    x = lon_ticks,   
    y = lat_ticks,
    labels.size = 0.7,
    labels.format = list(digits = 3)
  ) +
  
  tm_title("Predicted malaria prevalence (posterior mean)") +
  
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.width = 10,
    legend.text.size = 0.9,
    legend.title.size = 1.0,
    outer.margins = c(0.05, 0.05, 0.02, 0.02),
    clip = TRUE
  )

m_prev_mean


################################################################################
# “smooth” surface map (spatial random effect)
################################################################################

# This produces a smooth surface without terra/stars masking issues.
# It maps the *spatial random effect* only

# Extract spatial random field summary 
sp_field <- result6_spillover$summary.random[["spatial_field"]]

# Fail check
if (is.null(sp_field)) {
  stop(
    "Could not find random effect 'spatial_field' in result6_spillover$summary.random.\n",
    "Run: names(result6_spillover$summary.random) to see the exact name."
  )
}

# Build a projector grid that covers the study area
bb <- sf::st_bbox(sf::st_union(study_areas) %>% sf::st_transform(4326))

proj <- INLA::inla.mesh.projector(
  mesh = mesh,
  xlim = c(as.numeric(bb["xmin"]), as.numeric(bb["xmax"])),
  ylim = c(as.numeric(bb["ymin"]), as.numeric(bb["ymax"])),
  dims = c(250, 250)                                         # increase for smoother surface
)

# Project posterior mean spatial field onto the grid
field_mean_vec <- sp_field$mean                              # mean at mesh vertices
field_grid <- INLA::inla.mesh.project(proj, field_mean_vec)  # grid matrix

# Convert projected grid to a tidy data.frame
df_field <- expand.grid(lon = proj$x, lat = proj$y)
df_field$field_mean <- as.vector(field_grid)

# Clip to study polygon
df_field_sf <- sf::st_as_sf(df_field, coords = c("lon","lat"), crs = 4326, remove = FALSE)
study_union <- sf::st_union(study_areas) %>% sf::st_transform(4326)

inside <- sf::st_within(df_field_sf, study_union, sparse = FALSE)[, 1]
df_field <- df_field[inside, , drop = FALSE]

# Smooth-looking spatial random effect surface
# Interpretation:
# - Positive field_mean => higher prevalence than expected after covariates
# - Negative field_mean => lower prevalence than expected after covariates
p_field <- ggplot(df_field, aes(x = lon, y = lat, fill = field_mean)) +
  geom_raster(interpolate = TRUE) +
  geom_sf(data = study_areas, inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.4) +
  coord_sf(expand = FALSE) +
  labs(
    title = "Spatial random effect surface (prevalence Model 6)",
    subtitle = "Posterior mean of the SPDE field (log-odds scale)",
    fill = "Spatial effect\n(mean)"
  ) +
  theme_minimal()

print(p_field)