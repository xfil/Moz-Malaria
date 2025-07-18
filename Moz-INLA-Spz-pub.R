# Clean CGPT

###############################################
## R Script for Spatial-Temporal Modeling
## Using INLA and Spatial Data (Mosquito Data)
## Author: Lopes LF
## Date: 14-vii-2025
## R version: Tested on 4.5.1
## INLA PAckage version: 25.6.7
###############################################

###############################################
## SECTION 1: Package Installation and Loading
###############################################

# Install necessary packages (run only once if not already installed)
install.packages(c(
  'tidyverse', 'ggplot2', 'readxl', 'spdep', 'terra', 'sf',
  'ggmap', 'leaflet', 'viridis', 'writexl', 'INLA',
  'lubridate', 'flextable', 'officer', 'rsatscan', 'tmap'
))

# Load required libraries
library(INLA)       # Bayesian spatial-temporal models
library(readxl)     # Excel data import
library(tidyverse)     # Data manipulation
library(ggplot2)    # Plotting
library(spdep)      # Spatial dependency analysis
library(sf)         # Spatial vector data handling
library(terra)      # Spatial raster/vector processing
library(ggmap)      # Map visualizations
library(leaflet)    # Interactive maps
library(viridis)    # Color scales
library(writexl)    # Excel data export
library(lubridate)  # Date/time manipulation
library(flextable)  # Tables for reporting
library(officer)    # Exporting reports
library(rsatscan)   # Spatial scan statistics (Kulldorff)
library(tmap)       # Thematic maps

# Check versions
packageVersion("INLA")
R.version.string

###############################################
## SECTION 2: Data Import and Preprocessing
###############################################

# Load dataset from Excel file
spz_file_path <- "BDSpz"
spzdf <- read_excel(spz_file_path)

# Check missing values
colSums(is.na(spzdf))

# Summarize missing values percentage by column
na_summary <- data.frame(
  Variable = names(spzdf),
  NA_Percentage = (colSums(is.na(spzdf)) / nrow(spzdf)) * 100
)
print(na_summary[na_summary$NA_Percentage > 0, ])

# Replace NA or empty values in quantCpMosq with 0
spzdf$quantCpMosq[is.na(spzdf$quantCpMosq) | spzdf$quantCpMosq == ""] <- 0
spzdf$quantCpMosq <- as.numeric(spzdf$quantCpMosq)

# Remove unidentified species (sp = 7) and An. tenebrosus (sp = 4)
spzdf <- spzdf[!(spzdf$sp %in% c(4, 7)), ]
spzdf$sp <- factor(spzdf$sp, levels = c(1, 2),
                   labels = c("An_funestus", "An_gambiae"))

# Convert coordinates and factors
spzdf <- spzdf %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    inOut = factor(inOut, levels = c(1, 2), labels = c("Indoor", "Outdoor"))
  ) %>%
  filter(!is.na(lon) & !is.na(lat))

###############################################
## SECTION 3: Spatial Data Preparation
###############################################

# Convert to sf object (WGS84) and transform to UTM CRS
spzdf_sf <- st_as_sf(spzdf, coords = c("lon", "lat"), crs = 4326)
spzdf_sf <- st_transform(spzdf_sf, crs = 32736)  # UTM zone 36S (Mozambique)

# Extract coordinates for modeling
coords <- st_coordinates(spzdf_sf)

# Create DiD variables
spzdf_sf <- spzdf_sf %>%
  mutate(
    post = ifelse(Ncol > 1, 1, 0),
    trat = tratPulv,
    treated_time = post * trat,
    Ncol = as.integer(Ncol)
  )

###############################################
## SECTION 4: SPDE Mesh Construction (INLA)
###############################################

# Construct spatial mesh (triangle edges 0.5–1km, cutoff at 100m)
mesh <- inla.mesh.2d(
  loc = coords,
  max.edge = c(500, 1000),
  cutoff = 100
)

# Plot mesh
plot(mesh)
points(coords, col = "red", cex = 0.5, pch = 16)

# Create spatial matrices for modeling
A <- inla.spde.make.A(mesh = mesh, loc = coords)

# Define SPDE model (prior parameters chosen based on domain knowledge)
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(10, 0.5),
  prior.sigma = c(1, 0.01)
)

# Spatial index for modeling
s.index <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

###############################################
## SECTION 5: INLA Model Implementation
###############################################

# Define model formula (base model DiD + RW1 temporal + spatial effect)
formula <- pfSpz ~ 0 + intercept + post + trat + treated_time +
  f(Ncol, model = "rw1") +
  f(spatial.field, model = spde)

# Create data stack for INLA
stack <- inla.stack(
  data = list(pfSpz = spzdf_sf$pfSpz),
  A = list(A, 1),
  effects = list(
    c(s.index, list(intercept = 1)),
    list(post = spzdf_sf$post,
         trat = spzdf_sf$trat,
         treated_time = spzdf_sf$treated_time,
         Ncol = spzdf_sf$Ncol)
  ),
  tag = "base_model"
)

# Run INLA model
result <- inla(
  formula,
  family = "binomial",
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Summarize base model results
summary(result)

###############################################
## Zero Inflated Model
###############################################

result_zib <- inla(
  formula,
  family = "zeroinflatedbinomial1",
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

summary(result_zib)







###############################################
## SECTION 6: Advanced Models and Interactions
###############################################

# Example: Model adding mosquito species (sp) and interaction with treatment-time
formula_sp_interaction <- pfSpz ~ 0 + intercept + post + trat + treated_time +
  sp + treated_time:sp +
  f(Ncol, model = "rw1") +
  f(spatial.field, model = spde)

stack_sp <- inla.stack(
  data = list(pfSpz = spzdf_sf$pfSpz),
  A = list(A, 1),
  effects = list(
    c(s.index, list(intercept = 1)),
    list(
      post = spzdf_sf$post,
      trat = spzdf_sf$trat,
      treated_time = spzdf_sf$treated_time,
      sp = spzdf_sf$sp,
      Ncol = spzdf_sf$Ncol
    )
  ),
  tag = "sp_interaction"
)

# Run the species interaction model
result_sp <- inla(
  formula_sp_interaction,
  family = "binomial",
  data = inla.stack.data(stack_sp),
  control.predictor = list(A = inla.stack.A(stack_sp), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

summary(result_sp)

###############################################
## SECTION 7: Results Export & Visualization
###############################################

# Example: Save summary to Excel file
write_xlsx(as.data.frame(result_sp$summary.fixed),
           path = "model_sp_interaction_results.xlsx")

# Example: Interactive spatial visualization with leaflet (WGS84)
spzdf_leaflet <- st_transform(spzdf_sf, 4326)
leaflet(spzdf_leaflet) %>%
  addTiles() %>%
  addCircleMarkers(radius = 4, color = ~viridis_pal()(pfSpz + 1))

###############################################
## END OF SCRIPT
###############################################
