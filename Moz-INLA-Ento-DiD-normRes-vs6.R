library(INLA)
library(dplyr)
library(readxl)
library(sf)
library(sp)

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

# Create spatial mesh
coords <- st_as_sf(agg_df, coords = c("lon", "lat"), crs = 4326)
mesh <- inla.mesh.2d(loc = coords, max.edge = c(0.02, 0.2), cutoff = 0.01)

# SPDE and A matrix
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(0.02, 0.01),
                            prior.sigma = c(1, 0.01))
s_index <- inla.spde.make.index("space", spde$n.spde)

# Convert to spatial points for A matrix

# Extract numeric coordinates from sf object
coords_xy <- st_coordinates(coords)
agg_df$lon <- coords_xy[, "X"]
agg_df$lat <- coords_xy[, "Y"]

# Create agg_sp with numeric lon/lat
agg_sp <- agg_df
coordinates(agg_sp) <- ~lon + lat
proj4string(agg_sp) <- CRS("+proj=longlat +datum=WGS84")


A <- inla.spde.make.A(mesh = mesh, loc = coordinates(agg_sp))

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


##########################################################
### MAPS
##########################################################


library(sf)
library(tmap)
library(INLA)
library(dplyr)

# ----------------------------
# 1. Load Study Area Shapefiles
# ----------------------------
control_area_sf <- st_read("D:/OneDrive - IHMT-NOVA/IHMT/Projectos/Moçambique-HS-Joaquim/KMZ/studyAreaControl.shp")
treated_area_sf <- st_read("D:/OneDrive - IHMT-NOVA/IHMT/Projectos/Moçambique-HS-Joaquim/KMZ/studyAreaTreat.shp")

# Ensure geometries are valid
control_area_sf_poly <- st_make_valid(control_area_sf)
treated_area_sf_poly <- st_cast(treated_area_sf, "POLYGON") %>% st_make_valid()

# Add labels
control_area_sf_poly$area_type <- "Control"
treated_area_sf_poly$area_type <- "Repellent"

# Combine
study_areas <- rbind(control_area_sf_poly[, c("geometry", "area_type")],
                     treated_area_sf_poly[, c("geometry", "area_type")])

# ----------------------------
# 2. Prepare Prediction Data
# ----------------------------

# Extract predicted values for only the observed data rows
n_obs <- nrow(agg_df)
lambda_mean <- result_combined$summary.fitted.values$mean[1:n_obs]

# Add predictions to agg_df and convert to sf
lambda_df <- agg_df %>%
  mutate(lambda_mean = lambda_mean) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


# ----------------------------
# 3. Plot Map
# ----------------------------

tmap_mode("plot")

tm_shape(study_areas) +
  tm_polygons(
    fill = "area_type",
    fill.alpha = 0.3,
    fill.scale = tm_scale(values = c("Control" = "skyblue", "Repellent" = "orange")),
    fill.legend = tm_legend(title = "Study Area")
  ) +
  tm_shape(lambda_df) +
  tm_symbols(
    fill = "lambda_mean",
    fill.scale = tm_scale_intervals(style = "quantile", values = "viridis"),
    size = 0.6,
    fill.legend = tm_legend(title = expression(Predicted~lambda))
  ) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_title("Predicted Mosquito Density in the Study Areas") +
  tm_layout(legend.outside = TRUE)

# ----------------------------
# 4. Save Map
# ----------------------------

tmap_save(
  filename = "D:/OneDrive - IHMT-NOVA/IHMT/Projectos/Moçambique-HS-Joaquim/predicted_mosquito_density_map.png",
  width = 10, height = 8, units = "in", dpi = 300
)

                     


library(leaflet)
library(leaflet.extras)
library(sf)
library(htmlwidgets)


# -----------------------------------------
# 1. Prepare aggregated data (agg_df)
# -----------------------------------------
agg_df$lon <- as.numeric(agg_df$lon)
agg_df$lat <- as.numeric(agg_df$lat)
agg_df <- agg_df[!is.na(agg_df$lon) & !is.na(agg_df$lat), ]

# Optionally filter to post-treatment only or a specific species
# agg_df <- subset(agg_df, treated_time == 1 & sp == "An_funestus")

# -----------------------------------------
# 2. Load shapefiles for treated & control areas
# -----------------------------------------
control_area_sf <- st_read("D:/OneDrive - IHMT-NOVA/IHMT/Projectos/Moçambique-HS-Joaquim/KMZ/studyAreaControl.shp")
treated_area_sf <- st_read("D:/OneDrive - IHMT-NOVA/IHMT/Projectos/Moçambique-HS-Joaquim/KMZ/studyAreaTreat.shp")

# -----------------------------------------
# 3. Build leaflet heatmap using count or normalized count
# -----------------------------------------
heatmap_map <- leaflet(agg_df) %>%
  addProviderTiles(providers$OpenStreetMap) %>%
  
  # 🔥 Heatmap based on normalized mosquito counts
  addHeatmap(lng = ~lon, lat = ~lat,
             intensity = ~norm_count,   # or use ~count
             blur = 20,
             max = 0.2,                # adjust based on data scale
             radius = 30,
             gradient = c(
               "0.1" = "blue",
               "0.3" = "#00FF00",
               "0.6" = "orange",
               "1.0" = "red"
             )) %>%
  
  # 🟪 Control area
  addPolygons(data = control_area_sf,
              fill = FALSE,
              color = "purple",
              weight = 3,
              opacity = 1,
              popup = "Control Area") %>%
  
  # 🟧 Treated area
  addPolygons(data = treated_area_sf,
              fill = FALSE,
              color = "orange",
              weight = 3,
              opacity = 1,
              popup = "Treated Area") %>%
  
  addLegend("bottomleft",
            colors = c("purple", "orange"),
            labels = c("Control Area", "Treated Area"),
            title = "Study Areas",
            opacity = 1)

# -----------------------------------------
# 4. Export map as HTML
# -----------------------------------------
saveWidget(heatmap_map,
           file = "D:/OneDrive - IHMT-NOVA/IHMT/Projectos/Moçambique-HS-Joaquim/maps/heatmap_agg_model_combined.html",
           selfcontained = TRUE)

# View in RStudio
heatmap_map

