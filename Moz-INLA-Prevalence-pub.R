
# ----------------------------#
# Required Packages (Install Once)
# ----------------------------#

# Uncomment and run if needed
# install.packages(c("tidyverse", "ggplot2", "readxl", "sf", "terra", "spdep", "INLA", 
#                    "spatstat.geom", "spatstat.core", "spatstat.explore", 
#                    "viridis", "leaflet", "ggmap", "htmlwidgets", "lubridate", "fmesher"))

# ----------------------------#
# Load Libraries
# ----------------------------#

library(tidyverse)
library(ggplot2)
library(readxl)
library(sf)
library(terra)
library(spdep)
library(INLA)
library(spatstat.geom)
library(spatstat.core)
library(spatstat.explore)
library(viridis)
library(leaflet)
library(ggmap)
library(htmlwidgets)
library(lubridate)
library(fmesher)

# ----------------------------#
# Load Prevalence Data
# ----------------------------#

file_path <- "prevalence.xlsx"
prevdf <- read_excel(file_path)

# ----------------------------#
# Data Cleaning & Transformation
# ----------------------------#

# Convert wide to long format
prevdf_long <- prevdf %>%
  mutate(across(matches("^P[1-3]"), as.character)) %>%
  pivot_longer(cols = matches("^P[1-3]"), names_to = c("timepoint", "variable"),
               names_pattern = "(P[1-3])(.*)", values_to = "value") %>%
  pivot_wider(names_from = variable, values_from = value)

# Prepare clean version
prevdflg <- prevdf_long %>%
  mutate(
    timepoint = factor(timepoint, levels = c("P1", "P2", "P3")),
    Mal = factor(na_if(str_trim(Mal), ""), levels = c("0", "1")),
    Idade = as.numeric(na_if(str_trim(Idade), "")),
    CriancasAF = as.numeric(CriancasAF),
    sex = factor(sex),
    Afnovo = factor(Afnovo),
    UltimaVezMalaria = as.character(UltimaVezMalaria),
    unique_id = paste0(ordem5abr25, "_", timepoint),
    year = recode(timepoint, P1 = 2021, P2 = 2022, P3 = 2023),
    month = recode(timepoint, P1 = 6, P2 = 6, P3 = 7)
  ) %>%
  rename(AFnovo = Afnovo, idade = Idade)

# Remove NAs in malaria outcome
prevdfclean <- filter(prevdflg, !is.na(Mal))
prevdfclean$Mal <- as.numeric(as.character(prevdfclean$Mal))

# Recode interaction and time indicators
prevdfclean <- prevdfclean %>%
  mutate(
    post = ifelse(timepoint == "P1", 0, 1),
    treated_time = as.numeric(tratPulv) * post,
    time_index = as.numeric(factor(timepoint, levels = c("P1", "P2", "P3"))),
    tratPulv = as.numeric(tratPulv),
    sex = as.numeric(as.character(sex)),
    idade = as.numeric(idade),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  filter(!is.na(idade), CriancasAF != 0, UltimaVezMalaria != 99)

# Group malaria history
prevdfclean <- prevdfclean %>%
  mutate(
    UltMal_grouped = factor(case_when(
      UltimaVezMalaria %in% c("1", "2") ~ "recent",
      UltimaVezMalaria == "3" ~ "past",
      UltimaVezMalaria == "4" ~ "unknown",
      UltimaVezMalaria == "5" ~ "never"
    ), levels = c("recent", "past", "unknown", "never"))
  )

# ----------------------------#
# Build SPDE Mesh and Index
# ----------------------------#

coords <- as.matrix(prevdfclean[, c("lon", "lat")])
mesh <- inla.mesh.2d(loc = coords, max.edge = c(0.02, 0.2), cutoff = 0.01)
spde <- inla.spde2.pcmatern(mesh, alpha = 2, prior.range = c(0.1, 0.01), prior.sigma = c(1, 0.01))
A <- inla.spde.make.A(mesh, loc = coords)
spatial_index <- inla.spde.make.index("spatial_field", n.spde = spde$n.spde)

# ----------------------------#
# INLA Stack Function Builder
# ----------------------------#

build_stack <- function(data, vars) {
  inla.stack(
    data = list(Mal = data$Mal),
    A = list(A, 1),
    effects = list(
      spatial_field = spatial_index,
      as.data.frame(vars)
    ),
    tag = "est"
  )
}

# ----------------------------#
# Model 5: Best Full Model
# ----------------------------#

stack5 <- build_stack(prevdfclean, list(
  intercept = 1,
  tratPulv = prevdfclean$tratPulv,
  post = prevdfclean$post,
  treated_time = prevdfclean$treated_time,
  sex = prevdfclean$sex,
  idade = prevdfclean$idade,
  CriancasAF = prevdfclean$CriancasAF,
  time_index = prevdfclean$time_index
))

formula5 <- Mal ~ 1 + tratPulv + post + treated_time + sex + idade + CriancasAF +
  f(time_index, model = "rw1") +
  f(spatial_field, model = spde)

result5 <- inla(
  formula5, family = "binomial",
  data = inla.stack.data(stack5),
  control.predictor = list(A = inla.stack.A(stack5), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

summary(result5)

# ----------------------------#
# GEE
# ----------------------------#

install.packages("geepack")
library(geepack)

gee_model <- geeglm(
  Mal ~ tratPulv + post + treated_time + idade + sex,
  data = prevdfclean,
  id = AFnovo,              # Cluster ID (e.g., village, household)
  family = binomial(link = "logit"),
  corstr = "exchangeable"        # Correlation structure (other options: "independence", "ar1")
)

summary(gee_model)

# ----------------------------#
# Visualize Predicted Risk
# ----------------------------#

index_est <- inla.stack.index(stack5, tag = "est")$data
prevdfclean$predicted_prob <- result5$summary.fitted.values$mean[index_est]

ggplot(prevdfclean, aes(x = lon, y = lat)) +
  geom_point(aes(color = predicted_prob), size = 3) +
  scale_color_viridis_c(name = "Predicted
malaria risk", option = "plasma") +
  facet_wrap(~timepoint) +
  theme_minimal()

# ----------------------------#
# Leaflet Interactive Map
# ----------------------------#

pal <- colorNumeric(palette = "plasma", domain = prevdfclean$predicted_prob)

leaflet(prevdfclean) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(
    lng = ~lon, lat = ~lat, radius = 5,
    color = ~pal(predicted_prob), stroke = FALSE,
    fillOpacity = 0.8,
    label = ~paste("Risk:", round(predicted_prob, 3),
                   "<br>Timepoint:", timepoint)
  ) %>%
  addLegend("bottomright", pal = pal, values = ~predicted_prob,
            title = "Predicted malaria risk", opacity = 1)





