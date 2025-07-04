# Charger les données
load(here::here("data", "derived_data", "P_mat.Rdata"))

library(INLA)
library(sf)
library(sp)
library(ggplot2)

## 1. Créer un maillage spatial
# Créer un SpatialPointsDataFrame
coordinates(P_mat) <- ~ longitude + latitude

# Créer le maillage triangulé (meshing de l'espace)
mesh <- inla.mesh.2d(
  loc = coordinates(P_mat),
  max.edge = c(0.5, 5),  # ajustable selon ton échelle spatiale
  cutoff = 0.2           # pour éviter les points trop proches
)
plot(mesh)

save(mesh, file = here::here("outputs", "mesh.Rdata"))

## 2. Construire le modèle spde
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(1, 0.5),   # P(rangecorr > 1) = 0.5
  prior.sigma = c(1, 0.01)   # P(sigma > 1) = 0.01
)

save(spde, file = here::here("outputs", "mesh.Rdata"))

## 3. Indexer les effets spatiaux
# Créer l'index spatial pour le GMRF
A <- inla.spde.make.A(mesh = mesh, loc = coordinates(P_mat))
spatial_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

## 4. Stack des données
stack <- inla.stack(
  data = list(y = P_mat$mean_p_maturity),
  A = list(1, A),
  effects = list(
    data.frame(
      Intercept = 1,
      protection_status = P_mat$protection_status,
      gdp = P_mat$log_gdp,
      TotalGravity = P_mat$log_TotalGravity,
      distance_to_port = P_mat$log_distance_to_port,
      depth = P_mat$log_depth,
      sst_av_7d = P_mat$sst_av_7d,
      chlorophyll_av = P_mat$log_chlorophyll_av,
      Source = P_mat$Source,
      survey_date = P_mat$survey_date
    ),
    spatial = spatial_index
  ),
  tag = "est"
)

## 5. Modélisation

formula <- y ~ 1 + protection_status + log_gdp + log_TotalGravity + log_distance_to_port +
  log_depth + sst_av_7d + log_chlorophyll_av + 
  f(Source, model = "iid") +
  f(survey_date, model = "iid") +
  f(spatial.field, model = spde)

model <- inla(
  formula,
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  family = "gaussian",  # ajuster selon ta variable réponse
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

## 6. Sauvegarde du modèle
save(model, file = here::here("outputs", "INLA_model.Rdata"))
