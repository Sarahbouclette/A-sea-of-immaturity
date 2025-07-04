
# Charger les données
load(here::here("data", "derived_data", "P_mat.Rdata"))

# Charger les packages
library(dplyr)
library(lubridate)
library(spaMM)
library(here)

# Transformation des données
P_mat <- P_mat %>%
  mutate(across(where(is.character), as.factor),
         survey_date = year(survey_date))

# Ajustement du modèle spaMM avec structure spatiale
model_spatial <- spaMM::fitme(
  formula = mean_p_maturity ~ protection + gdp + TotalGravity + distance_to_port +
    depth + sst_av_7d + chlorophyll_av + distance_to_coast +
    (1 | Source) + (1 | survey_date) + Matern(1 | longitude + latitude),
  data = P_mat,
  method = "ML",
  nb_cores = 4,
  verbose = TRUE
)

# Sauvegarde du modèle
save(model_spatial, file = here::here("outputs", "spaMM_model1.Rdata"))

