##-----------------Loading packages-------------------
library(tidyverse)
library(INLA)
library(lubridate)
library(ggplot2)

##-------------------------loading data-----------------------

load("~/Sea_of_immaturity/data/derived_data/P_mat_versions/P_mat_new.Rdata")

P_mat <- P_mat %>%
  mutate(across(where(is.character), as.factor), survey_date= lubridate::year(survey_date)) 

##---------------Load model----------------------------------

model_spatial <- readRDS(here::here("Sea_of_immaturity", "outputs", "05_explicative_modeling", "Model nb mat indiv", "INLA_model_nbinomial.rds"))

# Résumé et diagnostic

model_spatial$summary.fixed

signif_table <- model_spatial$summary.fixed %>%
  mutate(
    Variable = rownames(.),
    Significance = case_when(
      `0.025quant` > 0 ~ "Positive",
      `0.975quant` < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  ) %>%
  select(Variable, mean, `0.025quant`, `0.975quant`, Significance)

coeffs <- model_spatial$summary.fixed %>%
  mutate(
    lower = mean - 1.96 * sd,
    upper = mean + 1.96 * sd
  ) %>%
  mutate(var=rownames(coeffs))

# Supprimer intercept pour affichage
coeffs <- coeffs %>% filter(var != "(Intercept)")

# Affichage ggplot
ggplot(coeffs, aes(x = reorder(var, mean), y = mean)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  labs(
    x = "Variable explicative",
    y = "Effet marginal estimé (± IC 95%)",
    title = "Effets marginaux du modèle INLA",
    caption = "Barres = intervalles de crédibilité 95 %"
  ) +
  theme_minimal()

##---------------Accuracy du model---------------------------------------
# Y_obs : la variable réponse observée
# fitted : prédictions marginales du modèle INLA

fitted_vals <- model_spatial$summary.fitted.values$mean
Y_obs <- P_mat$Nb_mat_indiv  # ou autre réponse

# Var expliquée / Var totale
r2_bayes <- var(fitted_vals) / (var(fitted_vals) + mean(model_spatial$summary.fitted.values$sd^2))
print(r2_bayes)

##--------------Plot effet marginal----------------------------
library(INLA)

#liste de variables numériques continues
num_vars <- c("depth", "TotalGravity", "gdp", "chlorophyll_av", "sst_av_7d", "distance_to_port")

# Fonction pour transformer chaque marginal en dataframe avec quantiles

extract_marginal_linear <- function(var_name) {
  # Extraire les valeurs marginales de l'effet
  marginal <- model_spatial$marginals.fixed[[var_name]]
  
  # Estimations de l'effet
  mean_val <- inla.emarginal(function(x) x, marginal)
  lower <- inla.qmarginal(0.025, marginal)
  upper <- inla.qmarginal(0.975, marginal)
  
  # Générer une séquence de x autour des valeurs observées
  x_seq <- seq(min(P_mat[[var_name]], na.rm = TRUE), 
               max(P_mat[[var_name]], na.rm = TRUE), length.out = 100)
  
  # Estimer la relation linéaire
  df_pred <- data.frame(
    x = x_seq,
    fit = x_seq * mean_val,
    lower = x_seq * lower,
    upper = x_seq * upper,
    Variable = var_name
  )
  
  

  return(df_pred)
}

library(purrr)

# Appliquer la fonction à toutes les marginals.fixed
combined_eff_df <- map_dfr(num_vars, extract_marginal_linear)

labels_df <- combined_eff_df %>%
  group_by(Variable) %>%
  summarise(
    x = max(x, na.rm = TRUE),
    y = max(fit, na.rm = TRUE),
    SignifDirection = case_when(
      all(lower > 0) ~ "positive",
      all(upper < 0) ~ "negative",
      TRUE ~ "NS"
    )
  )

my_colors <- c(
  # Variables environnementales
  "sst_av_7d" = "#94A428",  "chlorophyll_av" = "#94A428", "depth" = "#94A428",
  
  # Variables anthropiques 
  "distance_to_port" = "#9EB3C6", "TotalGravity" = "#9EB3C6", "gdp" = "#9EB3C6",
  
  # Protection 
  "protection_status" = "#E03C3C"  # utilisé si une seule valeur
)

plot_num <- ggplot(combined_eff_df, aes(x = x, y = density, fill = Variable)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_text(data = labels_df, aes(x = x, y = y, label = SignifDirection),
            inherit.aes = FALSE, size = 3, fontface = "bold") +
  facet_wrap(~ Variable, scales = "free_x") +
  scale_fill_manual(values = my_colors, name = "Numerical Variables") +
  ylab("Effet marginal estimé") +
  xlab("") +
  theme_minimal()

##----------------Plot marginal predictions--------------------

extract_marginal <- function(var_name) {

 x_seq <- seq(min(P_mat[[var_name]], na.rm = TRUE), max(P_mat[[var_name]], na.rm = TRUE), length.out = 100)

# On fixe les autres variables à leur moyenne (ou médiane)
  new_data <- data.frame(
    var_name = x_seq,
    gdp = mean(P_mat$gdp, na.rm = TRUE),
    TotalGravity = mean(P_mat$TotalGravity, na.rm = TRUE),
    distance_to_port = mean(P_mat$distance_to_port, na.rm = TRUE),
    depth = mean(P_mat$depth, na.rm = TRUE),
    sst_av_7d = mean(P_mat$sst_av_7d, na.rm = TRUE),
    chlorophyll_av = mean(P_mat$chlorophyll_av, na.rm = TRUE),
    protection_status = P_mat$protection_status # ou la modalité de référence
  )
  
  Xmat <- model.matrix(~ gdp + TotalGravity + distance_to_port + depth + sst_av_7d + chlorophyll_av + protection_status, data = new_data)
  
  samples <- inla.posterior.sample(n = 1000, result = model_spatial)
  
  # Identifier les positions des effets fixes dans les samples
  fixed_idx <- model_spatial$misc$configs$contents$start[model_spatial$misc$configs$contents$tag == "Beta"]
  
  # Pour chaque sample, calculer la prédiction
  pred_matrix <- sapply(samples, function(s) {
    beta <- s$latent[fixed_idx:(fixed_idx + ncol(Xmat) - 1)]
    as.vector(Xmat %*% beta)
  })  # matrice [n_points x n_samples]
  
  pred_mean <- apply(pred_matrix, 1, mean)
  pred_lower <- apply(pred_matrix, 1, quantile, probs = 0.025)
  pred_upper <- apply(pred_matrix, 1, quantile, probs = 0.975)
  
  gdp_effect_df <- data.frame(
    gdp = gdp_seq,
    fit = pred_mean,
    lower = pred_lower,
    upper = pred_upper
  )
  
  ggplot(gdp_effect_df, aes(x = gdp, y = fit)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
    labs(x = "GDP", y = "Effet marginal estimé (échelle linéaire)", title = "Effet marginal de GDP") +
    theme_minimal()
  
}

##---------------Comparaison des distributions Obs vs Pred----------------

ggplot(data.frame(fitted = model_spatial$summary.fitted.values$mean), aes(x = fitted)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "#377EB8", alpha = 0.6) +
  geom_density(color = "black", size = 1) +
  theme_minimal() +
  labs(
    title = "Distribution des valeurs ajustées (fitted values)",
    x = "Valeurs ajustées",
    y = "Densité"
  )

ggplot(data.frame(
  fitted = model_spatial$summary.fitted.values$mean[1:32026],
  observed = P_mat$Nb_mat_indiv
)) +
  geom_density(aes(x = observed, color = "Observées"), size = 1) +
  geom_density(aes(x = fitted, color = "Ajustées"), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("Observées" = "#E41A1C", "Ajustées" = "#377EB8")) +
  theme_minimal() +
  labs(
    title = "Comparaison des distributions : Observées vs Ajustées",
    x = "Valeur de la variable réponse",
    y = "Densité",
    color = ""
  )

##----------------Distribution des résidus------------------

residus <- P_mat$mean_p_maturity - model_spatial$summary.fitted.values$mean[1:32026]

ggplot(data.frame(residus), aes(x = residus)) +
  geom_histogram(bins = 30, fill = "#6699CC", color = "white") +
  theme_minimal() +
  labs(x = "Résidus", y = "Fréquence", title = "Distribution des résidus")

qqnorm(residus)
qqline(residus, col = "red")

# Test de corrélation spatiale 
library(spdep)
library(sf)

##------------------Test corrélation spatiale-----------------
# 1. Charger les coordonnées spatiales
# (remplace lon et lat par les vrais noms de tes colonnes)
P_mat_sf <- st_as_sf(P_mat, coords = c("longitude", "latitude"), crs = 4326)

# 2. Créer une matrice de voisinage spatiale (basée sur k plus proches voisins)
coords <- st_coordinates(P_mat_sf)
nb <- knn2nb(knearneigh(coords, k = 5))  # k=5 voisins
listw <- nb2listw(nb, style = "W")       # pondérations spatiales

# 3. Extraire les résidus du modèle INLA
resid_spatial <- residuals(model_spatial)

# 4. Test de Moran sur les résidus
moran_test <- moran.test(resid_spatial$deviance.residuals, listw)

# 5. Affichage des résultats
print(moran_test)
