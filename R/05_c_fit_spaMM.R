rm(list=ls())

##-----------------Loading packages-------------------
library(tidyverse)
library(spaMM)
library(lubridate)
library(ggplot2)

##-------------------------loading data-----------------------

load("~/Sea_of_immaturity/data/derived_data/P_mat.Rdata")

P_mat <- P_mat %>%
  mutate(across(where(is.character), as.factor), survey_date= lubridate::year(survey_date)) %>%
 group_by(longitude, latitude) %>%
 slice_max(order_by = survey_date, n = 1, with_ties = FALSE) %>%
 ungroup()
   
 cols <- c("chlorophyll_av", "gdp", "TotalGravity", "distance_to_port", "depth")
 
 if (all(cols %in% colnames(P_mat))) {
   P_mat <- P_mat %>%
     mutate(across(all_of(cols), ~ log10(. + 1), .names = "log_{.col}")) %>%
     select(-all_of(cols))  # Supprime les anciennes colonnes
 }
  

##---------------Selection des variables-----------------------------------

##---------------Load model----------------------------------

model_spatial <- readRDS(here::here("Sea_of_immaturity", "outputs", "05_explicative_modeling", "Model Shannon", "spaMM_model_shannon.rds"))

# Résumé et diagnostic

summary(model_spatial)
ranefs <- ranef(model_spatial)
fixefs <- fixef(model_spatial)
AIC(model_spatial)

##---------------Tentative visualisation 2 -----------------------------

library(ggeffects)
library(effects)

# all effects
plot(allEffects(model_spatial))

# Exemple : effet marginal de la profondeur
#plot(ggpredict(model_spatial, terms = "depth[:]"))

vars_to_plot <- c("log_gdp", "log_TotalGravity", "sst_av_7d", "log_distance_to_port", "log_depth", "log_chlorophyll_av")

#Version élaborée
all_eff <- allEffects(model_spatial, xlevels=4)
# Convert all effects to a list of data frames
all_eff_dfs <- lapply(all_eff, as.data.frame)

# Optionally, add a column with the variable name
# Nettoyer chaque effet
all_eff_dfs_num <- Map(function(df, name) {
  df <- as.data.frame(df)
  
  # Renommer la colonne de la variable prédictive en "x"
  var_col <- setdiff(names(df), c("fit", "lower", "upper", "se", "residual.scale"))
  names(df)[names(df) == var_col] <- "x"
  
  df$Variable <- name
  df
}, all_eff[vars_to_plot], names(all_eff[vars_to_plot]))


# Combine all into a single data frame
combined_eff_df <- do.call(rbind, all_eff_dfs_num)


my_colors <- c(
  # Variables environnementales
  "sst_av_7d" = "#94A428",  
  "log_chlorophyll_av" = "#94A428", 
  "log_depth" = "#94A428",
  
  # Variables anthropiques 
  "log_distance_to_port" = "#9EB3C6",
  "log_TotalGravity" = "#9EB3C6", 
  "log_gdp" = "#9EB3C6",
  
  # Protection 
  "protection_status" = "#E03C3C"  # utilisé si une seule valeur
)

protection_colors <- c(
  "out" = "darkgrey",        # rouge très clair
  "low" = "#ff6666",    # rouge moyen
  "medium" = "#e60000",   # rouge vif
  "high" = "#990000"   # rouge foncé
)

all_eff_dfs$protection_status$protection_status <- factor(all_eff_dfs$protection_status$protection_status,
                               levels = c("out", "low", "medium", "high"),
                               ordered = TRUE)

combined_eff_df$Variable <- factor(combined_eff_df$Variable, levels = c("log_depth", "log_chlorophyll_av", "sst_av_7d", "log_distance_to_port", "log_TotalGravity", "log_gdp"), ordered = TRUE)
combined_eff_df <- combined_eff_df %>% left_join(confint_df  %>% rename(Variable = Covariate), by = "Variable")
combined_eff_df <- combined_eff_df %>%
  mutate(
    SignifDirection = case_when(
      Lower95 > 0 & Upper95 > 0 ~ "Positive",
      Upper95 < 0 & Lower95 < 0 ~ "Negative",
      TRUE ~ "NS"
    )
  )

# Choisir un point par facette pour afficher le label
labels_df <- combined_eff_df %>%
  group_by(Variable) %>%
  summarise(
    x = mean(x),                           # Moyenne des x pour positionner le texte horizontalement
    y = max(upper) + 0.1,                  # Légèrement au-dessus de l'intervalle
    SignifDirection = first(SignifDirection)
  )

library(ggsignif)
# Séparer les comparaisons
annotations_df <- contrast_results %>%
  tidyr::separate(comparison, into = c("xmin", "xmax"), sep = " vs ") %>%
  mutate(
    annotation = significance,
    y_position = pmax(
      all_eff_dfs$protection_status$fit[match(xmin, all_eff_dfs$protection_status$protection_status)],
      all_eff_dfs$protection_status$fit[match(xmax, all_eff_dfs$protection_status$protection_status)]
    ) + 0.05
  )

annotations_df <- annotations_df %>%
  mutate(
    base_y = pmax(
      all_eff_dfs$protection_status$fit[match(xmin, all_eff_dfs$protection_status$protection_status)],
      all_eff_dfs$protection_status$fit[match(xmax, all_eff_dfs$protection_status$protection_status)]
    ),
    y_position = base_y + 0.25 + row_number() * 0.05  # ⬅️ décalage progressif
  )
annotations_df <- as.data.frame(annotations_df)

# Example: plot them all in facets
plot_num <- ggplot(combined_eff_df, aes(x = x, y = fit, fill=Variable)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_text(data = labels_df, aes(x = x, y = y, label = SignifDirection),
            inherit.aes = FALSE, size = 3, fontface = "bold")  + 
  facet_wrap(~ Variable, scales = "free_x") +
  scale_fill_manual(values = my_colors, name="Numerical Variables") +
  ylab("CV Maturity Probability") +
  xlab("") +
  theme_minimal()

plot_cat <- ggplot(all_eff_dfs$protection_status, aes(x = protection_status, y = fit, color = protection_status)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 1.1) +
  theme_minimal() +
  scale_color_manual(values = protection_colors, name = "Protection Status") +
  ylab("CV Maturity Probability") +
  xlab("Protection status") +
  geom_signif(
    inherit.aes = FALSE,
    data = annotations_df,
    manual = TRUE,
    mapping = aes(xmin = xmin, xmax = xmax, annotations=annotation, y_position = y_position),
    tip_length = 0.01,
    textsize = 4,
    vjust = 0.5, 
  )

library(patchwork)
plot_num / plot_cat

##----------------Visualisation spatiale avec effets fixes-----------------------

# Adapter ces valeurs à ton jeu de données
lon_range <- range(P_mat$longitude)
lat_range <- range(P_mat$latitude)

# Grille régulière
grid <- expand.grid(
  longitude = seq(lon_range[1], lon_range[2], length.out = 100),
  latitude  = seq(lat_range[1], lat_range[2], length.out = 100)
)

# Fixer les autres variables (moyennes ou valeurs typiques)
grid <- grid %>%
  mutate(
    protection = "restricted.take",  # ou autre modalité par défaut
    gdp = mean(P_mat$gdp, na.rm = TRUE),
    TotalGravity = mean(P_mat$TotalGravity, na.rm = TRUE),
    distance_to_port = mean(P_mat$distance_to_port, na.rm = TRUE),
    depth = mean(P_mat$depth, na.rm = TRUE),
    sst_av_7d = mean(P_mat$sst_av_7d, na.rm = TRUE),
    chlorophyll_av = mean(P_mat$chlorophyll_av, na.rm = TRUE),
    Source = 'rls'  # s'il n'y a pas de valeur Source associée
  )

grid$predicted <- predict(model_spatial, newdata = grid, re.form = NULL)

ggplot(grid, aes(x = longitude, y = latitude, fill = predicted)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(fill = "Predicted maturity", title = "Effet spatial estimé par spaMM")


## -------------------Visualisation spatiale sans effets fixes-----------

grid$spatial_effect <- predict(model_spatial, newdata = grid, which = "rand")

library(ggplot2)

ggplot(grid, aes(x = longitude, y = latitude, fill = spatial_effect)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +
  coord_fixed() +
  labs(fill = "Effet spatial seul",
       title = "Champ spatial latent estimé par spaMM") +
  theme_minimal()

##---------------Partition de la variance------------------------

## Variance des effets randoms
ranpars <- spaMM::get_ranPars(model_spatial)
phi <- ranpars$phi  # variance résiduelle
lambda_spatial <- ranpars$lambda[3]
lambda_source <- ranpars$lambda[1]
lambda_time <- ranpars$lambda[2]

# Total de la variance modélisée
total_var <-  lambda_source + lambda_time + lambda_spatial + phi

# Partition de la variance
var_partition <- c(
  spatial = lambda_spatial / total_var,
  source = lambda_source / total_var,
  time = lambda_time / total_var,
  residual = phi / total_var
)

round(var_partition, 3)


# Variance expliquée par les effets fixes
# Approximation : var(predictions fixed only)
pred_fixed <- predict(model_spatial, newdata = P_mat, which = "fixed")
var_fixed <- var(pred_fixed)

# R² marginal (effets fixes)
R2_marginal <- var_fixed / (var_fixed + total_var)

# R² conditionnel (effets fixes + aléatoires)
pred_full <- predict(model_spatial, newdata = P_mat, which = "response")
var_full <- var(pred_full)
R2_conditional <- var_full / (var_full + phi)  # phi = bruit inexpliqué

c(R2_marginal = round(R2_marginal, 3), R2_conditional = round(R2_conditional, 3))


##---------------Reproduction figure-------------------------------

# Obtenir le tableau des coefficients estimés
beta_table <- spaMM::summary.HLfit(model_spatial)$beta_table %>%
  as.data.frame() %>%
  rownames_to_column(var = "Covariate")

library(parallel)
library(spaMM)

# Détection automatique des cœurs (optionnel : tu peux fixer le nombre)
n_cores <- 3  

# Parallélisation du calcul des intervalles de confiance
confint_list <- mclapply(beta_table$Covariate, function(var) {
  confint(model_spatial, parm = var, level = 0.95)
}, mc.cores = n_cores)

# Combine les résultats en data.frame (si besoin)
confint_fixef <- do.call(rbind, confint_list)

#save(confint_fixef, file="confint_fixef_model2.Rdata")

# Obtenir les intervalles de confiance
confint_fixef <- confint(model_spatial, parm = beta_table$Covariate, level = 0.95) 
confint_df <- do.call(rbind, lapply(confint_fixef, function(x) x$interval)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Covariate") %>%
  rename(Lower95 = `lower (Intercept)`, Upper95=`upper (Intercept)`)

# Fusionner les deux pour un tableau complet
coef_table <- beta_table %>%
  left_join(confint_df, by = "Covariate") %>%
  rename(Estimate = Estimate, CondSE = `Cond. SE`) %>%
  select(Covariate, Estimate, CondSE, Lower95, Upper95)  %>%
  mutate(Color = case_when(Lower95 < 0 & Upper95 < 0 ~ "red",
                           Lower95 > 0 & Upper95 > 0 ~ "blue",
                           TRUE ~ "gray")) %>% 
  mutate(Covariate = factor(x = Covariate, levels = rev(Covariate)))

# Créer un graphique des tailles d'effet
effSizeLims <- c(min(coef_table$Lower95), max(coef_table$Upper95))

effPlot <- ggplot() +
  geom_pointrange(coef_table, 
                  mapping = aes(x = Covariate,
                                y = Estimate,
                                colour = Color,
                                ymin = Lower95,
                                ymax = Upper95)) +
  scale_colour_manual(values = c("blue" = "blue",
                                 "gray" = "gray",
                                 "red" = "red")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "dark gray") +
  scale_y_continuous(limits = effSizeLims) +
  xlab("") +
  ylab("") +
  coord_flip()

effPlot

# Multicomp protection_status

# Créer un vecteur de contraste (même longueur que coef(model_spatial))
contrast <- rep(0, length(coefs))
names(contrast) <- names(coefs)

# Tous les niveaux du facteur
levels_prot <- c("high", "low", "medium", "out")

# Extraire les coefficients et la matrice de variance-covariance
coefs <- model_spatial$fixef
vcov_mat <- vcov(model_spatial)

# Fonction pour créer les contrastes
contrast_results <- combn(levels_prot, 2, simplify = FALSE) %>%
  map_df(~{
    l1 <- .x[1]; l2 <- .x[2]
    
    # Vecteur de contraste
    contrast <- rep(0, length(coefs))
    names(contrast) <- names(coefs)
    
    # Si l’un des niveaux est la référence ("high"), alors sa contribution = 0
    if (l1 != "high") contrast[paste0("protection_status", l1)] <- 1
    if (l2 != "high") contrast[paste0("protection_status", l2)] <- -1

    # Calculs
    est <- sum(contrast * coefs)
    se <- sqrt(t(contrast) %*% vcov_mat %*% contrast)
    z <- est / se
    p <- 2 * (1 - pnorm(abs(z)))
    
    tibble(
      comparison = paste(l1, "vs", l2),
      estimate = est,
      se = se,
      z = z,
      p_value = p,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  })

##----------Visualisation des résidus---------------------

library(ggplot2)

# Extraire les valeurs ajustées (prédictions du modèle sur les données d'origine)
fitted_vals <- fitted(model_spatial)

# Les tracer avec un histogramme + courbe de densité
ggplot(data.frame(fitted = fitted_vals), aes(x = fitted)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "#377EB8", alpha = 0.6) +
  geom_density(color = "black", size = 1) +
  theme_minimal() +
  labs(
    title = "Distribution des valeurs ajustées (fitted values)",
    x = "Valeurs ajustées",
    y = "Densité"
  )

epsilon <- 1e-6
P_mat$shannon_adj <- pmin(pmax(P_mat$shannon_index, epsilon), 1 - epsilon)
P_mat$shannon_logit <- log(P_mat$shannon_adj / (1 - P_mat$shannon_adj))

# Comparaison aux valeurs observées
ggplot(data.frame(
  fitted = fitted(model_spatial),
  observed = model_spatial$y
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

# Résidus bruts
residus <- residuals(model_spatial)
residus_std <- residuals(model_spatial, type = "pearson")

ggplot(data.frame(residus_std), aes(x = residus)) +
  geom_histogram(bins = 30, fill = "#6699CC", color = "white") +
  theme_minimal() +
  labs(x = "Résidus", y = "Fréquence", title = "Distribution des résidus")

qqnorm(residus)
qqline(residus, col = "red")

fitted_vals <- fitted(model_spatial)

ggplot(data.frame(fitted = fitted_vals, residus = residus), aes(x = fitted, y = residus)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(x = "Valeurs ajustées", y = "Résidus", title = "Résidus vs valeurs ajustées")

# Visualisation spatiale
coords <- data.frame(longitude = P_mat$longitude,
                     latitude = P_mat$latitude)

# Combiner pour visualiser
residus_df <- cbind(coords, residus = residus)

library(sf)

ggplot(residus_df, aes(x = longitude, y = latitude, color = residus)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(color = "Résidus", title = "Résidus spatiaux")



# Test de corrélation spatiale 
library(spdep)
library(sf)

# 1. Charger les coordonnées spatiales
# (remplace lon et lat par les vrais noms de tes colonnes)
P_mat_sf <- st_as_sf(P_mat, coords = c("longitude", "latitude"), crs = 4326)

# 2. Créer une matrice de voisinage spatiale (basée sur k plus proches voisins)
coords <- st_coordinates(P_mat_sf)
nb <- knn2nb(knearneigh(coords, k = 5))  # k=5 voisins
listw <- nb2listw(nb, style = "W")       # pondérations spatiales

# 3. Extraire les résidus du modèle spaMM
resid_spatial <- residuals(model_spatial)

# 4. Test de Moran sur les résidus
moran_test <- moran.test(resid_spatial, listw)

# 5. Affichage des résultats
print(moran_test)

##--------Visualisation des résidus avec DHARMa-------------

library(DHARMa)

# Génère des simulations
sim <- simulateResiduals(fittedModel = model_spatial, n = 250)

# Visualisation simple
plot(sim)

# Tests de normalité, dispersion, etc.
testDispersion(sim)
testUniformity(sim)
testZeroInflation(sim)


