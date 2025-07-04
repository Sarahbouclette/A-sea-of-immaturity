rm(list=ls())

##-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "ggplot2", "slam", "readr",
          "pbmcapply", "patchwork", "ggplot2",
          "tibble", "stringr", "hrbrthemes", "randomForest", "ranger", "caret", "ggpubr", "pheatmap", "tydr", "dbscan", "geosphere", "corrplot", "rstatix")
library(geosphere)
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(ggplot2)

##-------------------------loading data-----------------------

load("~/Sea_of_immaturity/data/derived_data/P_mat.Rdata")

#-----------log_transformation----------------------

numeric_cols <- colnames(P_mat)[sapply(P_mat, class) == "numeric"]
df <- tidyr::pivot_longer(P_mat, cols = all_of(numeric_cols) ,
                          names_to = "trait", values_to = "value") |> 
  dplyr::filter(!is.na(value))

distrib_traits <- ggplot(data = df)+
  aes(x=value, group=trait, fill=trait) +
  geom_histogram(aes(y = ..density.., fill = trait), bins = 20,
                 color = "grey40", alpha = 0.2) +
  facet_wrap(~trait, scales = "free") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )

cols <- c("chlorophyll_min", "chlorophyll_max", "chlorophyll_av", "gdp", "TotalGravity", "distance_to_port", "depth")

if (all(cols %in% colnames(P_mat))) {
  P_mat <- P_mat %>%
    mutate(across(all_of(cols), ~ log10(. + 1), .names = "log_{.col}")) %>%
    select(-all_of(cols))  # Supprime les anciennes colonnes
}

##--------------Sélection par Vif-------------------------
library(car)      # pour la fonction vif()

vif_selection_plot <- function(data, response, predictors, vif_threshold = 5) {
  remaining_vars <- predictors
  vif_history <- list()
  iteration <- 1
  
  formula_str <- function(vars) paste(response, "~", paste(vars, collapse = " + "))
  
  repeat {
    form <- as.formula(formula_str(remaining_vars))
    lm_model <- lm(form, data = data)
    
    vifs <- car::vif(lm_model)
    vif_df <- data.frame(variable = names(vifs), VIF = as.numeric(vifs), iteration = iteration)
    vif_history[[iteration]] <- vif_df
    
    max_vif <- max(vifs)
    
    if (max_vif < vif_threshold) break
    
    var_to_remove <- names(which.max(vifs))
    message("Retrait de la variable : ", var_to_remove, " (VIF = ", round(max_vif, 2), ")")
    remaining_vars <- setdiff(remaining_vars, var_to_remove)
    iteration <- iteration + 1
  }
  
  # Combine les VIF de chaque itération
  vif_all <- do.call(rbind, vif_history)
  
  # Plot
  vif_plot <- ggplot(vif_all, aes(x = variable, y = VIF, fill = factor(iteration))) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = vif_threshold, color = "red", linetype = "dashed") +
    labs(title = "Évolution des VIF par variable et itération",
         x = "Variable", y = "VIF", fill = "Itération") +
    theme_minimal() +
    coord_flip()
  
  print(vif_plot)
  
  return(remaining_vars)
}

selected_vars <- vif_selection_plot(P_mat, "mean_p_maturity", c("log_distance_to_port", "log_gdp", "log_TotalGravity", "sst_min_7d", "sst_av_7d", "sst_max_7d","sst_min_3y", "sst_av_3y", "sst_max_3y", "log_depth", "log_chlorophyll_min", "log_chlorophyll_max", "log_chlorophyll_av"), vif_threshold = 5)
selected_vars_cv <- vif_selection_plot(P_mat, "cv_p_maturity", c("log_distance_to_port", "log_gdp", "log_TotalGravity", "sst_min_7d", "sst_av_7d", "sst_max_7d", "log_depth", "log_chlorophyll_min", "sst_min_3y", "sst_av_3y", "sst_max_3y", "log_chlorophyll_max", "log_chlorophyll_av"), vif_threshold = 5)
selected_vars <- vif_selection_plot(P_mat, "mean_p_maturity", c("distance_to_port", "gdp", "TotalGravity", "sst_min", "sst_av", "sst_max", "depth", "chlorophyll_min", "chlorophyll_max", "chlorophyll_av"), vif_threshold = 5)

##--------------Selection des variables numériques----------

# 1. Sélection des variables numériques à corréler
vars_corr <- P_mat %>%
  select(where(is.numeric), -latitude, -longitude) %>%  # garder seulement les variables numériques
  select(-mean_p_maturity) %>%   # enlever la variable réponse
  drop_na()                      # retirer les lignes avec NA (sinon cor() plante)

# 2. Calcul de la matrice de corrélation de Pearson
cor_matrix <- cor(vars_corr, method = "pearson")

# 3. Visualisation de la matrice
corrplot::corrplot(cor_matrix, method = "color", type = "upper", 
                   tl.col = "black", tl.srt = 45, 
                   addCoef.col = "black", number.cex = 0.7)

##----------------Sélection des variables catégorielles-------------

# 1. Charger les packages
library(rstatix)    # pour cramers_v

# 2. Sélection des colonnes catégorielles
vars_cat <- P_mat %>%
  dplyr::select(protection, protection_status, protection_status2) %>%
  drop_na()

# 3. Convertir en facteurs si ce n'est pas déjà fait
vars_cat <- vars_cat %>% mutate(across(everything(), as.factor))

# 4. Fonction pour calculer Cramér's V entre chaque paire de variables
get_cramers_v_matrix <- function(df) {
  var_names <- colnames(df)
  n <- length(var_names)
  result_matrix <- matrix(NA, n, n, dimnames = list(var_names, var_names))
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        result_matrix[i, j] <- 1
      } else {
        tbl <- table(df[[i]], df[[j]])
        v <- suppressWarnings(rstatix::cramer_v(tbl))
        result_matrix[i, j] <- v
      }
    }
  }
  
  return(result_matrix)
}

# 5. Calcul de la matrice
cramer_matrix <- get_cramers_v_matrix(vars_cat)

# 6. Visualisation
corrplot(cramer_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7)

#È. Visualisation avec la corrélation des données numériques

par(mfrow = c(1, 2))  # 1 ligne, 2 colonnes

corrplot(cor_matrix, method = "color", type = "upper",
         title = "Corrélation de Pearson", mar = c(0,0,2,0),
         tl.col = "black", addCoef.col = "black", number.cex = 0.7)

corrplot(cramer_matrix, method = "color", type = "upper",
         title = "Cramér's V (Catégorielles)", mar = c(0,0,2,0),
         tl.col = "black", addCoef.col = "black", number.cex = 0.7)

par(mfrow = c(1, 1))  # reset


##-------------------Fit RF------------------------

## Fit modèle avec composantes spatiales 

categorize_latitude <- function(lat) {
  lat_abs <- abs(lat)
  
  if (lat_abs < 23.5) {
    return("Tropical")
  } else if (lat_abs < 35) {
    return("Subtropical")
  } else if (lat_abs < 60) {
    return("Temperate")
  } else {
    return("Polar/Subpolar")
  }
}

# Application
zone_climatique <- sapply(P_mat$latitude, categorize_latitude)

# Résultat
P_mat <- P_mat %>%
  mutate(Latitudinale_zone = zone_climatique, survey_date = lubridate::year(survey_date))

library(ranger)
library(pdp)
library(ggplot2)

set.seed(123)  # Reproductibilité
num_folds <- 5  # Nombre de folds

## Fit le modèle 

# folds <- createFolds(P_mat$mean_p_maturity, k = num_folds, list = TRUE)
# 
# all_predictions <- list()  # Stocker les prédictions par fold
# cv_metrics <- data.frame(Fold = integer(), R2 = numeric(), RMSE = numeric(), Type = character())  # Stocke les R² et RMSE
# 
# for (i in seq_along(folds)) {
#   train_idx <- unlist(folds[-i])  # Les autres folds servent pour le train
#   test_idx <- folds[[i]]  # Le fold i sert pour le test
#   
#   train_fold <- P_mat[train_idx, ]
#   test_fold <- P_mat[test_idx, ]
#   
#   rf_model2 <- ranger::ranger(mean_p_maturity~ ., 
#                              data = train_fold %>% select(all_of(selected_vars), mean_p_maturity, protection_status, Source, Latitudinale_zone, latitude, longitude, survey_date),
#                              num.trees = 600, importance = "permutation",
#                              mtry = 3, splitrule = "variance", min.node.size = 1)
#   
#   # Prédictions
#   pred_train <- predict(rf_model, train_fold %>% select(select(all_of(selected_vars), mean_p_maturity, protection_status, Source, Latitudinale_zone, latitude, longitude, survey_date)))$predictions
#   pred_test <- predict(rf_model, test_fold %>% select(select(all_of(selected_vars), mean_p_maturity, protection_status, Source, Latitudinale_zone, latitude, longitude, survey_date)))$predictions
#   
#   # Calcul du R² pour Train
#   ss_res_train <- sum((train_fold$mean_p_maturity - pred_train)^2)
#   ss_tot_train <- sum((train_fold$mean_p_maturity - mean(train_fold$mean_p_maturity)^2))
#   r2_train <- 1 - (ss_res_train / ss_tot_train)
#   
#   # Calcul du RMSE pour Train
#   rmse_train <- RMSE(pred_train, train_fold$log_Lm_max)
#   
#   # Calcul du R² pour Test
#   ss_res_test <- sum((test_fold$mean_p_maturity - pred_test)^2)
#   ss_tot_test <- sum((test_fold$mean_p_maturity - mean(test_fold$mean_p_maturity))^2)
#   r2_test <- 1 - (ss_res_test / ss_tot_test)
#   
#   # Calcul du RMSE pour Test
#   rmse_test <- RMSE(pred_test, test_fold$mean_p_maturity)
#   
#   # Stockage des résultats
#   cv_metrics <- rbind(cv_metrics, 
#                       data.frame(Fold = i, R2 = r2_train, RMSE = rmse_train, Type = "Train"),
#                       data.frame(Fold = i, R2 = r2_test, RMSE = rmse_test, Type = "Test"))
#   
#   # Calculer les résidus individuels et leur valeur absolue pour chaque prédiction
#   residuals_train <- train_fold$mean_p_maturity - pred_train
#   residuals_test <- test_fold$mean_p_maturity - pred_test
#   
#   # Stocker les prédictions en gardant les noms d'espèces et les résidus RMS
#   all_predictions[[i]] <- list(
#     train = data.frame(Observed = train_fold$mean_p_maturity, 
#                        Predicted = pred_train, Fold = i, Set = "Train", 
#                        Residual = residuals_train, Abs_Residual = abs(residuals_train), Class = train_fold$Class, Order = train_fold$Order),
#     test = data.frame(Species = test_fold$species_name, Observed = test_fold$mean_p_maturity, 
#                       Predicted = pred_test, Fold = i, Set = "Test", 
#                       Residual = residuals_test, Abs_Residual = abs(residuals_test), Class = test_fold$Class, Order = test_fold$Order)
#   )
# }
# 
# # Fusionner les prédictions
# df_predictions <- do.call(rbind, lapply(all_predictions, function(x) rbind(x$train, x$test)))

rf_model <- ranger::ranger(mean_p_maturity~ ., 
                           data = P_mat %>% select(all_of(selected_vars), mean_p_maturity, protection_status, Source, latitude, longitude, survey_date, Latitudinale_zone),
                           num.trees = 700, importance = "permutation",
                           mtry = 3, splitrule = "variance", min.node.size = 1)


rf_model2 <- ranger::ranger(mean_p_maturity~ ., 
                            data = P_mat %>% select(all_of(selected_vars), mean_p_maturity, protection_status, Source, Realm, latitude, longitude, survey_date),
                            num.trees = 700, importance = "permutation",
                            mtry = 3, splitrule = "variance", min.node.size = 1)

rf_model3 <- ranger::ranger(mean_p_maturity~ ., 
                            data = P_mat %>% select(all_of(selected_vars), mean_p_maturity, protection_status, Source, Province, survey_date),
                            num.trees = 700, importance = "permutation",
                            mtry = 3, splitrule = "variance", min.node.size = 1)

rf_model_cv <- ranger::ranger(cv_p_maturity~ ., 
                              data = P_mat %>% select(all_of(selected_vars_cv), cv_p_maturity, protection_status, Source, latitude, longitude, survey_date) %>% na.omit(cv_p_maturity),
                              num.trees = 700, importance = "permutation",
                              mtry = 3, splitrule = "variance", min.node.size = 1)


##-----------Importance des variables-------------

my_colors <- c("latitude"="#9EB3C2", "longitude"="#9EB3C2", "Source"="#9EB3C2","Latitudinal_zone"="#9EB3C2", "survey_date"="#9EB3C2",
               "sst_av_7d"="#94A428", "sst_min_3y"="#94A428", "log_chlorophyll_min"="#94A428", "log_chlorophyll_max"="#94A428", "log_depth" = "#94A428",
               "log_distance_to_port"= "#923A54", "log_TotalGravity"= "#923A54", "log_gdp"= "#923A54",
               "protection_status"= "#87CDE8", "Province"="#9EB3C2")

importance_df <- as.data.frame(rf_model_cv$variable.importance) %>%
  rownames_to_column(var = "Variable") %>%
  rename( importance = `rf_model_cv$variable.importance`)%>%
  arrange(desc(importance))  # Trier par importance décroissante

Variables_importance <- ggplot(importance_df, aes(x = reorder(Variable, importance), y = importance, fill=Variable)) +
  geom_col() +
  scale_fill_manual(values = my_colors) +
  coord_flip() +
  labs(title = "Importance of the variables", x = "Variable", y = "Importance") 


##------------Plot la relation entre mean_p_mat et chaque variable---------

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(purrr)

# Compute la relation
partial_dep_function <- function(model, covariates, data, reponse_name) {
  partial_i <- pbmcapply::pbmclapply(seq_along(covariates), function(i) {
    var_name <- covariates[i]
    
    # Calcul de la dépendance partielle
    partial_dep <- pdp::partial(
      object = model,
      train = data,
      pred.var = var_name,
      chull = TRUE, 
      ice=TRUE,
    )
    
    return(partial_dep)
  }, mc.cores = 3)
  
  names(partial_i) <- covariates
  return(partial_i)
}


plot_partial_dep <- function(partial_dep_list, color_used) {
  
  # Créer une colonne "type" : numérique ou catégorielle
  df_list <- imap(partial_dep_list, function(df, varname) {
    df$varname <- varname
    df$type <- if (is.factor(df[[1]]) || is.character(df[[1]])) "categorical" else "numeric"
    df
  })
  
  num_df <- data.frame()
  cat_df <- data.frame()
  
  for (i in seq_along(df_list)) {
    colnames(df_list[[i]])[1] <- "val_var"
    
    if (df_list[[i]]$type[1] == "categorical") {
      cat_df <- bind_rows(cat_df, df_list[[i]])
    } else {
      num_df <- bind_rows(num_df, df_list[[i]])
    }
  }
  
  # --- Graphique pour variables numériques ---
  if (nrow(num_df) > 0) {
    summary_ice <- num_df %>%
      group_by(varname, across(1)) %>%
      summarise(
        mean_pred = mean(yhat, na.rm = TRUE),
        sd_pred = sd(yhat, na.rm = TRUE),
        n = n(),
        se_pred = sd_pred / sqrt(n),
        ci_low = mean_pred - 1.96 * se_pred,
        ci_high = mean_pred + 1.96 * se_pred,
        .groups = "drop"
      )
    
    plot_num <- ggplot(summary_ice, aes(x = val_var, y = mean_pred)) +
      geom_line(color = color_used, size = 1, stat='smooth') +
      geom_point(color = "darkgrey", size = 1.5) +
      facet_wrap(~varname, scales = "free_x") +
      theme_minimal() +
      labs(y = "Prédiction moyenne ± IC 95%", x = "Valeur de la variable")
  } else {
    plot_num <- NULL
  }
  
  # --- Graphique pour variables catégorielles ---
  if (nrow(cat_df) > 0) {
    plot_cat <- ggplot(cat_df, aes(x = .data[[names(cat_df)[1]]], y = yhat, fill = .data[[names(cat_df)[1]]])) +
      #stat_compare_means(method = "t.test") +
      geom_boxplot() +
      facet_wrap(~varname, scales = "free") +
      theme_minimal() +
      labs(y = "Valeur prédite (ICE)", x = "Catégorie") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    plot_cat <- NULL
  }
  
  return(list(
    plot_num,
    plot_cat
  ))
}

##-------------Plot results rf_model------------------------
partial_effects_hactivities <- partial_dep_function(rf_model, c("log_gdp", "log_distance_to_port", "log_TotalGravity"), P_mat, "mean_p_maturity")
partial_effect_env <- partial_dep_function(rf_model, c("log_chlorophyll_min", "log_chlorophyll_max", "sst_av_7d", "log_depth", "sst_min_3y"), P_mat, "mean_p_maturity")
partial_rnd_effect <- partial_dep_function(rf_model, c("latitude", "longitude"), P_mat, "mean_p_maturity")
partial_effect_categ <- partial_dep_function(rf_model, c("Source", "protection_status", "Latitudinale_zone"), P_mat, "mean_p_maturity")

plots <- plot_partial_dep(partial_effects)
# Afficher les deux figures combinées
plot_num <- plots[[1]]
plot_cat <- plots[[2]]

plot_effects_hactivities <- plot_partial_dep(partial_effects_hactivities, "#923A54")[[1]]
plot_effect_env <- plot_partial_dep(partial_effect_env, "#94A428")[[1]]
plot_effect_rnd <- plot_partial_dep(partial_rnd_effect, "#9EB3C2" )[[1]]
plot_effect_categ <- plot_partial_dep(partial_effect_categ, "white")[[2]]

(plot_effects_hactivities + plot_effect_env + plot_effect_rnd) | plot_effect_categ
##-----------Plot results cv_rf_model---------------------------------

partial_effects_hactivities <- partial_dep_function(rf_model_cv, c("log_gdp", "log_distance_to_port", "log_TotalGravity"), P_mat, "cv_p_maturity")
partial_effect_env <- partial_dep_function(rf_model_cv, c("log_chlorophyll_min", "log_chlorophyll_max", "sst_av_7d", "log_depth", "sst_min_3y"), P_mat, "cv_p_maturity")
partial_rnd_effect <- partial_dep_function(rf_model_cv, c("latitude", "longitude"), P_mat, "cv_p_maturity")
partial_effect_categ <- partial_dep_function(rf_model_cv, c("Source", "protection_status", "Latitudinale_zone"), P_mat, "cv_p_maturity")

plots <- plot_partial_dep(partial_effects)
# Afficher les deux figures combinées
plot_num <- plots[[1]]
plot_cat <- plots[[2]]

plot_effects_hactivities_cv <- plot_partial_dep(partial_effects_hactivities, "#923A54")[[1]]
plot_effect_env_cv <- plot_partial_dep(partial_effect_env, "#94A428")[[1]]
plot_effect_rnd_cv <- plot_partial_dep(partial_rnd_effect, "#9EB3C2" )[[1]]
plot_effect_categ_cv <- plot_partial_dep(partial_effect_categ, "white")[[2]]

(plot_effects_hactivities_cv + plot_effect_env_cv + plot_effect_rnd_cv) | plot_effect_categ_cv