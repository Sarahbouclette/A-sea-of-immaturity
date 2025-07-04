rm(list=ls())

##-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "ggplot2", "slam",
          "pbmcapply", "patchwork", "ggplot2",
          "tibble", "stringr", "hrbrthemes", "randomForest", "ranger", "caret", "ggpubr", "pheatmap", "tidyverse", "plotly")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(ggplot2)

##-------------loading data and functions-------------
# species traits 
load(file = here::here("Sea_of_immaturity","outputs", "02_fill_species_traits", "MF_species_traits_imputed.Rdata"))
load(file = here::here("Sea_of_immaturity","data", "derived_data", "01_gather_species_info", "species_traits_before_MF.Rdata"))

##-----------1.Preparing data and selection of explicative variables------------

#Add the phylogeny predictive variable doing an MCA

var_phylo <- c("species_name", "Class", "Order", "Family", "Genus") 

phylogeny <- species_traits_final |> 
  dplyr::select(all_of(var_phylo)) |> 
  tibble::column_to_rownames(var="species_name") |> 
  as.matrix()

dimensionality <- FactoMineR::MCA(phylogeny, ncp = 10, graph=T) 

phylo_space <- dimensionality[["ind"]][["coord"]]
colnames(phylo_space) <- gsub(" ", "", colnames(phylo_space))

traits_and_phylo <- cbind(phylo_space, species_traits_final)
Dims <- colnames(traits_and_phylo[,1:10])

# Add the already known length at maturity to the data inferred with Miss forest
traits_and_phylo <- traits_and_phylo %>%
  select(species_name, Lm_max, Lm_min, all_of(Dims)) %>%
  mutate(log_Lm_mean = log10(rowMeans(across(c(Lm_max, Lm_min)), na.rm = TRUE) + 1))

# LOG TRANSFORMATION highly right skewed predictors (already done)
#cols <- c("depth_min", "depth_max", "a", "b", "K", "Lm_mean")

#if (all(cols %in% colnames(species_traits_imputed))) {
#  species_traits_imputed <- species_traits_imputed %>%
#    mutate(across(all_of(cols), ~ log10(. + 1), .names = "log_{.col}")) %>%
#    select(-all_of(cols))  # Supprime les anciennes colonnes
#}
  
data_rf <- species_traits_imputed %>%
  left_join(traits_and_phylo, by = "species_name")

# Get rid of the "Rhincodus typus"
data_rf <- data_rf[data_rf$species_name!="Rhincodus typus",]

# Liste des prédicteurs potentiels basés sur des connaissances écologiques
selected_traits <- c( "TempPrefMean","Schooling", "TempPrefMin","Fertilization","ReproMode","TempPrefM","DemersPelag", "Vulnerability_fishing",
                      "IUCN_inferred_Loiseau23", "ClimVuln_SSP126","Troph", "Climate", "log_depth_min", "log_depth_max",  "log_a", "log_b", "log_LMax", "log_K", "Status")

# Filtrer les données sans valeurs manquantes pour l’entraînement
data_rf_known <- data_rf |>
  tidyr::drop_na(all_of(selected_traits), log_Lm_mean, all_of(Dims)) %>%  # Supprime les lignes avec NA uniquement dans ces colonnes
  select(all_of(selected_traits), log_Lm_mean, all_of(Dims), all_of(var_phylo), Lm_max, Lm_min) #Species with known LengthMat and known predictors


##-----------2. Evaluation des variables predictrices------------------------

# Sélection des variables optimales avec RFE (Recursive Feature Elimination)
cv_control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)

rfe_result <- rfe(data_rf_known %>% select(-log_Lm_mean, -all_of(var_phylo),-Lm_max, -Lm_min), data_rf_known$log_Lm_mean,
                  sizes = c(3, 5, 6), rfeControl = cv_control)

best_vars <- predictors(rfe_result)  # Variables retenues par RFE
print(best_vars)


##---------------3.Optimization Random Forest-------------
tune_grid <- expand.grid(mtry = c(2, 5, 10),
                         splitrule = c("variance",NULL),
                         min.node.size = c(1, 5, 10))

cv_control <- trainControl(method = "cv", number = 5)

bestTune<-data.frame(num_trees_values = c(100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050), MSE = NA, mtry = NA, splitrule=NA, min.node.size = NA)

for (n in bestTune$num_trees_values) {
  rf_tuned <- train(
    log_Lm_mean ~ ., 
    data = data_rf_known %>% select(all_of(best_vars), log_Lm_mean),
    method = "ranger", 
    trControl = cv_control, 
    tuneGrid = tune_grid,
    num.trees = n
  )
  
  bestTune[bestTune$num_trees_values==n,]$MSE <- rf_tuned$finalModel$prediction.error
  bestTune[bestTune$num_trees_values==n,]$mtry <- rf_tuned$bestTune$mtry
  bestTune[bestTune$num_trees_values==n,]$splitrule <- rf_tuned$bestTune$splitrule
  bestTune[bestTune$num_trees_values==n,]$min.node.size <- rf_tuned$bestTune$min.node.size
}


print(rf_tuned$bestTune)  # Meilleurs paramètres trouvés

# Erreur en fonction du nombre d'arbres
plot(bestTune$num_trees_values, bestTune$MSE,
     type='l',
     main = "Erreur en fonction du nombre d'arbres", 
     xlab = "Nombre d'arbres", 
     ylab = "Erreur")

##----------------4. Training Random Forest-------------------------
set.seed(123)  # Reproductibilité
num_folds <- 5  # Nombre de folds

folds <- createFolds(data_rf_known$log_Lm_mean, k = num_folds, list = TRUE)

all_predictions <- list()  # Stocker les prédictions par fold
cv_metrics <- data.frame(Fold = integer(), R2 = numeric(), RMSE = numeric(), Type = character())  # Stocke les R² et RMSE

models_list <- list()

for (i in seq_along(folds)) {
  train_idx <- unlist(folds[-i])  # Les autres folds servent pour le train
  test_idx <- folds[[i]]  # Le fold i sert pour le test
  
  train_fold <- data_rf_known[train_idx, ]
  test_fold <- data_rf_known[test_idx, ]
  
  rf_model <- ranger(log_Lm_mean ~ ., 
                     data = train_fold %>% select(all_of(best_vars), log_Lm_mean),
                     num.trees = 700, importance = "permutation",
                     mtry = 5, splitrule = "variance", min.node.size = 5)
  models_list[[i]] <- rf_model
  
  # Prédictions
  pred_train <- predict(rf_model, train_fold %>% select(all_of(best_vars)))$predictions
  pred_test <- predict(rf_model, test_fold %>% select(all_of(best_vars)))$predictions
  
  # Calcul du R² pour Train
  ss_res_train <- sum((train_fold$log_Lm_mean - pred_train)^2)
  ss_tot_train <- sum((train_fold$log_Lm_mean - mean(train_fold$log_Lm_mean))^2)
  r2_train <- 1 - (ss_res_train / ss_tot_train)
  
  # Calcul du RMSE pour Train
  rmse_train <- RMSE(pred_train, train_fold$log_Lm_mean)
  
  # Calcul du R² pour Test
  ss_res_test <- sum((test_fold$log_Lm_mean - pred_test)^2)
  ss_tot_test <- sum((test_fold$log_Lm_mean - mean(test_fold$log_Lm_mean))^2)
  r2_test <- 1 - (ss_res_test / ss_tot_test)
  
  # Calcul du RMSE pour Test
  rmse_test <- RMSE(pred_test, test_fold$log_Lm_mean)
  
  # Stockage des résultats
  cv_metrics <- rbind(cv_metrics, 
                      data.frame(Fold = i, R2 = r2_train, RMSE = rmse_train, Type = "Train"),
                      data.frame(Fold = i, R2 = r2_test, RMSE = rmse_test, Type = "Test"))
  
  # Calculer les résidus individuels et leur valeur absolue pour chaque prédiction
  residuals_train <- train_fold$log_Lm_mean - pred_train
  residuals_test <- test_fold$log_Lm_mean - pred_test
  
  # Stocker les prédictions en gardant les noms d'espèces et les résidus RMS
  all_predictions[[i]] <- list(
    train = data.frame(Species = train_fold$species_name, Observed = train_fold$log_Lm_mean, 
                       Predicted = pred_train, Fold = i, Set = "Train", 
                       Residual = residuals_train, Abs_Residual = abs(residuals_train), Class = train_fold$Class, Order = train_fold$Order),
    test = data.frame(Species = test_fold$species_name, Observed = test_fold$log_Lm_mean, 
                      Predicted = pred_test, Fold = i, Set = "Test", 
                      Residual = residuals_test, Abs_Residual = abs(residuals_test), Class = test_fold$Class, Order = test_fold$Order)
  )
}

# Fusionner les prédictions
df_predictions <- do.call(rbind, lapply(all_predictions, function(x) rbind(x$train, x$test)))


# Transformer les noms des folds en facteur pour l'affichage ggplot
cv_metrics$Fold <- as.factor(cv_metrics$Fold)
cv_metrics_long <- cv_metrics %>%
  pivot_longer(cols = c(R2, RMSE), names_to = "Metric", values_to = "value")


# Création du graphique
library(ggplot2)

ggplot(cv_metrics_long, aes(x = Fold, y = value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y") +  # Permet d'afficher R² et RMSE séparément
  labs(title = "Comparison R² and RMSE by fold (Train vs Test)", 
       x = "Fold", y = "Valeur") +
  scale_fill_manual(values = c("steelblue", "darkorange"))

ggplot(cv_metrics, aes(x = Fold, y = R2, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison R² by fold (Train vs Test)", 
       x = "Fold", y = "Valeur") +
  scale_fill_manual(values = c("steelblue", "darkorange"))


ggsave(file = here("Sea_of_immaturity", "figures","02_fill_species_traits", "RF_log_rsquared_by_fold_Lm_mean.png"), width = 8, height = 6)

df_test <- df_predictions %>% filter(Set == "Test")
df_test <- df_test %>% mutate(Predicted = Predicted, Observed = Observed, Residual = Residual, Abs_Residual= Abs_Residual)
df_train <- df_predictions %>% filter(Set == "Train")

r2_test_tot <- mean(cv_metrics[cv_metrics$Type=="Test",]$R2)
r2_train_tot <- mean(cv_metrics[cv_metrics$Type=="Train",]$R2)

# Visualisation pred train
ggplot(df_train, aes(x = Predicted, y = Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Comparaison Observed vs Predicted on the Train set", 
       x = "Observed values", y = "Predicted values") + 
  annotate("text", x = min(df_train$Observed), y = max(df_train$Predicted), 
           label = paste("R² =", round(r2_train_tot, 3)), 
           hjust = 0, size = 5, color = "blue") 
# Sauvegarde du graphique
ggsave(file = here("Sea_of_immaturity", "figures", "02_fill_species_traits", "RF_Lm_mean_log_regression_CV_train.png"), width = 8, height = 6)

# Visualisation pred test
ggplot(df_test, aes(x = Predicted, y = Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Comparaison Observed vs Predicted on the Test set", 
       x = "Observed values", y = "Predicted values") + 
  annotate("text", x = min(df_test$Observed), y = max(df_test$Predicted), 
           label = paste("R² =", round(r2_test_tot, 3)), 
           hjust = 0, size = 5, color = "blue") 

# Sauvegarde du graphique
ggsave(file = here("Sea_of_immaturity", "figures", "02_fill_species_traits", "RF_Lm_mean_log_regression_CV_test.png"), width = 8, height = 6)

## Interactif graphical

# Calculer R² dans l’échelle originale
ss_res_real <- sum((df_test$Observed - df_test$Predicted)^2)
ss_tot_real <- sum((df_test$Observed - mean(df_test$Observed))^2)
r2_real <- 1 - (ss_res_real / ss_tot_real)

print(r2_real)

# Créer une colonne avec les infos concaténées
df_test <- df_test %>%
  mutate(info_text = paste("Species:", Species, "<br>Residual:", round(Residual, 3), "<br>Order:", Order))

# Graphique ggplot pretty
my_colors <- c("Elasmobranchii" = "#EF6157", "Teleostei" = "#73A6AD")

p <- ggplot(df_test, aes(x = Observed, y = Predicted, text = info_text, color = Class)) + 
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = my_colors) +  # même palette personnalisée
  annotate("text", x = min(df_test$Observed), y = max(df_test$Predicted), 
           label = paste("R² =", round(r2_real, 3)), 
           hjust = 0, vjust = 1.5, size = 5, color = "red") +
  labs(
    title = "Observed vs Predicted values of log(Lm_mean)",
    x = "Observed log(Lm_mean)", 
    y = "Predicted log(Lm_mean)", 
    color = "Fish Class"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Rendre le graphique interactif avec plotly
ggplotly(p, tooltip = "text")

library(htmlwidgets)

# Sauvegarder le plot interactif
#htmlwidgets::saveWidget(ggplotly(p, tooltip = "text"), "interactif_plot_RF_Lm_mean.html")

##------------7. Prediction of the unknown maturity length----------------------
data_rf_unknown <- data_rf[is.na(data_rf$log_Lm_mean), ] %>% #Species with unknown maturity Length
  select(all_of(selected_traits), all_of(var_phylo), all_of(Dims))  %>%
  tidyr::drop_na(all_of(best_vars)) #but known predictors

# Créer une matrice pour stocker les prédictions
pred_matrix <- sapply(models_list, function(mod) {
  predict(mod, data_rf_unknown %>% select(all_of(best_vars)))$predictions
})

# Moyenne des prédictions
mean_pred <- rowMeans(pred_matrix)

# Calcul de l'écart-type pour chaque espèce
sd_pred <- apply(pred_matrix, 1, sd)

# IC 95% : Moyenne ± 1.96 * SD
ci_lower <- mean_pred - 1.96 * sd_pred
ci_upper <- mean_pred + 1.96 * sd_pred

##-------------8. Integration of the predictions into the initial dataset-------

# data_rf_unknown$log_Lm_mean_predicted <- mean_pred
# data_rf_unknown$log_Lm_max_lower <- ci_lower
# data_rf_unknown$log_Lm_min_upper <- ci_upper

# Pour revenir à l’échelle d’origine (si log10)
data_rf_unknown$Lm_mean<- 10^mean_pred-1
data_rf_unknown$Lm_min <- 10^ci_lower-1
data_rf_unknown$Lm_max <- 10^ci_upper-1


##-------------9. Plot results predictions--------------------------

# Extraire l'importance des variables
importance_df <- as.data.frame(rf_model$variable.importance) %>%
  rownames_to_column(var = "Variable") %>%
  rename( importance = `rf_model$variable.importance`)%>%
  arrange(desc(importance))  # Trier par importance décroissante

Predictors <- ggplot(importance_df, aes(x = reorder(Variable, importance), y = importance)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Relative Importance of Predictors",
    x = "Variable", 
    y = "Importance score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave(plot= Predictors, width = 8, height = 8, 
       file= here::here("Sea_of_immaturity", "figures", "02_fill_species_traits", "RF_Lm_max_W_whaleshark_log_importance_variables.png"))

figure_plot <- Predictors + p

save(figure_plot, file=here::here("Sea_of_immaturity", "figures", "02_fill_species_traits","Final Plot RF.png"))
##-----------10. Save------------------------------------------------------

data_rf_known <- data_rf_known %>% mutate(Lm_mean = 10^log_Lm_mean-1) %>% select(-log_Lm_mean)

LengthMat_predicted <- rbind(data_rf_known, data_rf_unknown)
LengthMat_predicted <- LengthMat_predicted %>%
  select(species_name, Class, Family, Order, Genus, Lm_min, Lm_max, Lm_mean, log_LMax)

LengthMat_predicted <- LengthMat_predicted %>% mutate(LMax = 10^log_LMax-1) %>% select(-log_LMax)

# LengthMax_sppredicted<-rbind(data_rf_known, data_rf_unknown)
# LengthMax_sppredicted <- LengthMax_sppredicted %>%
#   select(species_name, Class, Family, Order, Genus, log_LMax) %>%
#   mutate(log_LMax = 10^log_LMax - 1) %>%
#   rename(LMax = log_LMax)
# 
# # Identifier les espèces mal prédites
# species_replaced <- LengthMatMin_u_predicted$species_name[
#   LengthMatMin_u_predicted$Lm_max > LengthMax_sppredicted$LMax
# ]
# 
# Prediction_errors <- data_rf_unknown %>%
#   select(species_name, Class, Family, Order, Genus, log_Lm_mean, log_LMax) %>%
#   mutate(log_Lm_mean = 10^log_Lm_mean - 1, log_LMax = 10^log_LMax - 1) %>%
#   rename(Lm_max = log_Lm_mean, LMax = log_LMax)
# Prediction_errors <- Prediction_errors[Prediction_errors$species_name %in% species_replaced,]
# 
# save(Prediction_errors, file = here::here("Sea_of_immaturity", "outputs", "prediction_errors_RF4_Lm_max.Rdata" ))

save(LengthMat_predicted, file=here::here("Sea_of_immaturity", "outputs", "02_fill_species_traits", "RF_LengthMat_predictions.Rdata"))

##-----------11. Search for causality (optional)-----------------------------------
library(ape)
library(pheatmap)
library(dplyr)
library(vegan)

## Similarity based on maturity lengths
pred_matrix <- as.matrix(pred_unknown)
dist_pred <- as.matrix(dist(pred_matrix))
rownames(dist_pred) <- data_rf_unknown$species_name
colnames(dist_pred) <- data_rf_unknown$species_name
sim_pred <- 1 / (1 + dist_pred)

pheatmap(sim_pred, 
         clustering_method = "ward.D2", 
         main = "Clustering basé sur Random Forest",
         labels_row = rownames(sim_matrix),
         labels_col = colnames(sim_matrix), 
         fontsize_row = 6, fontsize_col = 6)

hc_rf <- hclust(as.dist(sim_matrix), method = "ward.D2")

# Découper en 5 clusters
clusters <- cutree(hc_rf, k = 5)

# Associer chaque espèce à son cluster
species_clusters <- data.frame(Species = rownames(sim_matrix), Cluster = clusters)
print(species_clusters)

## Similarity based on the phylogeny
taxo_df <- data_rf_unknown %>%
  select(all_of(var_phylo)) %>%
  na.omit()

taxo_df$species_name <- as.factor(taxo_df$species_name)

tree <- ape::as.phylo(~Class/Order/Family/Genus/species_name, data = taxo_df)

# Ajouter les noms des taxons supérieurs
tree$node.label <- unique(c(taxo_df$Class, taxo_df$Order, taxo_df$Family, taxo_df$Genus))

# Affichage de l'arbre phylogénétique
plot(tree, type = "cladogram", cex = 0.5, no.margin = TRUE, label.offset = 0.5)
nodelabels(tree$node.label, frame = "none", cex = 0.7, col = "red")

# Matrice de distances phylogénétiques
phylo_dist <- cophenetic(tree)  

# Matrice de similarité phylogénétique
phylo_sim <- 1 / (1 + phylo_dist)

# S'assurer que les matrices ont les mêmes espèces et le même ordre
common_species <- intersect(rownames(sim_pred), rownames(phylo_sim))

sim_pred_filtered <- sim_pred[common_species, common_species]
phylo_sim_filtered <- phylo_sim[common_species, common_species]

# Vérification des dimensions
print(dim(sim_pred_filtered))
print(dim(phylo_sim_filtered))

## Comparaison des matrices avec un test de Mantel
mantel_test_sim <- vegan::mantel(sim_matrix_filtered, phylo_sim_filtered, method = "spearman")

# Affichage des résultats
print(mantel_test_sim)
