library(ggplot2)
library(dplyr)
library(funbiogeo)

# Human pressure variables
load("~/Sea_of_immaturity/outputs/04_get_vars/human pressures/distance_to_port.Rdata")
load("~/Sea_of_immaturity/outputs/04_get_vars/human pressures/sites_with_gdp.RData")
load("~/Sea_of_immaturity/outputs/04_get_vars/human pressures/sites_with_gravity.Rdata")
var_d_to_coast <- read.csv("Sea_of_immaturity/outputs/04_get_vars/env/var_env_7days_4_stats.csv") %>%
  rename(`Sample ID`=id)

# Protection variables
load("~/Sea_of_immaturity/outputs/04_get_vars/protection/mpa_mondial_protection_status_v2.Rdata")

# Env variables
var_7days_chlorophyll <- read.csv("Sea_of_immaturity/outputs/04_get_vars/env/var_env_7days_2_stats.csv") %>%
  rename(`Sample ID`=id)
var_7days_sst <- read.csv("Sea_of_immaturity/outputs/04_get_vars/env/var_env_7days_1_stats.csv") %>%
  rename(`Sample ID`=id, sst_av_7d = sst_av, sst_min_7d = sst_min, sst_max_7d = sst_max)
var_3years_sst <- read.csv("Sea_of_immaturity/outputs/04_get_vars/env/var_env_3_years_1_stats.csv") %>%
  rename(`Sample ID`=id, sst_av_3y = sst_av, sst_min_3y = sst_min, sst_max_3y = sst_max)
load("~/Sea_of_immaturity/outputs/04_get_vars/env/mondial_ecoregion.Rdata")

# load Measurment file with the individuals probability of being mature
load("~/Sea_of_immaturity/data/derived_data/Measurment_data_clustered.Rdata")

##--------------------Join variables ----------------------------------------

site_data <- Measurment_data_clustered %>%
  group_by(latitude, longitude, survey_date, Source, `Sample ID`) %>%
  summarise(
  mean_p_maturity = mean(p_maturity, na.rm = TRUE)) %>%
  tidyr::drop_na(latitude, longitude, survey_date)%>%
  ungroup()


# Add explicative variables
site_data <- site_data %>%
  left_join(sites_with_d_to_port %>% select(`Sample ID`, distance_to_port), by = "Sample ID") %>%
  left_join(sites_with_gravity %>% select(`Sample ID`, TotalGravity), by = "Sample ID") %>%
  left_join(gdp %>% select(`Sample ID`, gdp), by = "Sample ID") %>%
  left_join(sites_w_mpa_status %>% rename(`Sample ID` = survey_id) %>% select(`Sample ID`, protection, protection_status, inside_mpa), by="Sample ID")%>%
  left_join(var_7days_chlorophyll %>% select(`Sample ID`, chlorophyll_av, chlorophyll_min, chlorophyll_max), by ="Sample ID") %>%
  left_join(var_7days_sst %>% select(`Sample ID`, sst_av_7d, sst_min_7d, sst_max_7d), by ="Sample ID") %>%
  left_join(var_3years_sst %>% select(`Sample ID`, sst_av_3y, sst_min_3y, sst_max_3y), by ="Sample ID") %>%
  left_join(mondial_ecoregion %>% select(`Sample ID`, Province, Realm), by ="Sample ID") %>%
  left_join(var_d_to_coast %>% select(`Sample ID`, distance.to.coast_av), by ="Sample ID") %>%
  mutate(distance.to.coast = as.numeric(distance.to.coast_av)) %>% select(-distance.to.coast_av)


Measurment_data_completed <- Measurment_data_clustered %>%
  left_join(site_data %>% dplyr::select(`Sample ID`, protection_status, protection, inside_mpa, Province, Realm, distance.to.coast,
                                        gdp, TotalGravity, distance_to_port, 
                                        sst_av_7d, sst_min_7d, sst_max_7d, chlorophyll_av, chlorophyll_min, chlorophyll_max, sst_av_3y, sst_min_3y, sst_max_3y ), by = "Sample ID") 

##---------------------P_maturity by site-------------------------------------

# Moyenne par cluster
mean_p_by_cluster <- Measurment_data_completed %>%
  group_by(cluster) %>%
  filter(n() >= 5) %>%  # Garde uniquement les clusters avec au moins 5 individus
  summarise(
    mean_p_maturity = mean(p_maturity, na.rm = TRUE),
    .groups = "drop"
  )

# Résumer les infos additionnelles par cluster
cluster_info <- Measurment_data_completed %>%
  group_by(cluster) %>%
  filter(n() >= 5) %>%
  summarise(
    depth = mean(depth),
    sst_av_7d = mean(sst_av_7d),
    sst_min_7d = mean(sst_min_7d),
    sst_max_7d = mean(sst_max_7d),
    sst_av_3y = mean(sst_av_3y),
    sst_min_3y = mean(sst_min_3y),
    sst_max_3y = mean(sst_max_3y),
    chlorophyll_av = mean(chlorophyll_av),
    chlorophyll_max = mean(chlorophyll_max),
    chlorophyll_min = mean(chlorophyll_min),
    survey_date = first(survey_date),
    protection_status = first(protection_status),
    protection = first(protection),
    gdp = mean(gdp, na.rm = TRUE),
    TotalGravity = mean(TotalGravity, na.rm = TRUE),
    distance_to_port = mean(distance_to_port, na.rm = TRUE),
    latitude = mean(latitude, na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE),
    Source = first(Source),
    Province = first(Province),
    Realm = first(Realm),
    distance.to.coast = mean(distance.to.coast, na.rm=TRUE),
    .groups = "drop"
  )

# Joindre avec les moyennes de maturité par cluster
P_mat <- mean_p_by_cluster %>%
  left_join(cluster_info, by = "cluster") %>%
  mutate(
    localisation = paste0(round(longitude, 5), ":", round(latitude, 5)),
    mean_p_maturity = ifelse(mean_p_maturity == 0, 1e-3, mean_p_maturity)
  ) %>%
  tidyr::drop_na(latitude, longitude, mean_p_maturity, survey_date)

P_mat <- P_mat %>%
  mutate(across(where(is.character), as.factor))

cv_p_by_cluster <- Measurment_data_completed %>%
  group_by(cluster) %>%
  filter(n() >= 5) %>%
  summarise(
    mean_p_maturity = mean(p_maturity, na.rm = TRUE),
    cv_p_maturity = ifelse(is.na(sd(p_maturity)), 1e-6/mean_p_maturity,sd(p_maturity, na.rm = TRUE)/mean_p_maturity),
    var_p_maturity = ifelse(is.na(sd(p_maturity)), 1e-6, sd(p_maturity, na.rm = TRUE)),
    .groups = "drop"
  )

P_mat <- P_mat %>%
  left_join(cv_p_by_cluster %>% select(cluster, cv_p_maturity, var_p_maturity), by ="cluster")

#-------------------------Variables completeness-------------------
P_mat <- P_mat %>%
  rename( cluster = species)

fb_plot_number_species_by_trait(P_mat)

## Correct unpossible values of sst
P_mat <- P_mat %>%
  mutate(
    sst_av_7d = if_else(!is.na(sst_av_7d) & (sst_av_7d < 250 | sst_av_7d >330), NA, sst_av_7d),
    sst_min_7d = if_else(!is.na(sst_min_7d) & (sst_min_7d < 250 | sst_min_7d > 330), NA, sst_min_7d),
    sst_max_7d = if_else(!is.na(sst_max_7d) & (sst_max_7d < 250 | sst_max_7d >330), NA, sst_max_7d),
    depth = abs(depth), 
    distance.to.coast = if_else(distance.to.coast < 0, NA, distance.to.coast)
  )


library(e1071)
## Skewed var
detect_skewed_variables <- function(df, skew_threshold = 0.5) {
  
  # Garder uniquement les colonnes numériques
  num_vars <- df[sapply(df, is.numeric)]
  
  # Calculer la skewness pour chaque variable
  skewness_values <- sapply(num_vars, skewness, na.rm = TRUE)
  
  # Créer un data frame avec les résultats
  result <- data.frame(
    variable = names(skewness_values),
    skewness = skewness_values,
    asymetrie = ifelse(skewness_values > skew_threshold, "right-skewed",
                       ifelse(skewness_values < -skew_threshold, "left-skewed", "symétrique")),
    log_transform_recommande = skewness_values > 1
  )
  
  # Trier par skewness décroissante
  result <- result[order(-result$skewness), ]
  
  return(result)
}

cols <- c("chlorophyll_av", "chlorophyll_max", "chlorophyll_min", "gdp", "TotalGravity", "distance_to_port", "depth", "distance.to.coast")

if (all(cols %in% colnames(P_mat))) {
  P_mat <- P_mat %>%
    mutate(across(all_of(cols), ~ log10(. + 1), .names = "log_{.col}")) %>%
    select(-all_of(cols))  # Supprime les anciennes colonnes
}

library(mice)
library(Metrics)

## MF to fill the gaps
var_to_fill <- c("sst_av_7d", "sst_min_7d", "log_depth", "sst_max_7d", "log_gdp", "log_chlorophyll_av", "log_chlorophyll_max", "log_chlorophyll_min", "Province", "Realm", "log_distance.to.coast", "log_TotalGravity")
var_to_fill <- c("sst_av_7d", "sst_min_7d", "depth", "sst_max_7d", "gdp", "chlorophyll_av", "chlorophyll_max", "chlorophyll_min", "Province", "Realm", "distance.to.coast")

# Résultats stockés ici
imputation_results <- data.frame(variable = character(), r2 = numeric(), stringsAsFactors = FALSE)

# Pour chaque variable à évaluer
for (var_test in var_to_fill) {
  
  cat("➡️ Traitement de la variable :", var_test, "\n")
  
  # 1. Sélectionner les lignes où la variable n'est pas manquante
  non_missing_rows <- which(!is.na(P_mat[[var_test]]))
  
  # Continuer uniquement s’il y a assez de données
  if (length(non_missing_rows) >= 10) {
    
    # 2. Masquer aléatoirement une partie des valeurs (20%)
    set.seed(42)
    n_test <- round(length(non_missing_rows) * 0.2)
    test_rows <- sample(non_missing_rows, n_test)
    
    # 3. Masquer les valeurs dans une copie
    P_mat_masked <- P_mat
    P_mat_masked[test_rows, var_test] <- NA
    
    # 4. Imputer avec MICE + méthode Random Forest
    imp_rf <- mice(P_mat_masked[, var_to_fill], method = "rf", m = 3, maxit = 5, seed = 123, printFlag = FALSE)
    P_mat_completed <- complete(imp_rf)
    
    # 5. Comparer les valeurs imputées et réelles
    true_vals <- P_mat[[var_test]][test_rows]
    imputed_vals <- P_mat_completed[[var_test]][test_rows]
    
    # Vérifier que la variable est numérique
    if (is.numeric(true_vals)) {
      ss_res_test <- sum((imputed_vals - true_vals)^2)
      ss_tot_test <- sum((true_vals - mean(true_vals))^2)
      r2_test <- 1 - (ss_res_test / ss_tot_test)
      
      # Stocker les résultats
      imputation_results <- rbind(imputation_results, data.frame(variable = var_test, r2 = round(r2_test, 3)))
    } else {
      warning(paste("⚠️ La variable", var_test, "n'est pas numérique : évaluation ignorée."))
    }
  } else {
    warning(paste("⚠️ Trop peu de valeurs non manquantes pour la variable", var_test))
  }
}

# Afficher le tableau final
print(imputation_results)

save(imputation_results, file=here::here("Sea_of_immaturity", "outputs", "Spamm explicative model", "Mice_imputation_results.Rdata"))

Mice_Pmat <- mice(P_mat[, var_to_fill_log], 
            method = "rf",  
            m = 3,           # nombre d'imputations
            maxit = 10,      # nombre d'itérations
            seed = 123)

plot(Mice_Pmat)

P_mat_imputed <- complete(Mice_Pmat)

#Complete P_mat
P_mat <- P_mat %>%
  mutate( log_gdp = P_mat_imputed$log_gdp, sst_av_7d = P_mat_imputed$sst_av_7d, log_depth = P_mat_imputed$log_depth, sst_max_7d = P_mat_imputed$sst_max_7d, log_TotalGravity = P_mat_imputed$log_TotalGravity, 
        sst_min_7d = P_mat_imputed$sst_min_7d, log_chlorophyll_av = P_mat_imputed$log_chlorophyll_av, log_chlorophyll_max = P_mat_imputed$log_chlorophyll_max, log_chlorophyll_min=P_mat_imputed$log_chlorophyll_min, Province=P_mat_imputed$Province, Realm=P_mat_imputed$Realm, log_distance.to.coast=P_mat_imputed$log_distance.to.coast)

##-------------Diversité ontogénétique--------------------------

shannon_index <- Measurment_data_clustered %>%
  mutate(class_mat = cut(p_maturity,
                         breaks = c(0, 0.1, 0.9, 1),
                         labels = c("Juvenile", "Subadult", "Adult"))) %>%
  group_by(cluster, class_mat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(p = n / sum(n),
         H = -sum(p * log(p))) %>%
  summarise(shannon_index = first(H))

P_mat <- P_mat %>% left_join(shannon_index, by="cluster")

##------------Nb mature individuals---------------------------

Nb_mat_indiv <- Measurment_data_clustered %>%
  group_by(cluster) %>%
  filter(p_maturity>0.9, !is.na(cluster)) %>%
  summarise(Nb_mat_indiv = n(), .groups = "drop") 

##------------Nb mature IUCN individuals----------------------

load("~/Sea_of_immaturity/outputs/02_fill_species_traits/MF_species_traits_imputed.Rdata")

Measurment_data_clustered <- Measurment_data_clustered %>% left_join(species_traits_imputed %>% select(species_name, IUCN_inferred_Loiseau23), by="species_name")

Nb_mat_indiv_IUCN <- Measurment_data_clustered %>%
  group_by(cluster) %>%
  filter(p_maturity>0.9, !is.na(cluster), IUCN_inferred_Loiseau23=="Threatened") %>%
  summarise(Nb_mat_indiv_IUCN = n(), .groups = "drop") 

Presence_mat_indiv_IUCN <- Measurment_data_clustered %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster) %>%
  summarise(
    Threatened_mature_present = as.integer(any(p_maturity > 0.9 & IUCN_inferred_Loiseau23 == "Threatened")),
    .groups = "drop"
  )

Presence_mat_indiv_IUCN <- Presence_mat_indiv_IUCN %>% na.omit(Threatened_mature_present)

##----------------Complete P_mat with those data------------------

P_mat <- P_mat %>% 
  left_join(Nb_mat_indiv, by="cluster") %>%
  left_join(Nb_mat_indiv_IUCN, by="cluster") %>%
  left_join(Presence_mat_indiv_IUCN, by="cluster")


P_mat <- P_mat %>%
  mutate(
    Nb_mat_indiv = replace_na(Nb_mat_indiv, 0),
    Nb_mat_indiv_IUCN = replace_na(Nb_mat_indiv_IUCN, 0),
    Threatened_mature_present = replace_na(Threatened_mature_present, 0),
         )

#-----------------Distribution of the variables-------------------

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

save(P_mat, file=here::here("Sea_of_immaturity", "data", "derived_data", "P_mat.Rdata"))