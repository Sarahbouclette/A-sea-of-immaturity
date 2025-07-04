##---------------------P_maturity by site-------------------------------------

library(ggplot2)
library(dplyr)
library(plotly)
library(sf)

# load Measurment file with the individuals probability of beiin mature
load("~/Sea_of_immaturity/outputs/03_process_data/Measurment_w_pmat.Rdata")

##-------------Clustering--------------------------------------------------

library(dbscan)
library(geosphere)  # Pour calculer les distances géographiques

site_data <- Measurment_data %>%
  group_by(latitude, longitude, survey_date, Source, `Sample ID`) %>%
  summarise(
    mean_p_maturity = mean(p_maturity, na.rm = TRUE),
    nb_individuals = n()
  ) %>%
  tidyr::drop_na(latitude, longitude, survey_date)%>%
  ungroup()

# On veut regrouper les données des 5 caméras d'un rig BRUVS
# Filtrer les sources concernées
filtered_data <- site_data %>%
  filter(Source %in% c("bruvs_benthic", "bruvs_pelagic"), 
         !is.na(latitude), 
         !is.na(longitude))

# Paramètres du clustering
eps_distance <- 1 # 1 km de distance maximale
min_samples <- 1   # Permet de considérer tous les points

# Appliquer le clustering par survey_date

# Générer un identifiant unique pour chaque date
survey_dates <- unique(filtered_data$survey_date)
date_to_id <- setNames(seq_along(survey_dates), survey_dates)

clustered_data <- filtered_data %>%
  group_by(survey_date) %>%
  mutate(
    cluster_raw = dbscan::dbscan(as.matrix(cbind(latitude, longitude)), 
                                 eps = eps_distance / 111, minPts = min_samples)$cluster,
    cluster = ifelse(cluster_raw == 0, `Sample ID`, 
                     paste0("C", date_to_id[survey_date], "_", cluster_raw))  # Ajoute un préfixe unique par date
  ) %>%
  ungroup() %>%
  select(-cluster_raw)  # Optionnel : Supprimer la colonne brute des clusters


# Identifier les clusters trop grands (> 5 points)
clustered_data <- clustered_data %>%
  group_by(survey_date, cluster) %>%
  mutate(nb_samples = n()) %>%
  ungroup()

clusters_too_big <- clustered_data %>% filter(nb_samples > 5)

# Liste pour stocker les clusters raffinés
refined_clusters <- list()

# Appliquer un deuxième DBSCAN aux clusters trop grands
for (cl in unique(clusters_too_big$cluster)) {
  
  subset_data <- clustered_data %>% filter(cluster == cl)
  
  # Vérifier que le subset n'est pas vide (sécurité)
  if (nrow(subset_data) > 0) {
    
    # DBSCAN avec un eps plus strict
    sub_clusters <- dbscan::dbscan(as.matrix(cbind(subset_data$latitude, subset_data$longitude)), 
                                   eps = (eps_distance / 111) * 0.5,  # Réduction de la distance
                                   minPts = min_samples)$cluster
    
    # Vérifier qu'on a bien des sous-clusters valides
    if (any(sub_clusters > 0)) {
      subset_data <- subset_data %>%
        mutate(cluster = paste0(cl, "_", sub_clusters))
    }
    
    # Ajouter à la liste des clusters raffinés
    refined_clusters[[cl]] <- subset_data
  }
}

# Convertir la liste en data frame
refined_clusters <- bind_rows(refined_clusters) %>%
  group_by(survey_date, cluster) %>%
  mutate(nb_samples = n()) %>%
  ungroup()

# Fusionner avec les clusters corrects (< 5 points)
final_clusters <- clustered_data %>%
  filter(nb_samples <= 5) %>%
  bind_rows(refined_clusters) %>%
  filter(nb_samples <= 5) %>%
  select(-nb_samples)

# Ajouter les clusters à site_data
site_data <- site_data %>%
  left_join(final_clusters %>% select(`Sample ID`, cluster), 
            by = c("Sample ID")) %>%
  mutate(cluster = ifelse(is.na(cluster), `Sample ID`, cluster))

Measurment_data_clustered <- Measurment_data %>%
  left_join(site_data %>% dplyr::select(cluster, `Sample ID`), by = "Sample ID") %>%
  group_by(cluster) %>%
  mutate(
    longitude = mean(longitude, na.rm = TRUE),
    latitude = mean(latitude, na.rm = TRUE) 
  ) %>%
  ungroup()

save(Measurment_data_clustered, file=here::here("Sea_of_immaturity", "outputs", "03_process_data", "Measurment_data_clustered.Rdata"))

##-------------Carte des probas de maturité par site-------------------------
# Agréger par cluster si la communauté possède plus de 5 individus et calculer la moyenne des probabilités de maturité

 site_data <- Measurment_data_clustered %>%
  group_by(latitude, longitude, survey_date, Source, cluster) %>%
  summarise(
    mean_p_maturity = ifelse(n() >= 5, mean(p_maturity, na.rm = TRUE), NA), # Moyenne uniquement si n >= 10
    coeff_var_maturity = ifelse(n() >= 5, sd(p_maturity, na.rm = TRUE)/mean_p_maturity, NA),
    var_p_maturity = ifelse(sd(p_maturity) ==0, 1e-6, sd(p_maturity, na.rm = TRUE)),
    nb_individuals = n()
    #all_p_maturity = list(p_maturity)  # Stocker toutes les valeurs de probabilité par site
  ) %>%
  ungroup()

# Transformer en objet spatial sf pour ggplot
site_sf <- st_as_sf(site_data, coords = c("longitude", "latitude"), crs = 4326, remove=FALSE)

#Carte intéractive pour chaque source

library(leaflet)

leaflet(site_sf) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  
  # Ajouter des marqueurs pour chaque source
  addCircleMarkers(
    data = site_sf %>% filter(Source == "rls"),
    lng = ~longitude, lat = ~latitude,
    radius = ~mean_p_maturity * 10,
    color =~colorNumeric("viridis", domain = site_sf$mean_p_maturity)(mean_p_maturity), 
    fillOpacity = 0.8,
    group = "rls mean p_mat",
    popup = ~paste("Maturity probability:", round(mean_p_maturity, 3), "<br>",
                   "Survey date:", survey_date, "<br>",
                   "Number of individuals", nb_individuals, "<br>",
                   "Variance:", var_p_maturity, "<br>")
  ) %>%
  
  addCircleMarkers(
    data = site_sf %>% filter(Source == "bruvs_pelagic"),
    lng = ~longitude, lat = ~latitude,
    radius = ~mean_p_maturity* 10,
    color =~colorNumeric("viridis", domain = site_sf$mean_p_maturity)(mean_p_maturity),
    fillOpacity = 0.8,
    group = "bruvs_pelagic mean p_mat",
    popup = ~paste("Maturity probability:", round(mean_p_maturity, 3), "<br>",
                   "Survey date:", survey_date, "<br>",
                   "Number of individuals", nb_individuals, "<br>",
                   "Variance:", var_p_maturity)
  ) %>%
  
  addCircleMarkers(
    data = site_sf %>% filter(Source == "bruvs_benthic"),
    lng = ~longitude, lat = ~latitude,
    radius = ~mean_p_maturity * 10,
    color =~colorNumeric("viridis", domain = site_sf$mean_p_maturity)(mean_p_maturity),
    fillOpacity = 0.8,
    group = "bruvs_benthic mean p_mat",
    popup = ~paste("Maturity probability:", round(mean_p_maturity, 3), "<br>",
                   "Survey date:", survey_date, "<br>",
                   "Number of individuals", nb_individuals, "<br>",
                   "Variance:", var_p_maturity)
  ) %>%
  
  addCircleMarkers(
    data = site_sf %>% filter(Source == "rls"),
    lng = ~longitude, lat = ~latitude,
    radius = ~coeff_var_maturity ,
    color =~colorNumeric("magma", domain = site_sf$coeff_var_maturity)(coeff_var_maturity), 
    fillOpacity = 0.8,
    group = "rls coeff var",
    popup = ~paste("Variation coefficient:", round(coeff_var_maturity, 3), "<br>",
                   "Survey date:", survey_date)
  ) %>%
  
  addCircleMarkers(
    data = site_sf %>% filter(Source == "bruvs_pelagic"),
    lng = ~longitude, lat = ~latitude,
    radius = ~coeff_var_maturity,
    color =~colorNumeric("magma", domain = site_sf$coeff_var_maturity)(coeff_var_maturity),
    fillOpacity = 0.8,
    group = "bruvs_pelagic coeff var",
    popup = ~paste("Variation coefficient:", round(coeff_var_maturity, 3), "<br>",
                   "Survey date:", survey_date)
  ) %>%
  
  addCircleMarkers(
    data = site_sf %>% filter(Source == "bruvs_benthic"),
    lng = ~longitude, lat = ~latitude,
    radius = ~coeff_var_maturity,
    color =~colorNumeric("magma", domain = site_sf$coeff_var_maturity)(coeff_var_maturity),
    fillOpacity = 0.8,
    group = "bruvs_benthic coeff var",
    popup = ~paste("Variation coefficient:", round(coeff_var_maturity, 3), "<br>",
                   "Survey date:", survey_date)
  ) %>%
  
  # Ajouter la légende des couleurs
  addLegend(
    "bottomright",
    pal = colorNumeric("viridis", domain = site_sf$mean_p_maturity),
    values = site_sf$mean_p_maturity,
    title = "Mean probability of maturity by site",
    opacity = 1
  ) %>%
  addLegend(
    "bottomright",
    pal = colorNumeric("magma", domain = site_sf$coeff_var_maturity),
    values = site_sf$coeff_var_maturity,
    title = "Variation coefficient by site",
    opacity = 1
  ) %>%
  
  # Ajouter le contrôle interactif pour sélectionner les sources à afficher
  addLayersControl(
    overlayGroups = c("rls mean p_mat", "bruvs_pelagic mean p_mat", "bruvs_benthic mean p_mat", 
                      "rls coeff var", "bruvs_pelagic coeff var", "bruvs_benthic coeff var"),
    options = layersControlOptions(collapsed = FALSE)
  )

#--------------- Carte intéractive partielle-------------------

library(leaflet)
library(sf)
my_palette <- colorNumeric(
  palette = colorRampPalette(c("#440154", "#31688E", "#35B779", "#FFD700"))(50),
  domain = P_mat$Nb_mat_indiv_IUCN/P_mat$n.x)

# Transformer en objet spatial sf pour ggplot
P_mat_sf <- st_as_sf(P_mat, coords = c("longitude", "latitude"), crs = 4326, remove=FALSE)

#Carte intéractive pour chaque source

leaflet(P_mat_sf) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  
  # Ajouter des marqueurs pour chaque source
  addCircleMarkers(
    data = P_mat_sf,
    lng = ~longitude, lat = ~latitude,
    radius = ~P_mat$Nb_mat_indiv_IUCN/P_mat$n.x*10,
    color =~my_palette(P_mat$Nb_mat_indiv_IUCN/P_mat$n.x),
    fillOpacity = 0.8,
    popup = ~paste("Maturity probability:", round(mean_p_maturity, 3), "<br>",
                   "Survey date:", survey_date, "<br>",
                   "Variance:", var_p_maturity, "<br>",
                   "Nb_mature_indivs_IUCN", Nb_mat_indiv_IUCN)
  ) %>%
  
  addLegend(
    "bottomright",
    pal = my_palette,
    values = P_mat$Nb_mat_indiv_IUCN/P_mat$n.x,
    title = "Proportion of mature Threatened IUCN individuals at each site",
    opacity = 1
  )

##-------------Graph de distribution de proba de maturité par site-----------------

# Créer un graphique des distributions par site
plot_distributions <- plot_ly()

for (i in 1:nrow(site_data[1:200,])) {
  plot_distributions <- plot_distributions %>%
    add_trace(
      x = site_data$all_p_maturity[[i]],
      type = "histogram",
      histnorm = "density",
      name = paste("Site", i),
      visible = ifelse(i == 1, TRUE, "legendonly")
    )
}

plot_distributions <- plot_distributions %>%
  layout(
    title = "Distribution des probabilités de maturité par site",
    xaxis = list(title = "p_maturity"),
    yaxis = list(title = "Densité"),
    updatemenus = list(
      list(
        type = "dropdown",
        buttons = list(
          list(method = "relayout", args = list("title.text", "Probabilité moyenne de maturité par site"), label = "Moyenne"),
          list(method = "relayout", args = list("title.text", "Distribution des probabilités de maturité par site"), label = "Distribution")
        ),
        direction = "down",
        x = 1.1,
        y = 1.1
      )
    )
  )

# Afficher le graphique
print(plot_distributions)  # Distribution des probabilités


##-----------Maturity spectra---------------------
library(plotly)

# 1. Création des classes de probabilité de maturité
Abundance_data <- Measurment_data %>%
  mutate(p_maturity_class = cut(p_maturity, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  group_by(site_name, survey_date, p_maturity_class) %>%
  summarise(Abundance = sum(n(), na.rm = TRUE), .groups = "drop")

# 2. Calculer la moyenne des probabilités par classe pour les régressions
Abundance_data <- Abundance_data %>%
  mutate(p_maturity_mean = (as.numeric((sub("[\\[\\(](.+),(.+)\\]", "\\1", p_maturity_class))) + 
                              as.numeric(sub("[\\[\\(](.+),(.+)\\]", "\\2", p_maturity_class))) / 2)


# 3. Regrouper les données par site
Abundance_site_data <- Abundance_data %>%
  group_by(site_name, survey_date) %>%
  summarise(
    all_abundance = list(Abundance),
    all_p_maturity = list(p_maturity_mean)
  ) %>%
  ungroup()

# 4. Création du graphique interactif
plot_abundance <- plotly::plot_ly()

# 5. Ajouter la régression linéaire et les points pour chaque site
for (i in 1:min(nrow(Abundance_site_data), 200)) {  # Limite à 200 sites pour éviter surcharge
  site_df <- data.frame(
    p_maturity = Abundance_site_data$all_p_maturity[[i]],
    Abundance = Abundance_site_data$all_abundance[[i]]
  )
  
  # Ajuster une régression linéaire si assez de données
  if (nrow(site_df) > 1) {
    lm_model <- lm(Abundance ~ p_maturity, data = site_df)
    
    # Générer des prédictions
    pred_data <- data.frame(p_maturity = seq(min(site_df$p_maturity, na.rm = TRUE), 
                                             max(site_df$p_maturity, na.rm = TRUE), length.out = 100))
    pred_data$Abundance <- predict(lm_model, newdata = pred_data)
    
    # Ajouter les points
    plot_abundance <- plot_abundance %>%
      add_trace(
        x = site_df$p_maturity,
        y = site_df$Abundance,
        type = "scatter",
        mode = "markers",
        marker = list(size = 6, opacity = 0.6),
        name = paste("Site", Abundance_site_data$site_name[i])
      ) %>%
      add_trace(
        x = pred_data$p_maturity,
        y = pred_data$Abundance,
        type = "scatter",
        mode = "lines",
        name = paste("Régression Site", Abundance_site_data$site_name[i])
      )
  }
}

# 6. Personnalisation du layout
plot_abundance <- plot_abundance %>%
  layout(
    title = "Relation entre abondance et probabilité de maturité par site (classées)",
    xaxis = list(title = "Probabilité de maturité (classes)"),
    yaxis = list(title = "Abondance totale"),
    showlegend = TRUE
  )

# 7. Afficher le graphique
print(plot_abundance)


##-----------Visualisation des stats-----------

site_data <- Measurment_data_clustered %>%
  group_by(latitude, longitude, survey_date, Source, cluster) %>%
  summarise(
    mean_p_maturity = ifelse(n()>5, mean(p_maturity, na.rm = TRUE), NA),
    var_p_maturity = ifelse(n()>5, sd(p_maturity, na.rm = TRUE), NA),
    coeff_var_maturity = ifelse(n()>5,var_p_maturity/mean_p_maturity, NA),
    all_p_maturity = list(p_maturity),
    nb_individuals = n(),
    cumul_p_maturity = sum(p_maturity),
    cumul_p_mat_norm = cumul_p_maturity / nb_individuals
  ) %>%
  ungroup()

# Calculer les moyennes par source
site_data_summary <- site_data %>%
  group_by(Source) %>%
  summarise(
    mean_overall = mean(mean_p_maturity, na.rm = TRUE),
    cv_overall = mean(coeff_var_maturity, na.rm = TRUE)
  )

my_colors <- c(
  "bruvs_benthic" = "#1976D2",
  "rls" = "#F57C00",  
  "bruvs_pelagic" = "#64B5F6"  
)

## Créer l'histogramme des moyennes de p_mat par source
plot_distrib_mean <- ggplot(site_data, aes(x = mean_p_maturity, fill = Source)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7) +
  scale_fill_manual(values = my_colors) + 
  geom_vline(data = site_data_summary, aes(xintercept = mean_overall), 
             color = "red", linetype = "dashed", size = 1) +
  geom_text(data = site_data_summary, aes(x = mean_overall, y = Inf, 
                                          label = sprintf("Moyenne: %.2f", mean_overall)), 
            vjust = -0.5, color = "red", size = 5) +
  theme_minimal() +
  labs(title = "Distribution of average maturity probabilities by site for each source", 
       x = "Average maturity probability by site", y = "Density") +
  facet_wrap(~ Source, scales = "free")  # Séparer par source

# Afficher le graphique
print(plot_distrib_mean)


## Distribution des coeffs de variation de p_mat par source
plot_distrib_coeff_var <- ggplot(site_data, aes(x = coeff_var_maturity, fill = Source)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7) +
  scale_fill_manual(values = my_colors) + 
  geom_vline(data = site_data_summary, aes(xintercept = cv_overall), 
             color = "red", linetype = "dashed", size = 1) +
  geom_text(data = site_data_summary, aes(x = cv_overall, y = Inf, 
                                          label = sprintf("Moyenne: %.2f", mean_overall)), 
            vjust = -0.5, color = "red", size = 5) +
  theme_minimal() +
  labs(title = "Distribution of coefficients of variation of maturity probabilities by site for each source", 
       x = "Variatipn coefficient of the maturity probabilities by site", y = "Density") +
  facet_wrap(~ Source, scales = "free")  # Séparer par source

# Afficher le graphique
print(plot_distrib_coeff_var)


## Boxplot nb individuals measured by survey for each source
plot_boxplot_indiv <- ggplot(site_data, aes(x = Source, y = nb_individuals, fill = Source)) +
  geom_boxplot(color = "black", alpha = 0.7, outliers = FALSE) +
  scale_fill_manual(values = my_colors) + 
  theme_minimal() +
  labs(title = "Boxplots showing the number of individuals measured per site for each source", 
       x = "Sources", y = "Number of individuals measured by survey") +
  facet_wrap(~ Source, scales = "free")  # Séparer par source

# Afficher le graphique
print(plot_boxplot_indiv)

##-----------Quelles classes de probabilités sont les plus représentées?-------

library(ggplot2)
library(dplyr)

# Définition des classes de probabilité
bins <- c(0, 0.1, 0.5, 0.9, 1)
labels <- c("0-0.1", "0.1-0.5", "0.5-0.9", "0.9-1")

# Attribution de la classe à chaque individu
Cumul_measurment_data <- Measurment_data_clustered %>%
  mutate(prob_class = cut(p_maturity, breaks = bins, labels = labels, include.lowest = TRUE))

# Calcul du nombre de clusters ayant au moins un individu dans chaque classe
cluster_counts <- Cumul_measurment_data %>%
  group_by(prob_class) %>%
  summarise(nb_clusters = n_distinct(cluster), .groups = "drop", 
            nb_individuals = n()) %>%
  mutate(perc_clusters = (nb_clusters / length(unique(Measurment_data_clustered$cluster))*100))

# Création de l'histogramme
ggplot(cluster_counts, aes(x = prob_class, y = perc_clusters, fill = prob_class)) +
  geom_col() +
  scale_fill_manual(values = c("#D7191C", "#FDAE61", "#ABDDA4", "#2B83BA")) + 
  theme_minimal() +
  labs(
    title = "Pourcentage de sites présentant au moins un individus dans une classe probabilité de maturité considérée",
    x = "Classe de probabilité de maturité",
    y = "Pourcentage de sites"
  ) +
  theme(legend.position = "none")

##---------------------------Cumulative analysis---------------------
library(dplyr)

# Filtrer les individus matures (p_maturity >= 0.90)
mature_indiv <- Measurment_data_clustered %>%
  filter(p_maturity >= 0.90, !is.na(cluster), !is.na(Source))

# Compter le nombre d'individus matures par cluster et par source
cluster_counts_mat <- mature_indiv %>%
  group_by(Source, cluster) %>%
  summarise(nb_indiv = n(), .groups = "drop") %>%
  arrange(Source, desc(nb_indiv))

# Calculer la somme totale d'individus matures par source
cluster_counts_mat <- cluster_counts_mat %>%
  group_by(Source) %>%
  mutate(
    total_mature_indiv = sum(nb_indiv),
    cum_sum = cumsum(nb_indiv),
    cum_perc = cum_sum / total_mature_indiv,
    cluster_rank = row_number(),
    cluster_perc = cluster_rank / n()
  ) %>%
  ungroup()

# Identifier la proportion de clusters nécessaires pour atteindre 80%
clusters_needed_mat <- cluster_counts_mat %>%
  filter(cum_perc <= 0.80) %>%
  group_by(Source) %>%
  summarise(
    threshold_x = max(cluster_perc),
    .groups = "drop"
  )

# Fusionner avec les données pour le tracé
plot_data_mat <- cluster_counts_mat %>%
  left_join(clusters_needed_mat, by = "Source")

# Tracer la courbe cumulative avec un plot par source
ggplot(plot_data_mat, aes(x = cluster_perc, y = cum_perc)) +
  geom_line(color = "blue", size = 1) +
  geom_point(aes(color = cum_perc >= 0.80), size = 1.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = threshold_x), linetype = "dashed", color = "red") +
  facet_wrap(~Source, scales = "free_x") +
  labs(
    title = "Répartition des individus matures (p_maturity ≥ 0.90) par cluster et par source",
    x = "Proportion de clusters (classés du + au - denses)",
    y = "Proportion cumulative des individus matures",
    color = "≥ 80% atteint"
  ) +
  theme_minimal()

# Same for immature indiv
immature_indiv <- Measurment_data_clustered %>%
  filter(p_maturity < 0.1)

# Compter le nombre d'individus par cluster
cluster_counts_immat <- immature_indiv %>%
  group_by(cluster) %>%
  summarise(nb_indiv = n(), .groups = "drop") %>%
  arrange(desc(nb_indiv))  # Trier du plus grand au plus petit

cluster_counts_immat <- cluster_counts_immat[!is.na(cluster_counts_immat$cluster),]

# Calculer le nombre total d'individus matures
total_immature_indiv <- sum(cluster_counts_immat$nb_indiv)

# Déterminer combien de clusters contiennent 80% des individus matures
cluster_counts_immat <- cluster_counts_immat %>%
  mutate(cum_sum = cumsum(nb_indiv),  # Somme cumulative des individus
         cum_perc = cum_sum / total_immature_indiv)  # Pourcentage cumulé

# Trouver le nombre minimum de clusters nécessaires pour atteindre 80%
clusters_needed_immat <- cluster_counts_immat %>%
  filter(cum_perc <= 0.80) %>%
  nrow()/nrow(cluster_counts_immat)

# Créer un dataframe pour la visualisation
plot_data_immat <- cluster_counts_immat %>%
  mutate(cluster_rank = row_number(), 
         cluster_perc = cluster_rank/n())  # Rang des clusters

# Tracer la courbe cumulative
ggplot(plot_data_immat, aes(x = cluster_perc, y = cum_perc)) +
  geom_line(color = "blue", size = 1) +
  geom_point(aes(color = cum_perc >= 0.80), size = 2) +  # Mettre en évidence le seuil 80%
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +  # Ligne à 80%
  geom_vline(xintercept = clusters_needed_immat, linetype = "dashed", color = "red") +  # Marqueur du seuil
  labs(
    title = "Répartition des individus immatures (p_maturity < 0.10) par cluster",
    x = "Nombre de clusters (classés du + au - denses)",
    y = "Proportion cumulative des individus immatures",
    color = "Seuil 80%"
  ) +
  theme_minimal()

Perc_indiv_mature <- Measurment_data_clustered %>%
  group_by(cluster, Source) %>%
  summarise(
    nb_indiv_mature = sum(p_maturity >= 0.90),   # Nombre d'individus matures
    total_indiv = n(),                           # Total d'individus dans le cluster
    perc_mature = nb_indiv_mature / total_indiv  # Proportion d'individus matures
  ) %>%
  ungroup()

plot_perc_mature_indiv <- ggplot(Perc_indiv_mature, aes(x = perc_mature, fill = Source)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7) +
  scale_fill_manual(values = my_colors)  +
  theme_minimal() +
  labs(title = "Distribution des sites selon le pourcentage d'individus matures qu'ils présentent", 
       x = "Pourcentage d'individus matures (p_maturity >= 0.90)", y = "Density = nombre de sites") +
  facet_wrap(~ Source, scales = "free")  # Séparer par source

Perc_indiv_immature <- Measurment_data_clustered %>%
  group_by(cluster, Source) %>%
  summarise(
    nb_indiv_immature = sum(p_maturity < 0.10),   # Nombre d'individus immatures
    total_indiv = n(),                           # Total d'individus dans le cluster
    perc_immature = nb_indiv_immature / total_indiv  # Proportion d'individus immatures
  ) %>%
  ungroup()

plot_perc_immature_indiv <- ggplot(Perc_indiv_immature, aes(x = perc_mature, fill = Source)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7) +
  scale_fill_manual(values = my_colors)  +
  theme_minimal() +
  labs(title = "Distribution des sites selon le pourcentage d'individus immatures qu'ils présentent", 
       x = "Pourcentage d'individus immatures (p_maturity < 0.10)", y = "Density = nombre de sites") +
  facet_wrap(~ Source, scales = "free")  # Séparer par source

  ## Courbe par in ou out mpa

# Ajouter une variable pour l'état MPA
mature_indiv <- Measurment_data_clustered %>%
  filter(p_maturity >= 0.90, !is.na(cluster), !is.na(Source), !is.na(inside_mpa)) %>%
  mutate(MPA_status = protection_status)

# Compter le nombre d'individus matures par cluster, par source et par MPA
cluster_counts_mat <- mature_indiv %>%
  group_by(MPA_status, cluster) %>%
  summarise(nb_indiv = n(), .groups = "drop") %>%
  arrange(MPA_status, desc(nb_indiv))

# Calculer la somme totale d'individus matures par source et MPA
cluster_counts_mat <- cluster_counts_mat %>%
  group_by(MPA_status) %>%
  mutate(
    total_mature_indiv = sum(nb_indiv),
    cum_sum = cumsum(nb_indiv),
    cum_perc = cum_sum / total_mature_indiv,
    cluster_rank = row_number(),
    cluster_perc = cluster_rank / n()
  ) %>%
  ungroup()

# Identifier la proportion de clusters nécessaires pour atteindre 80%
clusters_needed_mat <- cluster_counts_mat %>%
  filter(cum_perc <= 0.80) %>%
  group_by(MPA_status) %>%
  summarise(
    threshold_x = max(cluster_perc),
    .groups = "drop"
  )

# Fusionner pour le tracé
plot_data_mat <- cluster_counts_mat %>%
  left_join(clusters_needed_mat, by = c("MPA_status"))

# Tracer la courbe cumulative avec un plot par source et statut MPA
ggplot(plot_data_mat, aes(x = cluster_perc, y = cum_perc)) +
  geom_line(color = "blue", size = 1) +
  geom_point(aes(color = cum_perc >= 0.80), size = 1.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = threshold_x), linetype = "dashed", color = "red") +
  facet_wrap(~ MPA_status, scales = "free_x") +
  labs(
    title = "Répartition des individus matures (p_maturity ≥ 0.90) par cluster")

##----------------Effet de la réserve-----------------
library(ggpubr)

comparisons <- list(
  c("out", "low"),
  c("out", "medium"),
  c("out", "high"), 
  c("low", "medium"),
  c("low", "high"),
  c("medium","high")
)

comparisons <- list(
  c("fishing", "no.take"),
  c("no.take", "restricted.take"),
  c("fishing", "restricted.take")
)

# Création du boxplot avec les comparaisons paire par paire
plot_protect_status <- ggplot(site_data[!is.na(site_data$protection_status),], 
                              aes(x = as.factor(protection_status), y = shannon_index, fill= protection_status)) +
  geom_boxplot(alpha = 0.7) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif", size = 4) +  # Comparaisons paire par paire
  theme_minimal() +
  labs(title = "Comparison of maturity probabilities according to the protection status",
       x = "MPA protection status", y = "p_maturity mean by site") +
  scale_fill_brewer(palette = "Set2")

plot_protect_status

##----------------------Diversité ontogénétique------------------

library(dplyr)

# Calcul de l'index de shannon sur les diversité de stade de vie
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

P_mat <- P_mat %>%
  mutate(site_class = case_when(
    mean_p_maturity < 0.1 & cv_p_maturity < 1.75 ~ "Immature",
    mean_p_maturity > 0.9 & cv_p_maturity < 1 ~ "Mature",
    cv_p_maturity > 1.75 ~ "Highly Heterogeneous",
    between(mean_p_maturity, 0.1, 0.9) & cv_p_maturity <= 1.5 ~ "Intermediate",
    TRUE ~ "Uncertain"
  ))

# Test de la corrélation avec cv_p_mat 
# 1 . Visualisation
ggplot(P_mat, aes(x = shannon_index, y = cv_p_maturity, color=site_class)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C") +
  labs(
    x = "Shannon Index (Ontogenetic diversity)",
    y = "CV of maturity probability",
    title = "Relationship between ontogenetic diversity and maturity variability"
  ) +
  theme_minimal()

# 2. Tester la corrélation (de base Pearson, ou Spearman si non linéaire)
cor_test_result <- cor.test(P_mat$shannon_index, P_mat$cv_p_maturity, method = "pearson")

# 3. Afficher le résultat
cor_test_result

