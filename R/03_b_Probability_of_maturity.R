pkgs <- c("here", "dplyr", "ggplot2",
          "pbmcapply", "patchwork", "brms",
          "tibble", "stringr", "hrbrthemes", "pheatmap", "tiblr", "rfishbase")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


load("~/Sea_of_immaturity/outputs/03_process_data/data_bruvs_rls.Rdata")
load("~/Sea_of_immaturity/outputs/03_process_data/metadata_bruvs_rls.Rdata")
source(here::here("Sea_of_immaturity","R", "01_a_check_scientific_names.R"))
load("~/Sea_of_immaturity/data/derived_data/fishbase_sp_names.Rdata")
load("~/Sea_of_immaturity/outputs/02_fill_species_traits/RF_LengthMat_predictions.Rdata")

Measurment_data <- Measurment_data %>%
  select(-Genus) %>%
  rename(species_name = Binomial) 

sp_list <- sp_list %>% select(species_name, spec_code, fishbase_name) %>% distinct

Measurment_data <- Measurment_data %>% left_join(sp_list, by = "species_name")
  
Measurment_data <- Measurment_data[Measurment_data$fishbase_name %in% LengthMat_predicted$fishbase_name ,] #Keep species we know Lmat

# Add location info
Measurment_data <- Measurment_data %>%
  left_join(geographic_data %>% select(-Expedition, -Source), by = "Sample ID") %>%
  left_join(LengthMat_predicted %>% select(Family, Genus, Order, Class, fishbase_name), by = "fishbase_name") 

#-------------Formatage des données pour la régression---------------

Measurment_data <- Measurment_data %>%
  left_join(LengthMat_predicted %>% select(fishbase_name, Lm_min, Lm_max), by = "fishbase_name")

#paramètres
epsilon <- 1e-4 # Petite valeur pour éviter p_mat=1
p_min <- 0.1 # Probabilité maximale avant Lm_min
p_max <- 0.9
delta <- 0.01  # 1% de Lm_min, ajustable si nécessaire

# Fonction pour calculer beta assurant p_maturity(Lm_max) = 0.9
calculate_beta <- function(Lm_min, Lm_max, p_min, p_max) {
  log((p_max / (1 - p_max)) / (p_min / (1 - p_min))) / 0.5
}

Measurment_data <- Measurment_data %>%
  group_by(species_name) %>%
  mutate(
    # Vérifie si Lm_min == Lm_max avant de modifier
    same_Lm = (Lm_min == Lm_max),
    
    # Élargissement symétrique si nécessaire
    Lm_max = ifelse(same_Lm, Lm_min * (1 + delta), Lm_max),
    Lm_min = ifelse(same_Lm, Lm_min * (1 - delta), Lm_min),
    
    #define beta
    beta = calculate_beta(Lm_min, Lm_max, p_min, p_max),
    
    # Calcul de la probabilité de maturité avec la correction
    p_maturity = case_when(
      Length < Lm_min ~ p_min * (Length / Lm_min),  
      Length >= Lm_min & Length <= Lm_max ~ p_min + (1 - p_min) / (1 + exp(-beta * ((Length - Lm_min) / (Lm_max - Lm_min) - 0.5))),  
      Length > Lm_max ~ 1 - epsilon - (1 - (p_min + (1 - p_min) / (1 + exp(-beta * (0.5))))) * exp(-beta * (Length - Lm_max))
    )
  ) %>%
  ungroup()

save(Measurment_data, file=here::here("Sea_of_immaturity", "outputs", "03_process_data", "Measurment_w_pmat"))


##----------------------Visualisation--------------------------------------

library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)

# Calculer les distributions de probabilités
plot_data <- Measurment_data %>%
  group_by(Family, fishbase_name) %>%
  slice(1) %>%  # Étendre sur toute la plage de longueurs
  ungroup()

# Créer un graphique interactif avec histogrammes des distributions
plot <- plot_ly()

# Ajouter une trace par famille sous forme d'histogramme
for (fam in unique(plot_data$Family)) {
  subset_data <- plot_data %>% filter(Family == fam)
  
  plot <- plot %>%
    add_trace(
      x = subset_data$p_maturity,
      type = "histogram",
      histnorm = "density",
      name = fam,
      visible = ifelse(fam == unique(plot_data$Family)[1], TRUE, "legendonly")
    )
}

# Mise en forme du graphique
plot <- plot %>%
  layout(
    title = "Distribution des probabilités de maturité par famille",
    xaxis = list(title = "Probabilité de maturité"),
    yaxis = list(title = "Densité"),
    updatemenus = list(
      list(
        type = "dropdown",
        buttons = lapply(unique(plot_data$Family), function(fam) {
          list(
            method = "update",
            args = list("visible", plot_data$Family == fam),
            label = fam
          )
        }),
        direction = "down",
        x = 1.1,
        y = 1.1
      )
    )
  )

plot

# Vérification : histogramme des probabilités attribuées
hist(data_sampled$p_maturity, breaks = 50, main = "Distribution des probabilités de maturité", xlab = "p_maturity")

##-------------------Visualisation de la prior choisie---------------

# Paramètres de l'espèce test
Lm_min <- 30
Lm_max <- 40
p_min <- 0.1
p_max <- 0.9
epsilon <- 1e-4
delta <- 0.01

# Générer des longueurs de poissons
Length <- seq(20, 60, by = 0.1)

# Fonction pour calculer beta assurant p_maturity(Lm_max) = 0.9 et p_maturity(Lm_min)=0.1

calculate_beta <- function(Lm_min, Lm_max, p_min, p_max) {
  log((p_max / (1 - p_max)) / (p_min / (1 - p_min))) / 0.5
}
beta <- calculate_beta(Lm_min, Lm_max, p_min, p_max)


# Calcul de la probabilité de maturité avec la correction
p_maturity <- case_when(
  Length < Lm_min ~ p_min * (Length / Lm_min),  
  Length >= Lm_min & Length <= Lm_max ~ p_min + (1 - p_min) / (1 + exp(-beta * ((Length - Lm_min) / (Lm_max - Lm_min) - 0.5))),  
  Length > Lm_max ~ 1 - epsilon - (1 - (p_min + (1 - p_min) / (1 + exp(-beta * (0.5))))) * exp(-beta * (Length - Lm_max))
)

# Créer un data frame
df <- data.frame(Length, p_maturity) 

# Tracer la courbe
plot <- ggplot(df, aes(x = Length, y = p_maturity)) +
  geom_line(color = "blue", size = 1) +
  geom_vline(xintercept = Lm_min, linetype = "dashed", color = "red") +
  geom_vline(xintercept = Lm_max, linetype = "dashed", color = "red") +
  labs(x = "Length (L)", y = "Maturity probability", 
       title = "Maturity probability curve as a function of length") +
  theme_minimal()



