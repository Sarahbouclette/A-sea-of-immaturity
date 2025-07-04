################################################################################
##
##  Uses trophic guilds inferred by Parravicini et al. 2020 to extract trophic 
##   guild of each RLS species.
##
## 1b_other_traits_troph_morpho_range_vul.R
##
## 01/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "gtools")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

##-------------loading data-------------
# Extract trophic guild of tropical fishes inferred by Parravicini et al. 2020
#  (DOI: https://doi.or g/10.1371/journal.pbio.3000702) 
#  -> All species in S3 Table (https://doi.org/10.1371/journal.pbio.3000702.s004)
trophic_guilds <- read.csv(here::here("Sea_of_immaturity","data", "raw_data", "trophic_guilds_parravicini_2020_STable3.csv"))
trophic_guilds$species <- gsub("_", " ", trophic_guilds$species)

# Morphological data (from Friedman et al. 2021)
load(file = here::here("Sea_of_immaturity","data", "raw_data", "fishshape.RData"))

# Bathymetric and geographic range from Duhamet et al. 2023: https://doi.org/10.6084/m9.figshare.20403111
load(file = here::here("Sea_of_immaturity","data", "raw_data", "bathymetric_range_Duhamet2023", "Actinopterygii.Rdata"))
load(file = here::here("Sea_of_immaturity","data", "raw_data", "bathymetric_range_Duhamet2023", "Chondrichthyes.Rdata"))
Mat_Pa_actino <- readRDS(file = here::here("Sea_of_immaturity","data", "raw_data", "bathymetric_range_Duhamet2023", "Mat_Pa_teleo.RDS"))
Mat_Pa_chond <- readRDS(file = here::here("Sea_of_immaturity","data", "raw_data", "bathymetric_range_Duhamet2023", "Mat_Pa_chond.RDS"))

# Climatic Risk from Boyce et al. (2022), Dryad, Dataset, https://doi.org/10.5061/dryad.7wm37pvwr
clim_risk <- read.csv(file = here::here("Sea_of_immaturity","data","raw_data", "Boyce_etal_2022_NATCC", "Boyce_etal_2022_NATCC.csv"))

# Species home Range from Mouquet et al. 2023, Github, https://github.com/nmouquet/RLS_HUM_INT
gbif005 <- read.csv(here::here("Sea_of_immaturity","data","raw_data","species_range_Mouquet_2023" ,"species_range_size_005.csv"))
gbif01 <- read.csv(here::here("Sea_of_immaturity","data","raw_data", "species_range_Mouquet_2023","species_range_size_01.csv"))

# IUCN categories infered by Loiseau et al. 2023
load(file = here::here("Sea_of_immaturity","data", "raw_data", "all_predict_IUCN_Loiseau_2023.Rdata"))

# Length characteristics repertoried by Chen et al. 2021
Lengths_Chen <- read.csv(file = here::here("Sea_of_immaturity","data","raw_data", "Lm_Lmax_Chen_et_al2021.csv"))

# Length characteristics repertoried by Chu et Pauly 2021
Lengths_Pauly <- read.csv(file = here::here("Sea_of_immaturity","data","raw_data", "Lm_Lmax_chu_Pauly2021.csv"))

# LengthMax given by Jessica
Lengths_Jessica <- read_xlsx(here::here("Sea_of_immaturity","data","raw_data", "BRUVS taxa list 2025_01_6.xlsx"), sheet = 3)

# Other traits given by Jessica
Other_traits_Jessica <- read_xlsx(here::here("Sea_of_immaturity","data","raw_data", "BRUVS taxa list 2025_01_6.xlsx"), sheet = 2)

# Traits by Tsikliras
Lengths_Tsikliras <- read.csv(file = here::here("Sea_of_immaturity","data","raw_data", "Length_at_maturity_Tsikliras.csv"))

# Fonction pour nettoyer les valeurs numériques
clean_numeric_column <- function(col) {
  if (is.list(col)) {
    sapply(col, function(x) mean(as.numeric(x), na.rm = TRUE))
  } else {
    sapply(col, function(x) {
      x <- as.character(x)
      if (grepl(",", x, fixed = TRUE)) {
        mean(as.numeric(unlist(strsplit(x, ","))), na.rm = TRUE)
      } else {
        as.numeric(x)
      }
    })
  }
}

##-------------trophic guild -------------
sp_troph <- trophic_guilds$species
sp_rls <- species_traits_fishbase$species_name
sp_rls_fb <- species_traits_fishbase$fishbase_name
# 
length(which(sp_rls %in% sp_troph))
length(which(sp_rls_fb %in% sp_troph))

species_trophic_guild <- trophic_guilds |> 
  dplyr::select(species, trophic_guild_predicted_text) |> 
  dplyr::rename(trophic_guild = trophic_guild_predicted_text,
                species_name = species)


##-------------Fishshape (Morphology) -------------
# colnames(fishshape)
# unique(fishshape$method)
# unique(fishshape$taxo_scale)
# unique(fishshape$trait_name)
# unique(fishshape$unit)
# unique(fishshape$number_ind_measured)
#unique(fishshape$life_stage)
# summary(fishshape$trait_value)
# summary(fishshape$spec_code)

shape <- tidyr::pivot_wider(fishshape,
                            names_from = trait_name,
                            values_from = trait_value)

## Mean morphometric traits at the species scale
shape <- shape[-which(is.na(shape$standard_length)),]

shape_ratio <- shape |> 
  dplyr::mutate(species_names = gsub("_", " ", species_names),
                head_depth_ratio = head_depth / standard_length,
                lower_jaw_ratio = lower_jaw_length / standard_length,
                body_depth_ratio = max_body_depth / standard_length,
                body_width_ratio = max_fish_width / standard_length,
                caudalpeduncle_depth_ratio = min_caudalpeduncle_depth / standard_length,
                caudalpeduncle_width_ratio = min_caudalpeduncle_width / standard_length,
                mouth_width_ratio = mouth_width / standard_length) |> 
  dplyr::group_by(species_names) |> 
  dplyr::mutate(dplyr::across(head_depth_ratio:mouth_width_ratio,
                              mean)) |> 
  dplyr::select(-c(specimen_id:fieldsources, check:total_weight)) |> 
  unique() |> 
  dplyr::ungroup() |> 
  dplyr::rename(species_name = species_names)



##-------------Clean IUCN inference (Loiseau 2023)-------------
IUCN_inferred <- dat_network |> 
  dplyr::mutate(species_name = gsub("_", " ", species)) |> 
  dplyr::rename(IUCN_inferred_Loiseau23 = IUCN_final) |> #only for teleost species
  dplyr::select(species_name, IUCN_inferred_Loiseau23)



##-------------Clean bathymetric range (Duhamet 2023, from Albouy 2019)-------------
bathymetry <- gtools::smartbind(Actinopterygii, Chondrichthyes) |> 
  dplyr::select(Species:Depth_max) |> 
  dplyr::rename(species_name = Species, depth_min = Depth_min, depth_max = Depth_max)

summary(bathymetry)



##-------------Clean Climatic risk (Boyce 2023)-------------
table(clim_risk$Phylum)
clim_risk <- clim_risk[clim_risk$Phylum %in% "Chordata",]

ClimVuln_SSP126 <- clim_risk |> 
  dplyr::filter(Experiment == "SSP126") |> 
  dplyr::select(SPname, ClimVuln) |> 
  dplyr::rename(species_name = SPname, ClimVuln_SSP126 = ClimVuln)

ClimVuln_SS5856 <- clim_risk |> 
  dplyr::filter(Experiment == "SSP585") |> 
  dplyr::select(SPname, ClimVuln) |> 
  dplyr::rename(species_name = SPname, ClimVuln_SSP585 = ClimVuln)


ClimVuln <- dplyr::full_join(ClimVuln_SSP126, ClimVuln_SS5856)

summary(ClimVuln)
plot(ClimVuln$ClimVuln_SSP126~ClimVuln$ClimVuln_SSP585)
summary(lm(ClimVuln$ClimVuln_SSP126~ClimVuln$ClimVuln_SSP585)) #r-squared = 0.87 , estimate = 1.00



##-------------Clean geographic range-------------
## Mouquet et al. 2023
home_range <- dplyr::full_join(gbif01, gbif005) |> 
  dplyr::select(gbif_valid_name, GBIF_NC_n_cells_005, GBIF_NC_n_cells_01) |> 
  dplyr::rename(species_name = gbif_valid_name,
                range_n_cells_005 = GBIF_NC_n_cells_005, 
                range_n_cells_01 = GBIF_NC_n_cells_01)

summary(home_range)

## Duhamet et al. 2023
actino_range <- colSums(Mat_Pa_actino[, -c(1,2)])
chond_range <- colSums(Mat_Pa_chond[, -c(1,2, which(is.na(colnames(Mat_Pa_chond))))])

all_range <- data.frame(species_name = names(c(actino_range, chond_range)),
                        geographic_range_Albouy19 = c(actino_range, chond_range)) |> 
  
  dplyr::full_join(home_range) |> 
  dplyr::mutate(species_name = gsub("_", " ", species_name))


plot(log(all_range$range_n_cells_01) ~ log(all_range$range_n_cells_005))
summary(lm(log(home_range$range_n_cells_01) ~ log(home_range$range_n_cells_005))) #r-squared = 0.995 , estimate = 0.95

plot(log(all_range$geographic_range_Albouy19) ~ log(all_range$range_n_cells_005))

##-------------Clean Lengths from Chen et al. 2021---------------

Lengths_Chen <- Lengths_Chen %>%
  dplyr::select(Scientific.name, Sex, Lm..mm., Lmax..mm., Wmax..g.) %>%
  dplyr::rename(species_name = Scientific.name, LengthMatMin = Lm..mm., LengthMax = Lmax..mm., WeightMax = Wmax..g.) %>%
  dplyr::mutate(LengthMatMin = LengthMatMin/10, LengthMax  = LengthMax/10)

Lengths_Chen$LengthMax <- abs(Lengths_Chen$LengthMax)

# Harmoniser les valeurs de la colonne Sex
Lengths_Chen <- Lengths_Chen %>%
  mutate(Sex = case_when(
    Sex %in% c('U/M', 'U') ~ 'unsexed',
    Sex == 'F' ~ 'female',
    Sex == 'M' ~ 'male',
    TRUE ~ Sex  # Conserver les valeurs inchangées si non listées
  ))

# Pivot des longueurs
Lengths_Chen_pivot1 <- Lengths_Chen %>%
  pivot_wider(
    names_from = Sex,
    values_from = LengthMatMin,
    names_prefix = "LengthMatMin_"
  )%>%
  mutate(across(starts_with("LengthMatMin_"), clean_numeric_column)) %>%
  group_by(species_name) %>%   summarise(across(starts_with("LengthMatMin_"), ~ mean(., na.rm = TRUE)),
  .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))  # Remplacer NaN par NA
  
Lengths_Chen_pivot2 <- Lengths_Chen %>%
  pivot_wider(
    names_from = Sex,
    values_from = LengthMax,
    names_prefix = "LengthMax_"
  ) %>%
  mutate(across(starts_with("LengthMax_"), clean_numeric_column)) %>%
  group_by(species_name) %>%
  summarise(across(starts_with("LengthMax_"), ~ mean(., na.rm = TRUE)), 
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))  # Remplacer NaN par NA

Lengths_Chen_final<-left_join(Lengths_Chen_pivot1, Lengths_Chen_pivot2, by = 'species_name')
#save(Lengths_Chen_final, file = here::here("Documents", "Sea_of_immaturity", "data", "raw_data", "Lengths_Chen_bysex.Rdata"))


##-------------Clean Lengths from Chu et Pauly 2021---------------

Lengths_Pauly <- Lengths_Pauly %>%
  dplyr::select(Species, Sex, Lm..cm., Lmax..cm., Wmax..g.) %>%
  dplyr::rename(species_name = Species, LengthMatMin = Lm..cm., LengthMax = Lmax..cm., WeightMax = Wmax..g.)


Lengths_Pauly <- Lengths_Pauly %>%
  mutate(Sex = case_when(
    Sex %in% c('U', 'UI') ~ 'unsexed',
    Sex == 'F' ~ 'female',
    Sex == 'M' ~ 'male',
    TRUE ~ Sex  # Conserver les valeurs inchangées si non listées
  ))

# Pivot des longueurs
Lengths_Pauly_pivot1 <- Lengths_Pauly %>%
  pivot_wider(
    names_from = Sex,
    values_from = LengthMatMin,
    names_prefix = "LengthMatMin_"
  )%>%
  mutate(across(starts_with("LengthMatMin_"), clean_numeric_column)) %>%
  group_by(species_name) %>%   summarise(across(starts_with("LengthMatMin_"), ~ mean(., na.rm = TRUE)),
                                         .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))  # Remplacer NaN par NA

Lengths_Pauly_pivot2 <- Lengths_Pauly %>%
  pivot_wider(
    names_from = Sex,
    values_from = LengthMax,
    names_prefix = "LengthMax_"
  ) %>%
  mutate(across(starts_with("LengthMax_"), clean_numeric_column)) %>%
  group_by(species_name) %>%
  summarise(across(starts_with("LengthMax_"), ~ mean(., na.rm = TRUE)), 
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))  # Remplacer NaN par NA

Lengths_Pauly_final<-left_join(Lengths_Pauly_pivot1, Lengths_Pauly_pivot2, by = 'species_name')
#save(Lengths_Pauly_final, file = here::here("Documents", "Sea_of_immaturity", "data", "raw_data", "Lengths_Pauly_bysex.Rdata"))

##------------Clean Lengths given by Jessica----------------------------

Lengths_Jessica <- Lengths_Jessica %>%
  dplyr::rename(spec_code = SpecCode, species_name = Species, LengthMax_unsexed = Length)

Lengths_Jessica$LengthMax_unsexed <- as.double(Lengths_Jessica$LengthMax_unsexed)

##------------Clean other traits given by Jessica------------------------

Other_traits_Jessica <- Other_traits_Jessica %>%
  dplyr::select(Taxa, `Family +`, Group, Level, `TR (Troph)`, PD, `VUL-F`, `VUL-C`, IUCN, a, b, FLest ) %>%
  dplyr::rename(Family = `Family +`, Troph = `TR (Troph)`, PhyloDiv = PD, ForkLength = FLest, species_name = Taxa)

#-----------------Add maturity data from Tsikliras--------------------------

Lengths_Tsikliras <- Lengths_Tsikliras %>%
  rename(species_name = Species, LengthMax = Lmax, LengthMatMin = Lm) %>%
  select(species_name, LengthMax, LengthMatMin, Sex) %>%
  filter(Sex != "") %>%
  mutate(species_name = sub("\\*{1,3}$", "", species_name))  # Supprimer * ou ** ou *** à la fin

Lengths_Tsikliras[Lengths_Tsikliras$Sex%in%c('C', ' C', ' c'),]$Sex <-'unsexed'
Lengths_Tsikliras[Lengths_Tsikliras$Sex%in%c('F'),]$Sex <-'female'
Lengths_Tsikliras[Lengths_Tsikliras$Sex%in%c('M'),]$Sex <-'male'

#Pivot des longueurs
Lengths_Tsikliras_pivot1 <- Lengths_Tsikliras %>%
  pivot_wider(
    names_from = Sex,
    values_from = LengthMatMin,
    names_prefix = "LengthMatMin_"
  )%>%
  mutate(across(starts_with("LengthMatMin_"), clean_numeric_column)) %>%
  group_by(species_name) %>%   summarise(across(starts_with("LengthMatMin_"), ~ mean(., na.rm = TRUE)),
                                         .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))  # Remplacer NaN par NA

Lengths_Tsikliras_pivot2 <- Lengths_Tsikliras %>%
  pivot_wider(
    names_from = Sex,
    values_from = LengthMax,
    names_prefix = "LengthMax_"
  ) %>%
  mutate(across(starts_with("LengthMax_"), clean_numeric_column)) %>%
  group_by(species_name) %>%
  summarise(across(starts_with("LengthMax_"), ~ mean(., na.rm = TRUE)), 
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))  # Remplacer NaN par NA

Lengths_Tsikliras_final<-left_join(Lengths_Tsikliras_pivot1, Lengths_Tsikliras_pivot2, by = 'species_name')
#save(Lengths_Tsikliras_final, file = here::here("Documents", "Sea_of_immaturity", "data", "raw_data", "Lengths_Tsikliras_bysex.Rdata"))

##-------------Join all the data (without dealing with names)-------------

other_traits_raw <- species_trophic_guild |> 
  dplyr::full_join(shape_ratio) |> 
  dplyr::full_join(IUCN_inferred) |> 
  dplyr::full_join(bathymetry) |> 
  dplyr::full_join(ClimVuln) |> 
  dplyr::full_join(Lengths_Chen_final, copy=TRUE) |>
  dplyr::full_join(Lengths_Pauly_final, copy=TRUE) |>
  dplyr::full_join(Lengths_Jessica) |>
  dplyr::full_join(Other_traits_Jessica)|>
  dplyr::full_join(Lengths_Tsikliras_final)


##-------------Save data-------------
save(other_traits_raw, file = here::here("Sea_of_immaturity", "data", "derived_data", "01_c_other_traits_raw.Rdata"))

