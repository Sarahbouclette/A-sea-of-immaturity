################################################################################
##
##  Takes the Species x traits matrices from fishbase and other sources,
##   check all the species names, and merge the tables.
##
## 1c_merge_species_traits.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "gtools", "funbiogeo")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

 remotes::install_github("FRBCesab/funbiogeo")
##-------------loading data and functions-------------

source(here::here("Sea_of_immaturity", "R", "01_a_check_scientific_names.R"))

load(file = here::here("Sea_of_immaturity", "data", "derived_data", "01_b_species_traits_fishbase.Rdata"))
load(file = here::here("Sea_of_immaturity", "data", "derived_data", "01_c_other_traits_raw.Rdata"))

# traits obtained by Raphaël Seguin
load(file = here::here("Sea_of_immaturity", "data", "raw_data", "01_d_merge_species_traits", "fish_traits_RSeguin.Rdata"))
rownames(fish_traits) <- gsub("(?<=\\w)-(?=\\w)", " ", rownames(fish_traits), perl = TRUE)

##-------------Deal with scientific names-------------

species_with_spec_code <- other_traits_raw[!is.na(other_traits_raw$spec_code),]
species_to_find <- other_traits_raw[is.na(other_traits_raw$spec_code),] |> 
  dplyr::select(-c(spec_code, worms_id, fishbase_name))

# /!\ long time to run, and large RAM is needed, reduce mc.core if necessary.
#  put mc_cores = 1 if problem with fishbase
subset1 <- code_sp_check(species_to_find[c(1:5000),], mc_cores = 8)
save(subset1, file = here::here("Sea_of_immaturity", "data", "derived_data", "code_sp_check", "subset1.Rdata"))
subset2 <- code_sp_check(species_to_find[c(5001:10000),], mc_cores = 8) 
save(subset2, file = here::here("Sea_of_immaturity", "data", "derived_data","code_sp_check", "subset2.Rdata"))
subset3 <- code_sp_check(species_to_find[c(10001:nrow(species_to_find)),], mc_cores = 8) 
save(subset3, file = here::here("Sea_of_immaturity", "data", "derived_data","code_sp_check", "subset3.Rdata"))

all_subset <- rbind(subset1, subset2, subset3)
all_subset_delete <- all_subset[all_subset$check!=1,]
escaped <- all_subset_delete |> 
  dplyr::filter(if_any(dplyr::starts_with(c("LengthMatMin_", "LengthMax_")), ~ . != 0))
all_subset <- all_subset[all_subset$check==1,]
summary(all_subset$check) #should be only 1
NA_fishbase <- all_subset[is.na(all_subset$fishbase_name),] #Mostly non fishes -> OK

other_traits_large_matrix <- all_subset |> 
  dplyr::select(-check) |> 
  gtools::smartbind(species_with_spec_code) 

other_traits_large_matrix <- other_traits_large_matrix %>%
  dplyr::mutate(dplyr::across(c(a, b, ForkLength, `VUL-C`, `VUL-F`, PhyloDiv, Troph, LengthMatMin_female, LengthMatMin_male, LengthMatMin_unsexed, LengthMax_unsexed, LengthMax_male, LengthMax_female), as.numeric))

other_traits_large_matrix[] <- lapply(other_traits_large_matrix, function(x) {
  if (is.character(x)) return(factor(x))
  if (is.integer(x)) return(as.numeric(x))
  return(x)
})  


save(other_traits_large_matrix, file = here::here("Sea_of_immaturity", "data", "derived_data", "01_d_other_traits_large_matrix.Rdata") )
load(file = here::here("Documents", "Sea_of_immaturity", "data", "derived_data", "01_d_other_traits_large_matrix.Rdata") )


## Complete data for the same species
other_traits_with_SpecCode <- other_traits_large_matrix |> 
  dplyr::group_by(spec_code) |> 
  tidyr::fill(tidyr::everything(), .direction = 'updown') |> #fill identical lines
  dplyr::mutate(across(.cols = where(is.numeric), .fns = mean, .names = "{.col}")) |> #some rows are became duplicates -> means them.
  dplyr::ungroup() |> 
  dplyr::select( -range_n_cells_005, -worms_id) |> #-ClimVuln_SSP126,
  dplyr::distinct() |> 
  dplyr::filter(!is.na(spec_code))

save(other_traits_with_SpecCode, file = here::here("Sea_of_immaturity", "data", "derived_data", "01_d_other_traits_with_spec_code.Rdata") )
# load( file =here::here("data", "fishbase_data", "V0_other_traits_with_spec_code.Rdata") )

## Check duplicates in names
dup_names <- other_traits_with_SpecCode[duplicated(other_traits_with_SpecCode$spec_code) |
               duplicated(other_traits_with_SpecCode$spec_code, fromLast =T),]
#SOME DUPLICATED SPECIES HAVE SLIGHTLY DIFFERENT VALUES IN TRAITS

other_traits <- other_traits_with_SpecCode |> 
  dplyr::group_by(spec_code) |> 
  # Mean the duplicated numeric traits
  dplyr::mutate(across(.cols = where(is.numeric), .fns = ~mean(., na.rm = TRUE), .names = "{.col}")) |> # Mean numeric variables
 
  #Arrange IUCN categories: keep only the highest IUCN status
  dplyr::mutate(IUCN_inferred_Loiseau23 = 
                  ifelse(IUCN_inferred_Loiseau23 == "No Status", 0,
                         ifelse(IUCN_inferred_Loiseau23 == "Non Threatened", 1,
                                ifelse(IUCN_inferred_Loiseau23 == "Threatened", 2, NA)))) |>
  # dplyr::mutate(IUCN_inferred_Loiseau23 = max(IUCN_inferred_Loiseau23, na.rm = TRUE)) |> 
  dplyr::mutate(IUCN_inferred_Loiseau23 = if(all(is.na(IUCN_inferred_Loiseau23))){
    NA} else max(IUCN_inferred_Loiseau23, na.rm = TRUE)) |> 
  dplyr::mutate(IUCN_inferred_Loiseau23 = 
                  ifelse(IUCN_inferred_Loiseau23 == 0, "No Status",
                         ifelse(IUCN_inferred_Loiseau23 == 1, "Non Threatened",
                                ifelse(IUCN_inferred_Loiseau23 == 2, "Threatened", NA)))) |>
  

    dplyr::distinct(across(-species_name), .keep_all = TRUE)


other_traits[duplicated(other_traits$fishbase_name) | duplicated(other_traits$fishbase_name, fromLast =T),]
#S till problem with diet
other_traits$trophic_guild[other_traits$fishbase_name == "Brachygenys chrysargyreum"] <- "microinvertivore"


## Final data for other traits
other_traits <- other_traits |> 
  dplyr::mutate_all(~ifelse(. == "NaN", NA, .))|> #Replace NaN by NA
  dplyr::select(-species_name) |> 
  dplyr::distinct()

## Taxonomy data

taxonomy <- taxize::classification(species_traits_fishbase$fishbase_name, db="worms")

taxonomy_df <- do.call(rbind, taxonomy)


taxo <- do.call(rbind, taxonomy) |> 
  tibble::rownames_to_column("fishbase_name") |> 
  dplyr::select(-id) |> 
  dplyr::mutate(fishbase_name = gsub("\\.\\d+$", "", fishbase_name)) |> 
  dplyr::filter(rank %in% c("Phylum", "Class", "Order", "Family", "Genus")) |> 
  dplyr::distinct() |> 
  tidyr::pivot_wider(names_from = "rank", values_from = "name", values_fill = NA)

na <- setdiff(species_traits_fishbase$fishbase_name, taxo$fishbase_name)
taxonomy_na <- taxize::classification(na, db="gbif")
taxonomy_na_df <- do.call(rbind, taxonomy_na) |> 
   tibble::rownames_to_column("fishbase_name") |> 
   dplyr::select(-id) |> 
   dplyr::mutate(fishbase_name = gsub("\\.\\d+$", "", fishbase_name)) |> 
   dplyr::mutate(rank = dplyr::recode(rank, phylum = "Phylum", 
                                     cass = "Class", 
                                     order = "Order", 
                                     family = "Family", 
                                     genus = "Genus" )) |>
   dplyr::filter(rank %in% c("Phylum", "Class", "Order", "Family", "Genus")) |> 
   dplyr::distinct() |> 
   tidyr::pivot_wider(names_from = "rank", values_from = "name", values_fill = NA)


taxo <- taxo |>
  dplyr::bind_rows(taxonomy_na_df) 

save(taxo, file = here::here("Documents", "Sea_of_immaturity", "data", "derived_data", "01_d_species_taxonomy.Rdata"))
load(file = here::here("Sea_of_immaturity", "data", "derived_data", "01_d_species_taxonomy.Rdata"))


species_traits_fishbase <- species_traits_fishbase |> 
  dplyr::left_join(taxo, by = "fishbase_name")


##-------------Merge datasets-------------

library(dplyr)

# Supprimer les doublons
species_traits_fishbase <- species_traits_fishbase %>%
  group_by(spec_code) %>%
  arrange(desc(rowSums(!is.na(across(everything()))))) %>% # Trier par le nombre de valeurs non manquantes (descendant)
  slice(1) %>%  # Garder uniquement la première ligne (la plus complète)
  ungroup()

other_traits <- other_traits %>%
  distinct(spec_code, .keep_all = TRUE)

# Gérer les doublons potentiels dans other_traits
other_traits <- other_traits %>%
  group_by(spec_code) %>%
  summarise(across(everything(), ~ paste(unique(.x[!is.na(.x)]), collapse = ", "), .names = "{.col}")) %>%
  ungroup()

# Jointure et fusion des colonnes
merged_traits <- species_traits_fishbase %>%
  left_join(other_traits, by = "spec_code") %>%
  mutate(
    Family = Family.x,
    Troph = coalesce(Troph.x, Troph.y),
    a = dplyr::coalesce(a.x, a.y),
    b = dplyr::coalesce(b.x, b.y),
    LengthMatMin_female = dplyr::coalesce(LengthMatMin_female.x, LengthMatMin_female.y),
    LengthMatMin_male = dplyr::coalesce(LengthMatMin_male.x, LengthMatMin_male.y),
    LengthMatMin_unsexed = dplyr::coalesce(LengthMatMin_unsexed.x, LengthMatMin_unsexed.y),
    VulnerabilityClimate = dplyr::coalesce(VulnerabilityClimate, `VUL-C`)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"), -`VUL-F`, -`VUL-C`, -MLengthRef, -TLObserved) # Nettoyer les colonnes intermédiaires

fish_traits <- fish_traits %>%
  select(Max_length, Common_length, Climate, Repro.Mode, Repro.Fertil) %>%
  mutate(species_name = rownames(fish_traits))

merged_traits_final <- merged_traits %>%
  left_join(fish_traits, by = "species_name") %>%
  mutate(
    LengthMax_unsexed = coalesce(LengthMax_unsexed, Max_length),
    Common_length = coalesce(Length, Common_length),
    Climate = coalesce(Climate.x, Climate.y),
    ReproMode = coalesce(ReproMode, Repro.Mode),
    ReproFertil = Repro.Fertil
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"), -Max_length, -Common_length, -Repro.Mode, -Repro.Fertil)
  
species_traits_final <- merged_traits_final

##-------------Check and observe data-------------
## IUCN data: IUCN categories used to complete inference by Loiseau et al. (2023)
## mainly for Elasmobranch (that are not taken in Loiseau's study)
iucn <- species_traits_final |> 
  tibble::column_to_rownames("species_name") |> 
  dplyr::select(spec_code,
                iucn_inferred = IUCN_inferred_Loiseau23, 
                iucn_redlist = IUCN_category) 

missing_values <- dplyr::filter(iucn, is.na(iucn$iucn_inferred)) |>
  dplyr::mutate(iucn_redlist = dplyr::recode(iucn_redlist,  
                                             "NE" = "No Status",
                                             "DD" = "No Status",
                                             "LC" = "Non Threatened",
                                             "NT" = "Non Threatened",
                                             "VU" = "Threatened",
                                             "EN" = "Threatened",
                                             "CR" = "Threatened")) 

iucn[rownames(missing_values), "iucn_inferred"] <- missing_values$iucn_redlist
iucn$iucn_inferred[iucn$iucn_inferred == "No Status"] <- NA

species_traits_final[, "IUCN_inferred_Loiseau23"] <- iucn[species_traits_final$species_name, 
                                                          "iucn_inferred"]



#Check names
species_traits_final[is.na(species_traits_final$fishbase_name),] #no NA ?
species_traits_final[which(grepl("Genus Species",species_traits_final$fishbase_name)),] #OK
summary(species_traits_final$spec_code) #ok


length(unique(species_traits_final$fishbase_name))# No duplicates
dup_names<- species_traits_final[duplicated(species_traits_final$fishbase_name) |
                            duplicated(species_traits_final$fishbase_name, fromLast =T),]
# FEW NAs IN DUPLICATED NAMES ARE REMAINING IN FISHBASE DATA -> FILL THESE USELESS GAPS:
dup_names <- dup_names |> 
  dplyr::group_by(fishbase_name) |> 
  tidyr::fill(IUCN_category:last_col(), .direction = "downup") |> 
  dplyr::distinct(across(-species_name), .keep_all = TRUE)

#Same on the full dataframe:
species_traits_final <- species_traits_final |> 
  dplyr::group_by(species_name) |> 
  tidyr::fill(IUCN_category:last_col(), .direction = "downup") |> 
  dplyr::select(-trophic_guild)


## Explore data ##
library(funbiogeo)
library(ggplot2)
species_traits <- dplyr::rename(species_traits_final,species = species_name)

fb_plot_species_traits_completeness(species_traits)
ggsave(plot= last_plot(), file= here::here("Documents","Sea_of_immaturity","figures","01_gather_species_info", "traits_completedness.png"), width = 15, height = 7)
fb_plot_number_species_by_trait(species_traits, threshold_species_proportion = 1)
ggsave(plot= last_plot(), file= here::here("Documents","Sea_of_immaturity","figures", "01_gather_species_info", "percent_species_per_traits.png"), width = 8, height = 8)

get_column_classes <- function(data) {
  sapply(data, class)
}

# Get rid of the columns with less than 3% completness
species_traits <- species_traits %>%
  select(-LengthType, -TS, -AirBreathingRef, -MonogamyType, -ParentalCareRef, -SpawnModRef, -Circadian1, -Group, -Level)

# Utiliser la fonction sur ton tableau species_traits
column_classes <- get_column_classes(species_traits)
# Identifier les colonnes character
char_columns <- names(species_traits)[sapply(species_traits, is.character)]
# Vérifier si les colonnes contiennent des valeurs numériques uniquement
potential_numeric <- c("ParentalCareQ", "head_depth_ratio", "lower_jaw_ratio", "body_depth_ratio", "body_width_ratio", "caudalpeduncle_depth_ratio", "caudalpeduncle_width_ratio", "mouth_width_ratio", "depth_min",
                       "depth_max", "ClimVuln_SSP126", "ClimVuln_SSP585", "geographic_range_Albouy19", "range_n_cells_01", "WeightMax", "PhyloDiv", "ForkLength")
# Mutate the columns considered
species_traits <- species_traits %>%
  mutate(across(all_of(potential_numeric), as.numeric))


species_traits_clean <- species_traits |> 
  dplyr::select(where(is.numeric)) |>   #Garder les colonnes numériques
  dplyr::select(where(~ sum(!is.na(.)) > 0))  #Garder les colonnes avec plus d'une valeur non-NA


species_traits_clean <- species_traits_clean |> 
  dplyr::select(-spec_code, -HabitatsRef, -MGillnets, -ParentalCareQ, -FamCode, -GenCode, -LengthFemale, -LengthMatMin_female, -LengthMatMin_male, -LengthMatMin_unsexed, -Emblematic,
                -Abyssopelagic) |> 
  dplyr::mutate(species = species_traits_final$species_name) |>  # Ajout correct de la colonne species
  dplyr::filter(rowMeans(is.na(across(where(is.numeric)))) < 0.1)  # Filtre les lignes avec < 10% de NA dans les colonnes numériques


cor_matrix <- stats::cor(species_traits_clean[,colnames(species_traits_clean)!="species"], use='complete.obs')

fb_plot_trait_correlation(species_traits_clean|>dplyr::select(species, where(is.numeric)))
ggsave(plot= last_plot(), file= here::here("Sea_of_immaturity", "figures","01_gather_species_info", "traits_correlogram.png"), width = 10, height = 10)

correlated_traits_list <- function(cor_matrix, threshold = 0.7) {
  # Initialiser une liste vide
  trait_corr_list <- list()
  
  # Récupérer les noms des traits
  trait_names <- colnames(cor_matrix)
  
  # Parcourir chaque trait
  for (i in seq_len(ncol(cor_matrix))) {
    # Trouver les indices des traits corrélés avec le trait i
    correlated_indices <- which(abs(cor_matrix[i, ]) > threshold & seq_len(ncol(cor_matrix)) != i)
    
    # Ajouter les noms des traits corrélés à la liste, si non vide
    if (length(correlated_indices) > 0) {
      trait_corr_list[[trait_names[i]]] <- setNames(cor_matrix[i, correlated_indices], trait_names[correlated_indices])
    }
  }
  
  return(trait_corr_list)
}

traits_correlated <- correlated_traits_list(cor_matrix, threshold = 0.8)

##-------------Save final data-------------
species_traits_final<-species_traits
save(species_traits_final, file = here::here("Sea_of_immaturity","data", "derived_data", "01_d_species_traits_final.Rdata"))
load( file = here::here("Sea_of_immaturity", "data", "derived_data", "01_d_species_traits_final.Rdata"))
save(species_traits_final, file = here::here("Sea_of_immaturity", "outputs", "01_gather_species_info", "species_traits_before_MF.Rdata"))
