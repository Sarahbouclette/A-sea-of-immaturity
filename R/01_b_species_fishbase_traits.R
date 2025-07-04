################################################################################
##
##  Takes species list with fishbase name and make a Species x Traits matrix with
##   relevant traits from fishbase, IUCN category, and taxonomy.
##
## 1a_species_fishbase_trait.R
##
## 21/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "rredlist", "remotes")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


##-------------loading data-------------
load(file = here::here("Sea_of_immaturity","data", "raw_data", "01_b_species_fishbase_traits", "data_rfishbase.RData")) #/!\ large file

# Load the RData file and capture the name of the objects loaded
loaded_objects <- load(here::here("Sea_of_immaturity", "data", "derived_data", "00_list_species.RData"))

# Assuming the object in the RData file is named 'species_list', assign it to 'list_sp'
list_sp <- get(loaded_objects[1])


redlist_key <- "deaf690b2d1f6c50909040142e3e0fc540c1401f981996eaa222bd4ac97760f1"

source(here::here("Sea_of_immaturity","R", "01_a_check_scientific_names.R"))


##-------------Extract accepted names for the species list----------------###
colnames(list_sp) <- "species_name"
# Check names and synonyms (use the FISHBASE database)
species_list <- code_sp_check(list_sp, mc_cores = 1) # /!\ long time to run #if problem, try with only one core.*-


##-------------Check the data-------------
length(unique(species_list$spec_code)) # obs duplicated species in the data base 
duplications <- species_list[which(duplicated(species_list$spec_code)),] 
unique(species_list$check) #should be only 1
species_list[is.na(species_list$fishbase_name),] # No NA ?
summary(species_list$spec_code) 
head(species_list$spec_code[order(species_list$spec_code)]) 
species_list[which(grepl("Genus Species",species_list$fishbase_name)),]
#Check differences of naming 
diff <- species_list[
  which(species_list$species_name != species_list$fishbase_name),] #113 different names -> OK

species_list <- species_list %>%
  dplyr::filter(species_name != "name_checklist") %>%
  group_by(species_name) %>%
  filter(all(check == 1)) %>%
  ungroup()

##-------------Complete and save species list-------------
species_list <- dplyr::select(species_list, -check)

### Save Species list ###
save(species_list, file = here::here("Sea_of_immaturity","data", "derived_data", "01_a_species_list.Rdata"))
# load( file = here::here("data", "derived_data", "00_species_list.Rdata"))

##-------------Explore fishbase data / Extract traits-------------
variables <- colnames(data_rfishbase)
variables[order(variables)]
traits_list <- c("SpecCode", "Species", "Length", "Importance",
                 "PriceCateg", "UsedforAquaculture", "Aquarium", 
                 "Troph", "a",  "b", "K", "TempPrefMean", "Schooling", "VulnerabilityClimate", "Abyssopelagic", 
"Area", "Emblematic", "TLObserved", "Circadian1", "Status","CommonLength", "Weight", "Length", "LengthFemale", "Status", "Climate", "AverageDepth", "TempPrefMin", "AgeMin", "HabitatsRef",
"Fertilization", "Spawning", "ParentalCareQ", "NorthernLat", "SouthernLat", "WesternLat", "EasternLat", "n_T", "MatingSystem", "ParentalCareRef", "MGillnets", "Salinity", "MonogamyType", "SpawnModRef", "TS", "Troph", "LengthType",
"FoodTroph", "AirBreathingRef", "MLengthRef", "TempPrefMin", "ReproMode", "ParentalCare", "TempPrefM", "FamCode", "GenCode")
traits_list <- unique(traits_list)


fb_traits <- data_rfishbase |> 
  dplyr::select(dplyr::all_of(traits_list)) |> 
  dplyr::filter(SpecCode %in% species_list$spec_code) |> 
  dplyr::distinct() |> 
  dplyr::mutate(Schooling = factor(Schooling)) |>
  dplyr::mutate(Schooling= forcats::fct_recode(Schooling,yes="-1",no="0") )
  #dplyr::rename(spec_code = SpecCode)

#Check extraction
setdiff(fb_traits$SpecCode, species_list$spec_code) #OK
setdiff(species_list$spec_code, fb_traits$SpecCode) #OK
# rm(list=c("data_rfishbase"))

#Obs traits
length(unique(fb_traits$Species)) # one row per species ?
duplications <- fb_traits[which(duplicated(fb_traits$SpecCode) |
                                  duplicated(fb_traits$SpecCode, fromLast = T)),]
mistake <- which(fb_traits$SpecCode == 24 & fb_traits$Schooling == "no")
if(length(mistake) > 0) fb_traits <- fb_traits[-mistake, ]

unique(fb_traits$Importance)
fb_traits$Importance[which(fb_traits$Importance==" ")] <- NA
summary(fb_traits$Length)
unique(fb_traits$PriceCateg)
fb_traits$PriceCateg[which(fb_traits$PriceCateg=="unknown")] <- NA
unique(fb_traits$UsedforAquaculture)
unique(fb_traits$Aquarium)
table(fb_traits$Schooling)
summary(fb_traits$a)
summary(fb_traits$b)
summary(fb_traits$K)

fb_traits <- dplyr::rename(fb_traits, spec_code = SpecCode) |> 
  dplyr::select(- Species)

fb_traits <- fb_traits %>%
  dplyr::distinct(spec_code, .keep_all= TRUE)

# One row missing
new_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(fb_traits)))
colnames(new_row) <- colnames(fb_traits)
new_row$spec_code <- 70468
new_row$Length <- 11.6
fb_traits <- rbind(fb_traits, new_row)

##-------------Code water position + vulnerability to fishing-------------
species_traits_fb <- rfishbase::species() |> 
  dplyr::filter(SpecCode %in% species_list$spec_code)

colnames(species_traits_fb)[order(colnames(species_traits_fb))]

water_pos_vuln <- species_traits_fb |> 
  dplyr::select(SpecCode, DemersPelag, Vulnerability)
table(water_pos_vuln$DemersPelag)
water_pos_vuln$DemersPelag[which(water_pos_vuln$DemersPelag=="unknown")] <- NA
summary(water_pos_vuln$Vulnerability)

##-------------LengthMatMin-----------------------------------------------

# Load the RData file with the LengthMatMin data
loaded_objects <- load(here::here("Sea_of_immaturity", "data", "raw_data","01_b_species_fishbase_traits", "LengthMatMin.Rdata"))

# Assuming the object in the RData file is named 'species_list', assign it to 'list_sp'
LengthMatMin <- get(loaded_objects[1])
LengthMatMin[LengthMatMin$Sex%in%c('unsexed', 'Unsexed', 'mixed'),]$Sex <-'unsexed'
LengthMatMin[LengthMatMin$Sex%in%c('female', 'Female'),]$Sex <-'female'
LengthMatMin[LengthMatMin$Sex%in%c('male', 'Male'),]$Sex <-'male'

# Transformation des données pour créer une colonne pour chaque sexe
LengthMatMin <- LengthMatMin %>%
  # Pivot pour transformer les sexes en colonnes distinctes
  pivot_wider(
    names_from = Sex,
    values_from = LengthMatMin,
    names_prefix = "LengthMatMin_"
  ) %>%
  # Identifier les colonnes qui contiennent des listes et les traiter
  mutate(across(starts_with("LengthMatMin_"), ~ {
    if (is.list(.)) {
      # Si la colonne est une list-col, calculer la moyenne pour chaque cellule
      sapply(., function(x) mean(as.numeric(x), na.rm = TRUE))
    } else if (grepl(",", as.character(.), fixed = TRUE)) {
      # Si les valeurs sont des chaînes de caractères séparées par des virgules, les convertir en liste
      mean(as.numeric(strsplit(as.character(.), ",")[[1]]), na.rm = TRUE)
    } else {
      .  # Sinon, garder la valeur d'origine
    }
  })) %>%
  # Identifier et résoudre les doublons dans les colonnes pivotées
  group_by(SpecCode, Species) %>%
  summarise(across(starts_with("LengthMatMin_"), ~ mean(., na.rm = TRUE)), .groups = "drop")


LengthMatMin <- LengthMatMin %>%
  dplyr::mutate(across(where(is.numeric), ~ na_if(., NaN))) %>%
  dplyr::rename(spec_code = SpecCode) %>%
  dplyr::distinct(spec_code, .keep_all = TRUE) %>%
  dplyr::select(-Species)

#save(LengthMatMin, file=here::here("Documents", "Sea_of_immaturity", "data", "raw_data", "LengthMatMin_fishbase.Rdata"))

##-------------IUCN status-------------
#Here we extract all IUCN data from the redlist database. Some missing values are
# present. See '1d_bathymetry_and_vulnerability_climate.R' for the inferred data
# from Loiseau et al. 2023.

uicn <- rredlist::rl_sp(all=T, key=redlist_key)
all_iucn <- do.call(rbind, lapply(uicn, "[[", "result"))

#deal with names
iucn_rls_code <- dplyr::select(species_list,
                               species_name, fishbase_name, spec_code) |>
  dplyr::mutate(IUCN_category = NA)

for( i in c(1:nrow(iucn_rls_code))){
  sp_name <- iucn_rls_code$species_name[i]
  fb_name <- iucn_rls_code$fishbase_name[i]
  
  iucn_sp <- all_iucn$category[which(all_iucn$scientific_name == sp_name)]
  iucn_fb <- all_iucn$category[which(all_iucn$scientific_name == fb_name)]
  
  if(length(iucn_fb) == 1){ iucn_rls_code$IUCN_category[i] <- iucn_fb }
  if(length(iucn_sp) == 1){ iucn_rls_code$IUCN_category[i] <- iucn_sp }
  
  if(length(iucn_fb) == 1 & length(iucn_sp) == 1){
    if(iucn_fb != iucn_sp){ cat("Disagreement between", sp_name, ": ", iucn_sp,
                                " and ", fb_name, ": ", iucn_fb,
                                "\n", sp_name, ": ", iucn_sp," is took" )}}
}

length(unique(iucn_rls_code$species_name)) #no duplications
table(iucn_rls_code$IUCN_category)
sum(is.na(iucn_rls_code$IUCN_category)) #471 NAs

# Deplete duplications in the col species_code
iucn_rls_code <- iucn_rls_code |> 
  dplyr::distinct(spec_code, .keep_all = TRUE)


##-------------Merge fishbase data, water position, IUCN data-------------
water_pos_vuln <- dplyr::rename(water_pos_vuln, spec_code = SpecCode)

species_traits_fishbase <- species_list |> 
  dplyr::left_join(iucn_rls_code, by = "spec_code") |> 
  dplyr::left_join(fb_traits, by = "spec_code") |> 
  dplyr::left_join(water_pos_vuln, by = "spec_code") |>
  dplyr::left_join(LengthMatMin, by = "spec_code")

##-------------Save data-------------
species_traits_fishbase <- species_traits_fishbase %>%
  dplyr::select(-species_name.x, -fishbase_name.x)
  dplyr::rename(species_name = species_name.y, fishbase_name = fishbase_name.y, Vulnerability_fishing = Vulnerability)


save(species_traits_fishbase, file = here::here("Sea_of_immaturity","data", "derived_data", "01_b_species_traits_fishbase.Rdata"))
# load( file = here::here("data", "derived_data", "1a_species_traits_fishbase.Rdata"))