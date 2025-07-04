load("~/Sea_of_immaturity/outputs/04_get_vars/protection/in_out_mpa.Rdata")
mpa <- in_out_mpa
load("~/Sea_of_immaturity/outputs/03_process_data/metadata_bruvs_rls.Rdata")

library(lubridate)

mpa <- mpa %>%
  mutate(protection = case_when(
    removal_of_marine_life_is_prohibited %in% c(5) ~ "no.take",
    removal_of_marine_life_is_prohibited %in% c(4, 3) ~ "restricted",
    removal_of_marine_life_is_prohibited %in% c(2) ~ "less.restricted",
    removal_of_marine_life_is_prohibited == 1 ~ "least.restricted",
    TRUE ~ "fishing"  # Au cas où il y aurait des valeurs autres
  ))
mpa[mpa$iucn_cat %in% c("Not Reported","Unassigned"),]$iucn_cat <- NA

mpa <- mpa %>%
  rename(survey_id = `Sample ID`)

# on load proprement geographic_data
geographic_data <- geographic_data |> 
  rename(survey_id=`Sample ID`, site_code_x= site_code) %>%
  mutate(survey_date_ymd = ymd(geographic_data$survey_date))

geographic_data$survey_date <- lubridate::as_date(geographic_data$survey_date) 

library(sf)
library(dplyr)
library(tidyr)

# Charger les coordonnées des sites d’échantillonnage
sites <- geographic_data %>% dplyr::select(survey_id, latitude, longitude, survey_date, mpa_effectiveness) 

# filter to keep valid coords
sites <- sites %>%
  filter(latitude >= -90 & latitude <= 90, longitude >= -180 & longitude <= 180) %>%
  drop_na(latitude, longitude)  # Supprimer les lignes avec des NA

# On ajoute les dates des survey et l'effectiveness
mpa <- mpa %>% 
  left_join(sites %>% select(survey_id, survey_date, mpa_effectiveness), by= "survey_id") %>%
  rename(size = total_area, year_of_protection = year_est , Effectiveness = mpa_effectiveness) %>%
  mutate(year = as.numeric(lubridate::year(survey_date)))

mpa$year_of_protection <- as.numeric(mpa$year_of_protection)

mpa <- mpa |>
  dplyr::mutate(Effectiveness = ifelse(year_of_protection < year |
                                         is.na(year_of_protection),
                                       Effectiveness, NA),
                category_name= ifelse(year_of_protection < year |
                                         is.na(year_of_protection),
                                       category_name, NA),
                size = ifelse(year_of_protection < year |
                                is.na(year_of_protection),
                              size, NA),
                age = ifelse(year - year_of_protection >= 0,
                             year - year_of_protection, NA))

# assigned protection status
mpa_assigned <- mpa |>
    dplyr::mutate(
      protection_status = dplyr::case_when(
        # Hors AMP
        inside_mpa == FALSE ~ "out",
        
        # Haute protection
        #Effectiveness %in% c("High") |
          (protection == "no.take" &
          #category_name %in% c("IUCN MPA", "Jurisdictional Authority Area", "Fisheries Management Area", "OECM", "Water Quality/Human Health Area", NA) & 
          year > year_of_protection) ~ "high",
        
        # Moyenne protection
        #Effectiveness %in% c("Medium") |
          (protection == "restricted" &
          #category_name %in% c("IUCN MPA", "Jurisdictional Authority Area", "Fisheries Management Area", "OECM", "Water Quality/Human Health Area", NA) &
          (year > year_of_protection)) ~ "medium",
        
        # Faible protection
        (protection %in% c("less.restricted") & (year > year_of_protection) )  ~ "low",

        # Défaut
        TRUE ~ "out" )
    )
    

table(mpa_assigned$protection_status, useNA = "always")
table(mpa_assigned$protection_status2, useNA = "always")


sites_w_mpa_status <- mpa_assigned %>%
  mutate(
    inside_mpa = case_when(
      category_name %in% c("IUCN MPA") ~ TRUE,
      !(category_name %in% c("IUCN MPA")) ~ FALSE,
      TRUE ~ NA  # pour les cas où category_name est NA
    )
  )


sites_w_mpa_status <- sites_w_mpa_status %>%
  select(survey_id, protection_status, protection, category_name, inside_mpa, removal_of_marine_life_is_prohibited, iucn_cat)

save(sites_w_mpa_status, file=here::here("Sea_of_immaturity", "outputs", "04_get_vars", "protection", "mpa_mondial_protection_status.Rdata"))
