library(readxl)
library(dplyr)

##-------------Process metadata-------------------------

rls_data <- read.csv("~/Sea_of_immaturity/data/raw_data/03_a_process_raw_data/RLS data/2006-2025_Global_reef_fish_abundance_and_biomass.csv", skip=71, header=TRUE) %>%
  mutate(Source = "rls", survey_id = as.character(survey_id)) %>% # Ajouter la colonne Source
  rename(`Sample ID` = survey_id, Binomial = species_name, `Length (mm)` = size_class) %>%
  select(`Sample ID`, country, area, site_code, site_name, latitude, longitude, survey_date, depth, survey_latitude, survey_longitude, ecoregion, realm, location)

bruvs_benthic_metadata <- readxl::read_xlsx("~/Sea_of_immaturity/data/raw_data/03_a_process_raw_data/BRUVS data/MFL Benthic Meta 2024_11_21.xlsx", sheet = 'Benthic Meta') %>%
  mutate(Source = "bruvs_benthic")  # Ajouter la colonne Source

bruvs_benthic_metadata <- bruvs_benthic_metadata %>%
  rename(site_name = Site, `Sample ID`=NewOpcode, survey_date = Date, location= Location, latitude = Latitude, longitude = Longitude, depth= `Depth (m)`, Expedition = `New Exped`)%>%
  select(`Sample ID`, Source, latitude, longitude, site_name, survey_date, location, depth, Expedition, Status, Month, Year) %>%
  mutate(survey_date=as.character(survey_date), depth=as.double(depth), Month = as.character(Month))


bruvs_mesophotic_metadata <- readxl::read_xlsx("~/Sea_of_immaturity/data/raw_data/03_a_process_raw_data/BRUVS data/MFL Mesophotic Database 2025_02_22.xlsx", sheet = 'Meta') %>%
  mutate(Source = "bruvs_pelagic") # Ajouter la colonne Source

bruvs_mesophotic_metadata <- bruvs_mesophotic_metadata %>%
  rename(`Sample ID`= Sample , survey_date = Date, location= Location, latitude = lat, longitude = lon, depth= Depth, Expedition = Survey, realm = Sublocation)%>%
  select(`Sample ID`, Source, latitude, longitude, survey_date, location, depth, Expedition, realm) %>%
  mutate(latitude = as.double(latitude), longitude = as.double(longitude), depth = as.double(depth))

bruvs_midwater_metadata <- readxl::read_xlsx("~/Sea_of_immaturity/data/raw_data/03_a_process_raw_data/BRUVS data/MFL Midwater Database 2025_02_22.xlsx", sheet = 'Meta') %>%
  mutate(Source = "bruvs_pelagic") # Ajouter la colonne Source

bruvs_midwater_metadata <- bruvs_midwater_metadata %>%
  rename( survey_date = Date, location= Location, latitude = `Lat in`, longitude = `Long in`, depth= `Survey Mean Depth`, Expedition = Exped, country = Country)%>%
  select(`Sample ID`, Source, latitude, longitude, survey_date, location, depth, Expedition, country, Month, Year) %>%
  mutate(latitude = as.double(latitude), longitude = as.double(longitude), depth = as.double(depth))

geographic_data <- Measurment_data %>%
  select(`Sample ID`, Expedition, Source) %>%
  distinct()

geographic_data <- geographic_data[geographic_data$Source=="rls",]

geographic_data <- geographic_data %>% 
  left_join(rls_data, by = "Sample ID") %>%
  mutate(Status = NA, Month = NA, Year = NA) %>%
  distinct()

geographic_data <- geographic_data %>%
  bind_rows(bruvs_benthic_metadata)

geographic_data <- geographic_data %>%
  bind_rows(bruvs_mesophotic_metadata) 

geographic_data <- geographic_data %>%
  bind_rows(bruvs_midwater_metadata) %>%
  distinct()

save(Measurment_data, file = here::here("Sea_of_immaturity/outputs/03_process_data/metadata_bruvs_rls.Rdata"))


##-------------------Process raw measurments----------------------------

rls_data <- read.csv("~/Sea_of_immaturity/data/raw_data/03_a_process_raw_data/RLS data/2006-2025_Global_reef_fish_abundance_and_biomass.csv", skip=71, header = TRUE) %>%
  mutate(Source = "rls", survey_id = as.character(survey_id)) %>% # Ajouter la colonne Source
  rename(`Sample ID` = survey_id, Binomial = species_name, Length = size_class)%>%
  select(`Sample ID`, Binomial, Length, Source)

bruvs_benthic_data <- readxl::read_xlsx("~/Sea_of_immaturity/data/raw_data/03_a_process_raw_data/BRUVS data/MFL Benthic FL 2025_02_22.xlsx", sheet = 'FL') %>%
  mutate(Source = "bruvs_benthic", `Length (mm)` = `Length (mm)`/10) %>%
  rename (Length = `Length (mm)`)# Ajouter la colonne Source

bruvs_mesophotic_data <- readxl::read_xlsx("~/Sea_of_immaturity/data/raw_data/03_a_process_raw_data/BRUVS data/MFL Mesophotic Database 2025_02_22.xlsx", sheet = 'Lengths') %>%
  rename(`Sample ID` = OpCode, Length = `Length (mm)`) %>%
  mutate(Source = "bruvs_pelagic", Length = Length/10)  # Ajouter la colonne Source

bruvs_midwater_data <- readxl::read_xlsx("~/Sea_of_immaturity/data/raw_data/BRUVS data/MFL Midwater Database 2025_02_22.xlsx", sheet = 'FL') %>%
  rename(`Sample ID` = `New OpCode`, Expedition = Exped, Length = `Length (mm)`) %>%
  mutate(Source = "bruvs_pelagic", Length = Length/10) # Ajouter la colonne Source

conversion_TLFL <- readxl::read_xlsx("~/Sea_of_immaturity/data/raw_data/03_a_process_raw_data/BRUVS data/BRUVS taxa list 2025_01_6.xlsx", sheet = 'Taxa') %>%
  rename(Binomial = Taxa)

# Sélectionner les colonnes communes
common_cols <- intersect(colnames(bruvs_benthic_data), colnames(bruvs_mesophotic_data))
common_cols <- intersect(common_cols, colnames(bruvs_midwater_data))

bruvs_benthic_data <- bruvs_benthic_data %>% select(all_of(common_cols))
bruvs_mesophotic_data <- bruvs_mesophotic_data %>% select(all_of(common_cols))
bruvs_midwater_data <- bruvs_midwater_data %>% select(all_of(common_cols))

# Fusionner les deux bases avec la colonne Source
bruvs_data <- rbind(bruvs_benthic_data, bruvs_mesophotic_data)
bruvs_data <- rbind(bruvs_data, bruvs_midwater_data)

bruvs_data <- bruvs_data %>%
  left_join(conversion_TLFL %>% select(Binomial, TLFL), by= "Binomial") %>%
  mutate(TLFL = as.numeric(TLFL), Length = ifelse(!is.na(TLFL), Length*TLFL, Length))

Measurment_data <- bruvs_data %>%
  bind_rows(rls_data)

save(Measurment_data, file = here::here("Sea_of_immaturity/outputs/03_process_data/data_bruvs_rls.Rdata"))

