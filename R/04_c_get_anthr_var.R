
## load geographic_data
load("~/Sea_of_immaturity/outputs/03_process_data/metadata_bruvs_rls.Rdata")

##------------mpa protection status Graham-----------------------------------

mpa <- read_csv("Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/mpa_rls.csv")

mpa <- mpa %>%
  rename(`Sample ID`= survey_id) %>%
  mutate(`Sample ID`= as.character(`Sample ID`))

geographic_data <- geographic_data %>%
  left_join(mpa, by=c("Sample ID"))

geographic_data$Status[!(geographic_data$Status %in% c("MPA", "FISHED", "Reserve", "No-take"))] <- NA


##------------sites to consider--------------------------------------------
library(sf)
library(dplyr)
library(tidyr)

# Charger les coordonnées des sites d’échantillonnage
sites <- geographic_data %>% dplyr::select(`Sample ID`, latitude, longitude) %>% na.omit(latitude, longitude)

# filter to keep valid coords
sites <- sites %>%
  filter(latitude >= -90 & latitude <= 90, longitude >= -180 & longitude <= 180) %>%
  drop_na(latitude, longitude)  # Supprimer les lignes avec des NA

#sites$survey_date <- sub(" .*", "", sites$survey_date)

# Convertir `sites` en sf 
sites_sf <- st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326)
##----------distance to port-----------------------------------------

library(geosphere)
library(dplyr)
library(readr)
library(tidyr)
# Charger les données des ports
ports_df <- read.csv("~/Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/Positions_ports.csv")

# Garder uniquement les colonnes utiles (latitude et longitude des ports)
ports_coords <- ports_df %>% select(Latitude, Longitude, Main.Port.Name) %>% drop_na()

# Fonction pour calculer la distance minimale entre un site et les ports
min_distance_to_port <- function(lat, lon, ports_coords) {
  site_coords <- cbind(lon, lat)  
  port_coords <- cbind(ports_coords$Longitude, ports_coords$Latitude)  # Extraction propre des coordonnées
  
  # Calculer les distances en km
  distances <- distGeo(site_coords, port_coords) / 1000
  
  # Trouver l'index du port le plus proche
  min_index <- which.min(distances)
  
  # Extraire la distance et le nom du port correspondant
  min_distance <- distances[min_index]
  proximate_port <- ports_coords$Main.Port.Name[min_index]  # Assure-toi que cette colonne existe dans ports_coords
  
  return(list(distance_to_port = min_distance, proximate_port = proximate_port))
}

# Ajouter la colonne "distance_to_port"
sites <- sites %>%
  rowwise() %>%
  mutate(result = list(min_distance_to_port(latitude, longitude, ports_coords)),
         distance_to_port = result$distance_to_port,
         proximate_port = result$proximate_port) %>%
  select(-result) %>%  # Supprime la colonne intermédiaire
  ungroup()

distance_to_port <- sites %>%
  select(`Sample ID`, distance_to_port, proximate_port)

geographic_data <- geographic_data %>%
  left_join(distance_to_port, by=c("Sample ID"))

# Sauvegarder le fichier avec la nouvelle colonne
#write_csv(geographic_data, "~/chemin/vers/geographic_data_with_ports.csv")

# Afficher un aperçu des résultats
print(head(geographic_data))

##---------------gravity---------------------------------------

load("~/Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/gravity/dataGravity_201223.RData")

library(sf)

# Charger le raster de gravité
gravity_raster <- rast("~/Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/gravity/rastGravity.tif") 

dataGravity <- dataGravity %>%
  rename(longitude = lon, latitude = lat) %>%
  select(longitude, latitude, TotalGravity) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 3857) %>%  # CRS des données
  st_transform(crs = 4326) %>%  # Conversion en WGS84
  mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2]) %>%
  st_drop_geometry()  # Retirer la colonne géométrique pour revenir à un dataframe normal

# Convertir `dataGravity` en sf si ce n'est pas déjà fait
dataGravity_sf <- st_as_sf(dataGravity, coords = c("longitude", "latitude"), crs = 4326)

# Effectuer un joint spatial pour récupérer les valeurs de gravité les plus proches
sites_with_gravity <- st_join(sites_sf, dataGravity_sf, join = st_nearest_feature)

# Ajouter la valeur de gravité extraite à `sites`
sites$Gravity <- sites_with_gravity$TotalGravity

# Vérifier les NA
summary(sites$Gravity)

##------------Global protection status---------------------------------------

library(sf)
library(nngeo)  # Pour les distances aux polygones
library(lwgeom)# Pour les tests de validité
library(pbmcapply)
library(terra)
library(readxl)

# Charger les fichiers
mondial_mpas_metadata <- read.csv("Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/ProtectedSeas_Navigator_data_20241212.csv") %>% rename( SITE_NAME = site_name)
mondial_mpas <- st_read("Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/ProtectedSeas_Navigator_20241212_shp/ProtectedSeas_Navigator_20241212_shp.shp")

# Vérifier si tous les polygones sont valides
invalid_mpas <- mondial_mpas[!st_is_valid(mondial_mpas), ]

# Voir combien de polygones sont invalides
nrow(invalid_mpas)

mondial_mpas <- st_make_valid(mondial_mpas)

# Sont-ils tous valides ?
nrow(invalid_mpas)

# Keep valid polygons
mondial_mpas <- mondial_mpas[st_is_valid(mondial_mpas), ]

mondial_mpas$area <- as.numeric(st_area(mondial_mpas))

# Add fishing info
mondial_mpas <- merge(mondial_mpas, mondial_mpas_metadata[, c("SITE_NAME", "removal_of_marine_life_is_prohibited", "category_name")], by = "SITE_NAME", all.x = TRUE)

# Vérifier si les points sont dans une MPA
intersections <- st_intersects(sites_sf, mondial_mpas)

# Adapter ici les noms de colonnes dans le shapefile si différent
mpa_name_col <- "SITE_NAME"   # Remplace par le vrai nom de colonne
mpa_code_col <- "SITE_ID" # Remplace par le vrai nom de colonne

sites_sf$inside_mpa <- lengths(intersections) > 0

# test <- pbmcapply::pbmclapply(1:nrow(sites_sf), function(i) {
#   
#   if (sites_sf$inside_mpa[i] && length(intersections[[i]]) > 0) {
#     mpa_index = intersections[[i]][1]
#     # Prendre la première MPA dans laquelle le point tombe
#     mpa_dtf <- data.frame(
#                           mpa_name = mondial_mpas[[mpa_name_col]][mpa_index],
#                           mpa_code = mondial_mpas[[mpa_code_col]][mpa_index],
#                           mpa_dist = 0) # À l’intérieur
#   } else {
#     # Trouver la MPA la plus proche
#     nearest_index <- st_nearest_feature(sites_sf[i, ], mondial_mpas)
#     mpa_dtf <- data.frame(
#                           mpa_name = mondial_mpas[[mpa_name_col]][nearest_index],
#                           mpa_code = mondial_mpas[[mpa_code_col]][nearest_index],
#                           mpa_dist = st_distance(sites_sf[i, ], mondial_mpas[nearest_index, ]))
#   }
#   
# }, mc.cores = 6)
# test_bind <- do.call(rbind, test)


# # Charger les données en format vecteur terra
# mondial_mpas_vect <- terra::vect(mondial_mpas)
# sites_vect <- terra::vect(sites_sf)  # si tes points étaient un objet sf ; sinon adapte selon ton input
# 
# # Assurer la même projection
# mondial_mpas_vect <- project(mondial_mpas, "EPSG:4326")
# sites_vect <- project(sites_vect, mondial_mpas)
# 
# # Vérifier l'appartenance à une MPA
# intersections_vect <- relate(sites_vect, mondial_mpas_vect, "within")  # liste d’indices

# Adapter ici les noms de colonnes dans le shapefile si différent
mpa_name_col <- "SITE_NAME"
mpa_code_col <- "SITE_ID"

# Créer un vecteur logique inside_mpa
inside_mpa <- lengths(intersections) > 0

# Boucle parallèle
test <- pbmclapply(1:nrow(sites_sf), function(i) {
  
  if (inside_mpa[i] && length(intersections[[i]]) > 0) {
    ix <- intersections[[i]]
    mpa_index <- ix[which.max(mondial_mpas$removal_of_marine_life_is_prohibited[ix])]
    data.frame(
      mpa_name = mondial_mpas[[mpa_name_col]][mpa_index],
      mpa_code = mondial_mpas[[mpa_code_col]][mpa_index],
      mpa_dist = "in"
    )
  } else {
    # Trouver la MPA la plus proche
    #nearest_index <- which.min(distance(sites_vect[i, ], mondial_mpas_vect))
    data.frame(
      mpa_name = NA,
      mpa_code = NA,
      mpa_dist = "out"
    )
  }
  
}, mc.cores = 2)


# # Priority by category
# priority_order <- c("IUCN MPA", "OECM", "Jurisdictional Authority Area", "Fisheries Management Area", "Voluntary Conservation Measure Area",
#                     "Vessel Restricted Area", "Water Quality/Human Health Area", "Vessel Reporting Area", "Recreational Area",
#                     "Other Recreational Area", "Other")
# 
# # Fonction pour obtenir le rang de priorité avec gestion des valeurs inconnues
# get_priority <- function(cat) {
#   p <- match(cat, priority_order)
#   ifelse(is.na(p), Inf, p)  # Priorité très basse pour les non listées ou NA
# }
# 
# # Appliquer la logique avec la priorité des catégories
# # Application principale
# test <- pbmclapply(1:nrow(sites_sf), function(i) {
#   
#   if (inside_mpa[i] && length(intersections[[i]]) > 0) {
#     ix <- intersections[[i]]
#     
#     categories <- mondial_mpas$category_name[ix]
#     priorities <- sapply(categories, get_priority)
# 
#     #Si IUCN est présent
#     if ("IUCN MPA" %in% categories) {
#       iucn_ix <- ix[categories == "IUCN MPA"]
#       mpa_index <- iucn_ix[which.max(mondial_mpas$removal_of_marine_life_is_prohibited[iucn_ix])]
#     } else {
#       # Autres : combinaison priorité + max removal_of_marine_life_is_prohibited
#       best_priority <- min(priorities)
#       best_ixs <- ix[priorities == best_priority]
#       removal_vals <- mondial_mpas$removal_of_marine_life_is_prohibited[best_ixs]
#       mpa_index <- best_ixs[which.max(removal_vals)]
#     }
#     
#     data.frame(
#       mpa_name = mondial_mpas[[mpa_name_col]][mpa_index],
#       mpa_code = mondial_mpas[[mpa_code_col]][mpa_index],
#       mpa_category = mondial_mpas$category_name[mpa_index],
#       mpa_dist = "in"
#     )
#   } else {
#     data.frame(
#       mpa_name = NA,
#       mpa_code = NA,
#       mpa_category = NA,
#       mpa_dist = "out"
#     )
#   }
#   
# }, mc.cores = 2)

# Coller les résultats
test_bind <- do.call(rbind, test)

# Ajouter les résultats à l’objet `sites_vect`
sites_sf$inside_mpa <- inside_mpa
sites_sf$mpa_name <- test_bind$mpa_name
sites_sf$mpa_code <- test_bind$mpa_code
sites_sf$mpa_dist <- test_bind$mpa_dist

## import csv 
library(readr)
mondial_mpas_metadata <- read_csv("Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/ProtectedSeas_Navigator_data_20241212.csv") 
mondial_mpas_metadata <- mondial_mpas_metadata %>%
  rename(site_name=SITE_NAME)

in_out_mpa <- sites_sf %>%
  rename(site_name=mpa_name, site_id=mpa_code)

in_out_mpa <- in_out_mpa %>%
  left_join(mondial_mpas_metadata, by = c("site_id", "site_name"))

##-------------------GDP----------------------------------------------

library(readr)
GDP_mondial <- read_csv("Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/0_25deg/final_GDP_0_25deg_postadjust_pop_dens_0_adjust.csv") 
GDP_mondial <- GDP_mondial %>%
  select(latitude, longitude, predicted_GCP_const_2017_PPP, cell_id, year, iso) %>%
  rename(gdp=predicted_GCP_const_2017_PPP)
GDP_mondial <- GDP_mondial %>%
  group_by(longitude, latitude) %>%
  slice_max(order_by = year, n = 1, with_ties = FALSE) %>%
  ungroup() # Keep the most recent value

gdp_vect <- terra::vect(GDP_mondial, geom = c("longitude", "latitude"), crs = "WGS84", keepgeom = TRUE)# Creat a RLS spatial object
sites_vect <- terra::vect(sites, geom = c("longitude", "latitude"), crs = "WGS84", keepgeom = TRUE)# Assess the gdp closest value from each RLS pointsclosest_gdp <- pbmcapply::pbmclapply(1:nrow(rls_vect), function(i) {  nn_cell <- terra::as.data.frame(terra::nearest(rls_vect[i], gdp_vect))  if(nn_cell$distance < 100000) { # if the distance is greater to 100km (100000m), gdp value = NA, otherwise associate the gdp value to the closest point    data_value <- as.numeric(data.frame(gdp_vect[nn_cell$to_id,])[,"gdp"])


# Assess the gdp closest value from each RLS 
pointsclosest_gdp <- pbmcapply::pbmclapply(1:nrow(sites_vect), function(i) {  
  nn_cell <- terra::as.data.frame(terra::nearest(sites_vect[i], gdp_vect))  
  if(nn_cell$distance < 100000) { # if the distance is greater to 100km (100000m), gdp value = NA, otherwise associate the gdp value to the closest point    
    data_value <- as.numeric(data.frame(gdp_vect[nn_cell$to_id,])[,"gdp"])
    data_n <- nn_cell$distance
    data_objects <- c(data_value, data_n, terra::geom(sites_vect[i])[3:4], nn_cell$to_id, nn_cell$to_x, nn_cell$to_y)  
  }
  else{    
    data_value <- NA
    data_n <- nn_cell$distance
    data_objects <- c(data_value, data_n, terra::geom(sites_vect[i])[3:4], nn_cell$to_id, nn_cell$to_x, nn_cell$to_y)  
    }
  }, mc.cores = parallel::detectCores() - 1)

closest_gdp_bind <- do.call(rbind, pointsclosest_gdp) |>
  dplyr::as_tibble()

colnames(closest_gdp_bind) <- c("gdp", "distance", "longitude", "latitude", "n_cell", "x", "y")
closest_gdp_bind <- closest_gdp_bind[-which(is.na(closest_gdp_bind$gdp)),] # only 6 NA
closest_gdp_bind <- closest_gdp_bind[,c("gdp", "longitude", "latitude")]
gdp_final <- sites |>
  dplyr::inner_join(closest_gdp_bind, multiple = "first")
gdp <- gdp_final[,c("Sample ID", "longitude", "latitude", "gdp")]
save(gdp, file = "Sea_of_immaturity/outputs/04_get_vars/human_pressures/sites_with_gdp.RData")

##----------------Ecoregion-------------------------------------

mondial_ecoregion <- st_read("Sea_of_immaturity/data/raw_data/04_c_get_anthro_var/ecoregion_shp/meow_ecos.shp")

invalid_ecoregion <- mondial_ecoregion[!st_is_valid(mondial_ecoregion), ]

# Voir combien de polygones sont invalides
nrow(invalid_ecoregion)

# Vérifier si les points sont dans une MPA
intersections <- st_intersects(sites_sf, mondial_ecoregion)

# Adapter ici les noms de colonnes dans le shapefile si différent
Province <- "PROVINCE"   # Remplace par le vrai nom de colonne
Realm <- "REALM" # Remplace par le vrai nom de colonne

sites_sf$inside_ecoregion <- lengths(intersections) > 0

classify <- pbmcapply::pbmclapply(1:nrow(sites_sf), function(i) {
  
  if (sites_sf$inside_ecoregion[i] && length(intersections[[i]]) > 0) {
    ecoregion_index = intersections[[i]]
    # Prendre la première MPA dans laquelle le point tombe
    ecoregion_dtf <- data.frame(
      `Sample ID` = sites_sf$`Sample ID`[i],
      province = mondial_ecoregion[[Province]][ecoregion_index],
      realm = mondial_ecoregion[[Realm]][ecoregion_index]) # À l’intérieur
  } else {
    # Trouver la MPA la plus proche
    ecoregion_dtf <- data.frame(
      `Sample ID` = sites_sf$`Sample ID`[i],
      province = NA,
      realm = NA)
  }
  
}, mc.cores = 6)
classify_bind <- do.call(rbind, classify)

t$Province <- classify_bind$province
t$Realm <- classify_bind$realm

sites_sf <- sites_sf %>%
  left_join(classify_bind  %>% rename(`Sample ID`= Sample.ID) , by ="Sample ID")
