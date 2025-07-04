# Install and load necessary packages
install.packages("reticulate")
library(reticulate)
library(dplyr)

# Install the geoenrich Python package
py_install("geoenrich", pip = TRUE)

# Import Python modules
os <- import("os")
dataloader <- import("geoenrich.dataloader", convert = FALSE)
enrichment <- import("geoenrich.enrichment", convert = FALSE)
exports <- import("geoenrich.exports")

# Charger les données d'occurrence 
sites_data <- dataloader$import_occurrences_csv(path = 'Sea_of_immaturity/outputs/04_gather_vars/sites_data.csv', id_col = 'survey_id', lat_col='latitude', lon_col='longitude', date_col= 'survey_date') 

# Créer un fichier d'enrichissement pour stocker les métadonnées
dataset_ref <- 'var_env_7days'
enrichment$create_enrichment_file(sites_data, dataset_ref)

# Définir les variables environnementales à récupérer
variables <- c('sst',  'chlorophyll')  # Température et Chlorophylle-a
geo_buff <- 0  # Buffer spatial en km
time_buff <- c(-7L, 0L) # Moyenne sur 7 jours avant la date de relevé

#Enrichir les données avec les variables téléchargées

# Définir les slices de 10 000, et finir avec 20000-30000
slices <- list(
  c(0L, 10000L),
  c(10001L, 20000L),
  c(30001L, 40000L),
  c(40001L, 50000L),
  c(50001L, 52111L),
  c(20001L, 25500L)
)

# Appliquer la fonction enrich pour chaque slice
for (var_id in variables) {
  for (slice in slices) {
    enrichment$enrich(dataset_ref, "distance-to-coast", geo_buff, time_buff, slice = slice)
  }
}


for (var_id in variables){
  enrichment$enrich(dataset_ref, var_id, geo_buff, time_buff, slice = c(26100L, 27000L))
}

# Check enrichment status
enrichment$enrichment_status(dataset_ref)

# Export enriched data as raster and PNG for visualization

ids <- py_to_r(enrichment$read_ids(dataset_ref))

output 

for (occ_id in ids) {
  output <- exports$retrieve_data(dataset_ref, occ_id, var_id, shape = 'buffer')
}

exports$export_raster(dataset_ref, occ_id, var_id, path = './')

for (var_id in variables) {
  exports$produce_stats(dataset_ref, 'distance-to-coast')
}

library(ncdf4)
library(lubridate)

library(ncdf4)

# Fonction pour détecter les dates corrompues
detect_corrupt_dates_nc <- function(nc_path) {
  
  # Ouvrir le fichier NetCDF
  nc_data <- nc_open(nc_path)
  
  # Vérifier la présence de la dimension 'time'
  if (!"time" %in% names(nc_data$dim)) {
    stop("Aucune dimension 'time' trouvée dans le fichier NetCDF.")
  }
  
  # Extraire les valeurs de la dimension 'time'
  time_vals <- ncvar_get(nc_data, "time")
  
  # Récupérer les unités de temps pour la conversion (heure depuis 1950-01-01 00:00:00)
  time_units <- ncatt_get(nc_data, "time", "units")$value
  
  # Convertir les valeurs de temps en dates
  time_dates <- as.POSIXct(time_vals * 3600, origin = "1950-01-01", tz = "UTC")
  
  # Vérifier les dates corrompues (ex : dates avant 1950 ou dans un futur lointain)
  corrupt_dates <- time_dates[is.na(time_dates) | time_dates < as.POSIXct("1950-01-01", tz = "UTC")]
  
  # Fermer le fichier NetCDF
  nc_close(nc_data)
  
  return(corrupt_dates)
}

# Appel de la fonction avec un fichier NetCDF
corrupt_dates <- detect_corrupt_dates_nc("sat/sst.nc")




detect_corrupt_dates_nc("sat/sst.nc")
