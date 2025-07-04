library(openxlsx)

path = ("~/Sea_of_immaturity/data/raw_data/00_prep_lspecies")

#Taxa from BRUVS program

info_species_BRUVS <- read.xlsx(paste(path, sep= "/", "BRUVS taxa list 2025_01_6.xlsx"), sheet = 2)
lspecies_BRUVS <- info_species_BRUVS$Taxa

#Taxa from RLS program

raw_data_RLS <- read.csv(paste(path, sep= "/", "RLS_raw.csv"))
lspecies_RLS <- unique(raw_data_RLS$species_name)

prep_lspecies <- function(list1, list2, output_file = "~/Sea_of_immaturity/outputs/01_gather_species_info/raw_species_list.RData") {
  # Validate inputs
  if (!is.vector(list1) || !is.vector(list2)) {
    stop("Both inputs must be vectors.")
  }
  
  # Define a pattern to exclude species names containing "sp", "spp", or "sp1" followed by numbers
  exclusion_pattern <- "\\bsp(\\d+)?\\b|\\bspp\\b"
  
  # Filter out unwanted species names
  filtered_list1 <- list1[!grepl(exclusion_pattern, list1, ignore.case = TRUE)]
  filtered_list2 <- list2[!grepl(exclusion_pattern, list2, ignore.case = TRUE)]
  
  # Find the species in common between the two lists
  common_species <- intersect(filtered_list1, filtered_list2)
  
  # Merge both lists into one, ensuring no duplicates
  merged_list <- unique(c(filtered_list1, filtered_list2))
  merged_list <- data.frame(species_name = merged_list)
  
  # Display results
  message("Number of common species: ", length(common_species))
  
  # Save the merged list to an RData file
  save(merged_list, file = output_file)
  message("Merged list saved to: ", output_file)
  
  raw_lspecies <- merged_list
  
  # Return the merged list
  return(raw_lspecies)
}

save(raw_lspecies, file = here::here("Documents","Sea_of_immaturity","data", "derived_data", "00_raw_species_list.Rdata"))
