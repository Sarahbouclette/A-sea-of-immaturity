#' check_scientific_names.R
#'
#' @param data A data frame with species in rows and a column "species_name", with names in this format: "Mola mola"
#'
#' @return
#' @export
#'
#' @examples


#data = species_list
#var <- sample(unique(data$species_name), 30)
# check_ind <- synonyms(str_replace("Acanthobrama_terraesanctae", "_", " "),version = "19.04")
#i = 5908
#

code_sp_check <- function(data, mc_cores = 1) {
  var <- unique(data$species_name)
  
  check <- pbmcapply::pbmclapply(1:length(var), function(i) {
    if (i %% 500 == 0) {
      Sys.sleep(sample(1:60, 1)) # Avoid API request limits
    }
    
    # Prepare species name
    if (is.na(stringr::word(var[i], 3, sep = " "))) {
      sp <- var[i]
    } else {
      sp <- paste0(stringr::word(var[i], 1, sep = " "),
                   " ", stringr::word(var[i], 3, sep = " "))
    }
    
    print(paste0(i, "/", length(var), ", ", round(i / length(var), 3) * 100, "%"))
    
    # Query synonyms from rfishbase
    check_ind <- tryCatch({
      rfishbase::synonyms(sp)
    }, error = function(e) {
      NULL # Handle errors in API request
    })
    
    # Handle cases where check_ind is empty or null
    if (is.null(check_ind) || nrow(check_ind) == 0) {
      return(data.frame(
        species_name = var[i],
        spec_code = NA,
        fishbase_name = NA,
        check = 0
      ))
    }
    
    # Filter and process results
    check_ind <- check_ind[check_ind$Status %in% c("accepted name", "synonym",
                                                   "provisionally accept",
                                                   "provisionally accepted name",
                                                   "ambiguous synonym"), ]
    
    if (nrow(check_ind) > 1 && "accepted name" %in% check_ind$Status) {
      check_ind <- check_ind[check_ind$Status == "accepted name", ]
    } else if (nrow(check_ind) > 1 && "provisionally accept" %in% check_ind$Status) {
      check_ind <- check_ind[check_ind$Status == "provisionally accept", ]
    } else if (nrow(check_ind) > 1 && "synonym" %in% check_ind$Status) {
      check_ind <- check_ind[check_ind$Status == "synonym", ]
    }
    
    # Check if there are no accepted names after filtering
    if (nrow(check_ind) == 0) {
      return(data.frame(
        species_name = var[i],
        spec_code = NA,
        fishbase_name = NA,
        check = 0
      ))
    }
    
    # Ensure correct handling when there's only one row to access
    spec_code <- ifelse(nrow(check_ind) > 0, check_ind$SpecCode[1], NA)
    fishbase_name <- ifelse(nrow(check_ind) > 0, check_ind$Species[1], NA)
    
    return(data.frame(
      species_name = var[i],
      spec_code = spec_code,
      fishbase_name = fishbase_name,
      check = nrow(check_ind)
    ))
  }, mc.cores = mc_cores)
  
  # Combine the results from all iterations
  check <- do.call(rbind, check)
  
  # Merge results with the original data
  res <- merge(data, check, by = "species_name", all.x = TRUE)
  return(res)
}
