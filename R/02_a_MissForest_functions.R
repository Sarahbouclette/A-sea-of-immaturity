################################################################################
##
##  
##
## MissForest_functions.R
##
## 29/11/2022
##
## Ulysse Flandrin
##
################################################################################
 # ##-----------------Loading packages-------------------
 # pkgs <- c("here", "dplyr", "missForest", "pbmcapply", "patchwork", "ggplot2",
 #           "maditr", "slam")
 # nip <- pkgs[!(pkgs %in% installed.packages())]
 # nip <- lapply(nip, install.packages, dependencies = TRUE)
 # ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


##-------------preping dataset-------------
preping_data <- function(species_traits_df, dim, Dims){
  
  traits_data <- species_traits_df |>  
    dplyr::select(-spec_code) |> 
    tibble::column_to_rownames("species_name") |> 
    dplyr::mutate(across(where(is.character), as.factor)) 
  
  ## Check data
  factor_length = data.frame()
  for (i in 1:length(traits_data[,])){
    #If the trait is a factor the get the length
    if (is.factor(traits_data[,i])==TRUE){
      data_frame = data.frame(id = colnames(traits_data[i]),
                              length = length(unique(traits_data[,i])))
      factor_length = rbind(factor_length,data_frame)
    }
  } #OK
  
  # #IF any of the factors has over 53 categories filter it out of the data
  # if(any(factor_length>53)){
  #   
  #   over53 <- factor_length |> 
  #     dplyr::filter(length>=53)
  #   
  #   data_to_infer <- traits_data |> 
  #     dplyr::select(-(all_of(over53$id)))
  #   # warning("Some of your traits had more than 53 categories. These traits were
  #   #         filtered out during the missForest and then re-added to your data later") 
  #   
  #   #Else keep it as it is
  # }else{
  #   data_to_infer <- traits_data 
  # }
  data_to_infer <- traits_data 
  
  
  ## long format table ##
  
  traits_data_factors <- traits_data |> 
    #extract only categorial variables 
    dplyr::select(all_of(factor_length$id)) |>
    dplyr::select(-Class, -Order, -Family, -Genus) |>
    tibble::rownames_to_column("species_name")  |> 
    #long format table
    tidyr::pivot_longer(cols = IUCN_category:last_col(),
                        names_to = "variable",
                        values_to = "observed")
  
  ### Filtrer uniquement les variables numériques
  if (dim == 0) {
    var_num <- setdiff(colnames(traits_data), unique(traits_data_factors$variable))
  }
  else {
    var_num <- setdiff(colnames(traits_data), c(unique(traits_data_factors$variable), Dims)) 
  }
  
  traits_data_num <- traits_data |> 
    #extract only categorial variables 
    dplyr::select(all_of(var_num)) |>
    dplyr::select(-Class, -Order, -Family, -Genus) |>
    tibble::rownames_to_column("species_name")  |> 
    #long format table
    tidyr::pivot_longer(cols = c(TempPrefMean:last_col()),
                        names_to = "variable",
                        values_to = "observed")
  
  return(list(data_to_infer, traits_data_factors, traits_data_num, factor_length))
         
} # END of function preping_data



##-------------MCA: Evaluate inference with Missforest - Phylogeny with MCA-------------

fct_missforest_evaluation_MCA <- function(data_to_infer,
                                          traits_data_factors = traits_data_factors,
                                          traits_data_num = traits_data_num,
                                          factor_length = factor_length,
                                          model_iteration = 5,
                                          maxiter_mf = 3, #missforest parameter #10
                                          ntree_mf = 10,  #missforest parameter #100
                                          prop_NA = 0.2) { # Proportion of NA created to evaluate the model
  
  results_list <- vector("list", model_iteration) # Stocker les résultats de chaque itération
  
  for (N in 1:model_iteration) {
    message("Début de l'itération ", N)
     result <- tryCatch({
       categorials <- factor_length$id
       categorials <- categorials[!categorials %in% c("Class", "Order", "Family", "Genus")]
       Dims <- colnames(data_to_infer)[grep("Dim", colnames(data_to_infer))]
       
       phylo <- data_to_infer |> dplyr::select(all_of(Dims))
       traits <- data_to_infer |> 
         dplyr::select(-Class, -Order, -Family, -Genus, -all_of(Dims)) |>
         mutate(across(where(is.character), as.factor),
                across(where(is.integer), as.numeric))
       
       # Produce 20% of NA in the matrix
       data_withNA <- dplyr::bind_cols(
         phylo,
         missForest::prodNA(traits, noNA = prop_NA)
       )
       
       # Ensure valid columns
       valid_columns <- sapply(data_withNA, function(x) {
         if (is.factor(x)) {
           return(length(levels(x)) > 1)
         } else {
           return(TRUE)
         }
       })
       data_withNA <- data_withNA[, valid_columns]
       
       # Convert character to factor and integer to numeric
       data_withNA[] <- lapply(data_withNA, function(x) {
         if (is.character(x)) return(factor(x))
         if (is.integer(x)) return(as.numeric(x))
         if (is.numeric(x) && length(unique(x)) <= 5) return(as.factor(x)) # Convertir les entiers avec peu de valeurs uniques en facteurs
         return(x)
       })

       
       # Remove columns with all NA
       data_withNA <- data_withNA[, colSums(is.na(data_withNA)) < nrow(data_withNA)]
       
       message("Nombre de colonnes avant missForest : ", ncol(data_withNA))
       message("Nombre de lignes avant missForest : ", nrow(data_withNA))
       
       if (ncol(data_withNA) == 0 | nrow(data_withNA) == 0) {
         stop("Problème : Aucune donnée valide après filtrage des NA !")
       }
       
       #long format Data frame with NA
       data_withNA_factors <- data_withNA  |>
         dplyr::select(all_of(categorials)) |>
         tibble::rownames_to_column("species_name")  |> 
         tidyr::pivot_longer(cols = c(IUCN_category:last_col()),
                             names_to = "variable",
                             values_to = "obs_NA")
        
       ### Filtrer uniquement les variables numériques
       var_num <- setdiff(colnames(data_withNA), c(unique(data_withNA_factors$variable), Dims))
       
       
       data_withNA_num <- data_withNA  |>
         dplyr::select(all_of(var_num)) |>
         tibble::rownames_to_column("species_name")  |> 
         tidyr::pivot_longer(cols = c(TempPrefMean:last_col()),
                             names_to = "variable",
                             values_to = "obs_NA") 
       
       
       
       graphics.off()  # Ferme toutes les fenêtres graphiques ouvertes
       
       ## Imputing data ##
       impute <- missForest::missForest(data_withNA, 
                                        maxiter = maxiter_mf, ntree = ntree_mf,
                                        variablewise = TRUE, verbose = TRUE)
       
       message("Imputation terminée pour l'itération ", N)
       print(impute$error)
       
       # Long format Data frame with imputed values 
       imputed_factors <- impute$ximp |>
         dplyr::select(all_of(categorials)) |>
         tibble::rownames_to_column("species_name") |>
         tidyr::pivot_longer(cols = c(IUCN_category:last_col()),
                             names_to = "variable",
                             values_to = "imputed")
       
       imputed_num <- impute$ximp |>
         dplyr::select(all_of(var_num)) |>
         tibble::rownames_to_column("species_name") |>
         tidyr::pivot_longer(cols = c(TempPrefMean:last_col()),
                             names_to = "variable",
                             values_to = "imputed")
       
       ### Model evaluation ###
       ## Factors imputation
       eval_factors <- traits_data_factors |>
         dplyr::left_join(data_withNA_factors) |>
         dplyr::left_join(imputed_factors) |>
         dplyr::mutate(missForest = N) |>
         dplyr::filter(is.na(obs_NA) & !is.na(observed))
       
       ## Numeric imputation
       eval_num <- traits_data_num |>
         dplyr::left_join(data_withNA_num) |>
         dplyr::left_join(imputed_num) |>
         dplyr::mutate(missForest = N) |>
         dplyr::filter(is.na(obs_NA) & !is.na(observed))
       
       list(eval_factors, eval_num)
       
     }, error = function(e) {
       message("Erreur à l'itération ", N, ": ", e$message)
       return(NULL)
     })
     
     results_list[[N]] <- result
  }
  
  return(results_list)
}




## ----------------------- Missforest application -----------------

#' missforest_applied
#'
#' @param data_to_infer the species X traits dataframe with some missing values
#' @param factor_length a two column df with the names of factoral variables of 'data_to_infer'
#'                      and the number of categories inside each.
#' @param traits_data_factors same as 'data_to_infer' but in long format with only factor variables
#' @param traits_data_num same as 'data_to_infer' but in long format with only numerical variables
#' @param var_to_infer a list of the variable we want to infer
#' @param confidence_threshold a threshold of consistency between the different missforest: 
#' proportion of missforest with the same result for factors, and 1-sd/median for numericals.
#' @param model_iteration number of missforest to do, to test the consistency of the model
#' @param maxiter_mf max number of iteration in each missforest 
#' @param ntree_mf number of tree in each missforest 
#'
#' @return the species X traits dataframe with inferred values by missforest for
#'  the selected variables, and values for which most models (>confidence_threshold)
#'  converged towards the same result
#'
#' @examples


missforest_applied <- function(data_to_infer, 
                               factor_length,
                               traits_data_factors,
                               traits_data_num,
                               var_to_infer = c("Length", "K", "IUCN_category", "trophic_guild"),
                               confidence_threshold = 0.8,
                               model_iteration = 2,
                               maxiter_mf = 1, # missForest parameter
                               ntree_mf = 10) { # missForest parameter
  
  res_missforest <- vector("list", model_iteration)  # Stocker les résultats de chaque itération
  
  for (N in seq_len(model_iteration)) {
    message("Début de l'itération ", N)
    
    result <- tryCatch({
      
      Dims <- grep("Dim", colnames(data_to_infer), value = TRUE)
      categorials <- setdiff(factor_length$id, c("Class", "Order", "Family", "Genus"))
      
      ## Imputation des données ##
      impute <- missForest::missForest(data_to_infer, 
                                       maxiter = maxiter_mf, 
                                       ntree = ntree_mf,
                                       variablewise = TRUE, 
                                       verbose = TRUE)
      
      ## Transformation en format long ##
      imputed_factors <- impute$ximp %>%
        dplyr::select(all_of(categorials)) %>%
        tibble::rownames_to_column("species_name") %>% 
        tidyr::pivot_longer(cols = -species_name,
                            names_to = "variable",
                            values_to = paste0("imputed_", N))
      
      var_num <- setdiff(colnames(impute$ximp), c(unique(imputed_factors$variable), Dims))
      
      imputed_num <- impute$ximp %>%
        dplyr::select(all_of(var_num)) %>%
        tibble::rownames_to_column("species_name") %>% 
        tidyr::pivot_longer(cols = -species_name,
                            names_to = "variable",
                            values_to = paste0("imputed_", N))
      
      ## Évaluation des facteurs imputés ##
      eval_factors <- traits_data_factors %>%
        dplyr::left_join(imputed_factors, by = c("species_name", "variable")) %>%
        dplyr::filter(is.na(observed) & variable %in% var_to_infer) %>%
        dplyr::select(-observed)
      
      ## Évaluation des variables numériques imputées ##
      eval_num <- traits_data_num %>%
        dplyr::left_join(imputed_num, by = c("species_name", "variable")) %>%
        dplyr::filter(is.na(observed) & variable %in% var_to_infer) %>%
        dplyr::select(-observed)
      
      list(eval_factors = eval_factors, eval_num = eval_num)
      
    }, error = function(e) {
      message("Erreur lors de l'itération ", N, ": ", e$message)
      return(NULL)
    })
    
    res_missforest[[N]] <- result
  }
  
  ## Vérification des résultats
  if (all(sapply(res_missforest, is.null))) {
    stop("Toutes les itérations ont échoué, vérifiez les erreurs.")
  }
  
  ## Sauvegarde des résultats
  save(res_missforest, file = here::here("Sea_of_immaturity", "outputs", "Missforest3_application_raw_results.Rdata"))
  
  ## Extraction des données imputées
  estimates_factors <- Reduce(function(x, y) left_join(x, y, by = c("species_name", "variable")),
                              lapply(res_missforest, `[[`, "eval_factors"))
  estimates_num <- Reduce(function(x, y) left_join(x, y, by = c("species_name", "variable")),
                          lapply(res_missforest, `[[`, "eval_num"))
  
  ## Confiance des imputations MissForest
  # Facteurs
  imputed_factors_values <- estimates_factors %>%
    select(starts_with("imputed_"))
  
  norm_factors <- apply(imputed_factors_values, 1, function(x) names(which.max(table(x))))
  agreement_factors <- apply(imputed_factors_values, 1, function(x) mean(x == names(which.max(table(x)))))
  
  final_estimation_factors <- estimates_factors %>%
    mutate(norm = norm_factors, confidence = agreement_factors) %>%
    select(-starts_with("imputed_"))
  
  # Numériques
  imputed_num_values <- estimates_num %>%
    select(starts_with("imputed_"))
  
  median_num <- apply(imputed_num_values, 1, median)
  deviation_num <- apply(imputed_num_values, 1, function(x) 1 - sd(x) / mean(x))  # 1 - Coeff Variation
  
  final_estimation_num <- estimates_num %>%
    mutate(norm = median_num, confidence = deviation_num) %>%
    select(-starts_with("imputed_")) %>%
    mutate(norm = as.character(norm))
  
  ## Intégration dans les données initiales
  final_imputation <- bind_rows(final_estimation_factors, final_estimation_num)
  infered_data <- data_to_infer %>%
    select(-any_of(Dims)) %>% # Suppression des dimensions pour éviter les conflits
    mutate(species_name = rownames(data_to_infer))
  
  for (i in seq_len(nrow(final_imputation))) {
    if (final_imputation$confidence[i] > confidence_threshold) {
      sp <- final_imputation$species_name[i]
      var <- final_imputation$variable[i]
      infered_data[infered_data$species_name == sp, var] <- final_imputation$norm[i]
    }
  }
  
  infered_data <- infered_data %>%
    mutate(across(where(is.character), as.numeric))
  
  return(infered_data)
}

