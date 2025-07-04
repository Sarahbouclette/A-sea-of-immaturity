################################################################################
##
##  
##
## evaluation_prediction_model.R
##
## 29/11/2022
##
## Ulysse Flandrin
##
################################################################################
# ##-----------------Loading packages-------------------
# pkgs <- c("here", "dplyr", "missForest", "pbmcapply", "patchwork", "ggplot2")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(patchwork)
library(ggplot2)

##-------------1) Extract model result -------------

extract_model_perf <- function(raw_result ){
  flat_list <- unlist(raw_result, recursive = F)
  
  
  ### Raw data ###
  raw_estimates_factor <- do.call(rbind, flat_list[seq(1, length(flat_list), 2)])
  raw_estimates_num <- do.call(rbind, flat_list[seq(2, length(flat_list), 2)])
  
  ### Extract evaluation measure per traits ###
  traits_perf = list()
  
  ## Factors
  for(i in seq(1, length(flat_list), 2)){ 
    model_i <- list()
    temp_fact <- flat_list[[i]]
    
    for(var in unique(temp_fact$variable)){
      temp_fact_trait <- temp_fact |> dplyr::filter(variable == var)
      model_number <- unique(temp_fact_trait$missForest)
      
      m_test <- mean(as.character(temp_fact_trait$observed) == as.character(temp_fact_trait$imputed))
      # Accuracy
      
      na_created <- nrow(temp_fact_trait)
      
      model_i[[var]] <- c(m_test, na_created, model_number)
    }
    
    if (length(model_i) > 0) {
      res <- as.data.frame(do.call(rbind, model_i)) |> 
        dplyr::rename(estimate = V1, NA_created = V2, model = V3) |> 
        dplyr::mutate(method = "accuracy") |> 
        tibble::rownames_to_column(var = "variable")
      traits_perf[[i]] <- res
    }
  }
  
  
  ## Numeric
  for(i in seq(2, length(flat_list), 2)){ 
    model_i <- list()
    temp_num <- flat_list[[i]]
    
    for(var in unique(temp_num$variable)){
      temp_num_trait <- temp_num |> dplyr::filter(variable == var)
      model_number <- unique(temp_num_trait$missForest)
      
      lm_test <- NA
      if (length(unique(temp_num_trait$observed)) > 1 & length(unique(temp_num_trait$imputed)) > 1) {
        lm_test <- summary(lm(temp_num_trait$observed ~ temp_num_trait$imputed))[["r.squared"]]
      }
      
      na_created <- nrow(temp_num_trait)
      
      model_i[[var]] <- c(lm_test, na_created, model_number)
    }
    
    if (length(model_i) > 0) {
      res <- as.data.frame(do.call(rbind, model_i)) |> 
        dplyr::rename(estimate = V1, NA_created = V2, model = V3) |> 
        dplyr::mutate(method = "r_squared") |> 
        tibble::rownames_to_column(var = "variable")
      traits_perf[[i]] <- res
    }
  }
  
  traits_performance <- as.data.frame(do.call(rbind, traits_perf))
  
  ### Désactivé si non utilisé
  order_performance <- NULL
  
  ### Résultats finaux
  all_res <- list(traits_performance, order_performance, raw_estimates_factor, raw_estimates_num)
  return(all_res)
}



##-------------2) Plot estimates boxplot-------------

## BOXPLOT BY TRAITS
estimates_boxplot <- function(df_estimates = traits_performance){
  #Colors
  col <- rev(fishualize::fish(n = length(unique(df_estimates$variable)), 
                          option = "Ostracion_whitleyi", begin = 0.2, end = 0.9))
  
  #order boxplot
  order_trait <- df_estimates |>
    dplyr::group_by(variable) |>
    dplyr::summarise(median_estimate = median(estimate))
  
  order_boxplot <- order_trait$variable[order(order_trait$median_estimate)]
  
  #merge data
  data <- merge(df_estimates, order_trait)
  data$variable <- factor(data$variable, levels = order_boxplot)
  
  ggplot(data) +
    aes(x= variable, y= estimate, fill = variable, col = method)+
    geom_boxplot() +
    scale_color_manual(values = c("grey", "black"))+
    scale_fill_manual(values = col)+
    xlab("") + ylab("Assessement quality (R-squared, or Accuracy)") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 10))
  
} # END of estimates_boxplot (2a)

## BOXPLOT BY ORDER
estimates_boxplot_per_order <- function(df_estimates = order_performance){
  #Colors
  col <- rev(fishualize::fish(n = length(unique(df_estimates$order)), 
                              option = "Ostracion_whitleyi", begin = 0.2, end = 0.9))
  
  ### Continuous traits
  df_estimates_num <- dplyr::filter(df_estimates, method == "r_squared")
  #order boxplot
  order_taxa <- df_estimates_num |>
    dplyr::group_by(order) |>
    dplyr::summarise(median_estimate = median(estimate, na.rm = T))
  
  order_boxplot <- order_taxa$order[order(order_taxa$median_estimate)]
  
  #merge data
  data <- merge(df_estimates_num, order_taxa)
  data$order <- factor(data$order, levels = order_boxplot)
  
  plot_r_squared <-ggplot(data) +
    aes(x= order, y= estimate, fill = order)+
    geom_boxplot(col = "black") +
    scale_fill_manual(values = col)+
    xlab("") + ylab("Assessement quality (R-squared)") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 0))
  
  ### Categorial traits
  df_estimates_trait <- dplyr::filter(df_estimates, method == "accuracy")

  #merge data
  data <- merge(df_estimates_trait, order_taxa)
  data$order <- factor(data$order, levels = order_boxplot)
  
  plot_accuracy <- ggplot(data) +
    aes(x= order, y= estimate, fill = order)+
    geom_boxplot(col = "grey") +
    scale_fill_manual(values = col)+
    xlab("") + ylab("Assessement quality (Accuracy)") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 10))
  
  plot <- plot_r_squared / plot_accuracy
  plot
  
} # END of estimates_boxplot (2a)

##-------------3) Plot estimates histogram-------------
estimates_histogramm <- function(data = traits_performance){

  plot_distribution <- function(trait, data){
    col <- fishualize::fish(n = length(unique(data$variable)), option = "Ostracion_whitleyi", begin = 0.2, end = 0.9)
    names(col) <- unique(data$variable)[order(unique(data$variable))]
    
    data <- dplyr::filter(data, variable == trait)
    
    ggplot(data) +
      aes(x = estimate) +
      geom_histogram(bins = 30,
                     fill = col[trait][[1]],
                     col = "black") +
      xlim(0,1)+
      labs(title = trait) +
      xlab("") + ylab("") +
      theme_minimal() +
      theme( legend.position = "none")+
      geom_vline(xintercept = median( data$estimate),
                 linetype="dashed", 
                 color = "black", linewidth=1) +
      annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1,
               label=paste("median = ", round(median( data$estimate), 2) ),
               col="firebrick4")

  }
  
  
  plots <- lapply(unique(data$variable)[order(unique(data$variable))], FUN = plot_distribution, data )

  vect_plots <- c()
  for(i in c(1:length(plots))){ vect_plots <- c(vect_plots, paste0("plots[[", i, "]]")) }
  code_expr <- paste(vect_plots, collapse = "+")
  
  all_plot <- eval(parse(text=code_expr)) +
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  
  all_plot
} # END of estimates_histogramm (1b)



##-------------4) study relationship between r squared and NAs-------------
estimates_NA_plot <- function(df_estimates = traits_performance){

  plot_correlation <- function(trait,data){
    col <- fishualize::fish(n = length(unique(df_estimates$variable)), option = "Ostracion_whitleyi", begin = 0.2, end = 0.9)
    names(col) <- unique(df_estimates$variable)[order(unique(df_estimates$variable))]
    
    data <- dplyr::filter(df_estimates, variable == trait)
    ggplot() +
      geom_point(aes(y = data$estimate,
                     x = data$NA_created),
                 color = col[trait][[1]], alpha = 0.6, size = 1) +
      theme_bw() +
      labs(x = "NA created", y = "Estimate", title = trait) +
      theme(panel.grid.minor = element_blank(),
            axis.text = element_text(color = "black"), 
            axis.title = element_text(size = 10))
  } 
  
  plots <- lapply(unique(df_estimates$variable)[order(unique(df_estimates$variable))], FUN = plot_correlation, data )
 
  vect_plots <- c()
  for(i in c(1:length(plots))){ vect_plots <- c(vect_plots, paste0("plots[[", i, "]]")) }
  code_expr <- paste(vect_plots, collapse = "+")
  
  all_plot <- eval(parse(text=code_expr)) +
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  
  all_plot
} # END of estimates_NA_plot (4)


##-------------5) Influence of phylogeny on estimates------------- 
estimates_phylogeny <- function(data = model_eval_missforest,                   ## check the measures + number of order different per traits...?###
                                traits = unique(traits_performance$variable),
                                level_name = FALSE){
  #color
  col <- fishualize::fish(n = length(traits), option = "Ostracion_whitleyi", begin = 0.2, end = 0.9)
  names(col) <- traits[order(traits)]

  #data
  flat_list <- unlist(data, recursive = F)
  raw_estimates_factor <- do.call(rbind,flat_list[seq(1,length(flat_list),2)])
  raw_estimates_num <- do.call(rbind,flat_list[seq(2,length(flat_list),2)])  
  
  boxplot_error_phylo <- function(trait){
    
    #trait = "IUCN_category"
    ### the trait is categorical ###
    if( trait %in% raw_estimates_factor$variable){
      
      temp_trait <- raw_estimates_factor |>
        dplyr::filter(variable == trait) |>
        dplyr::mutate(accuracy = (observed == imputed)) |> 
        dplyr::group_by(order, missForest) |>
        dplyr::summarise(prop_accuracy = mean(accuracy)*100,
                         n_sp = dplyr::n_distinct(species_name))
      
      #order boxplot
      order_taxa <- temp_trait |>
        dplyr::group_by(order) |>
        dplyr::summarise(median_accuracy = median(prop_accuracy, na.rm = T))
      
      order_boxplot <- order_taxa$order[order(order_taxa$median_accuracy)]
      
      #merge data
      data <- merge(temp_trait, order_taxa)
      data$order <- factor(data$order, levels = order_boxplot)
      

      #New labels
      replacement_labels <- lapply(order_boxplot, function(x) {
        x <- as.character(x)
        paste0(x, " (", median(temp_trait[temp_trait$order == x, "n_sp"][[1]]), ")")
      })
      
      #Boxplot
      plot <- ggplot(data)+
        geom_boxplot(aes(x = order, y = prop_accuracy),
                     color = col[trait][[1]], alpha = 0.6, linewidth = 1)+
        theme_bw() +
        labs(x = "", y = "Model accuracy (%)", title = trait) +
        theme(panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"),
              axis.title = element_text(size = 10),
              axis.text.x = element_text(angle = 45, hjust = 1,
                                         size= ifelse(level_name,10,0)))+
        
        #Add median error
        geom_hline(yintercept = median(data$prop_accuracy),
                   linetype="dashed",
                   color = "firebrick4", linewidth=1)+
        # annotate(geom = "text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
        #          label=paste("median error = ", round(median(temp_trait$prop_error), 2) ),
        #          col="firebrick4")+
        #Add the number of species in each class
        scale_x_discrete(labels = setNames(replacement_labels, order_boxplot))
    } ### end of if categorical
    
    
    #trait = "K"
    ### the trait is numerical ##
    if( trait %in% raw_estimates_num$variable){
      
      temp_trait <- raw_estimates_num |>
        dplyr::filter(variable == trait) |>
        dplyr::mutate(accuracy = 1 - abs(observed - imputed)/(max(observed)-min(observed)) ) |> 
        dplyr::group_by(order, missForest) |>
        dplyr::summarise(prop_accuracy = mean(accuracy)*100,
                         n_sp = dplyr::n_distinct(species_name))
      
      #order boxplot
      order_taxa <- temp_trait |>
        dplyr::group_by(order) |>
        dplyr::summarise(median_accuracy = median(prop_accuracy, na.rm = T))
      
      order_boxplot <- order_taxa$order[order(order_taxa$median_accuracy)]
      
      #merge data
      data <- merge(temp_trait, order_taxa)
      data$order <- factor(data$order, levels = order_boxplot)
      
      
      #New labels
      replacement_labels <- lapply(order_boxplot, function(x) {
        x <- as.character(x)
        paste0(x, " (", median(temp_trait[temp_trait$order == x, "n_sp"][[1]]), ")")
      })
      
      #Boxplot
      plot <- ggplot(data)+
        geom_boxplot(aes(x = order, y = prop_accuracy),
                     color = col[trait][[1]], alpha = 0.6, linewidth = 1)+
        theme_bw() +
        labs(x = "", y = "Model precision (%)", title = trait) +
        theme(panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"),
              axis.title = element_text(size = 10),
              axis.text.x = element_text(angle = 45, hjust = 1,
                                         size= ifelse(level_name,10,0)))+
        
        #Add median error
        geom_hline(yintercept = median(data$prop_accuracy),
                   linetype="dashed",
                   color = "firebrick4", linewidth=1)+
        # annotate(geom = "text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
        #          label=paste("median error = ", round(median(temp_trait$prop_error), 2) ),
        #          col="firebrick4")+
        #Add the number of species in each class
        scale_x_discrete(labels = setNames(replacement_labels, order_boxplot))
    } ### end of if numerical
    
    plot
  } #END OF FUNCTION boxplot_error_phylo

  plots <- lapply(traits[order(traits)], FUN = boxplot_error_phylo)
  
  vect_plots <- c()
  for(i in c(1:length(plots))){ vect_plots <- c(vect_plots, paste0("plots[[", i, "]]")) }
  code_expr <- paste(vect_plots, collapse = "+")
  
  all_plot <- eval(parse(text=code_expr)) +
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  return(all_plot)
} # END of estimates_phylogeny (5)


##-------------6) Effect of NAs on estimates ------------- 
estimates_initial_NA <- function(NA_proportion = NA_proportion){
  data <- NA_proportion[[1]]
  
  #Colors
  col <- rev(fishualize::fish(n = 2, option = "Ostracion_whitleyi", begin = 0.2, end = 0.9))
  
  P1 <- ggplot(data) +
    geom_point(aes(x= prop_NA, y= precision, col = type)) +
    geom_smooth(aes(x= prop_NA, y= precision, col = type))+
    scale_color_manual(values = col)+
    xlab("Proportion of NAs") + ylab("Assessement quality (R-squared, or Accuracy)") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 10))
  
  data_traits <- NA_proportion[[2]]
  
  P2 <-ggplot(data_traits) +
    geom_point(aes(x= prop_NA, y= precision, col = type)) +
    geom_smooth(aes(x= prop_NA, y= precision, col = type))+
    scale_color_manual(values = col)+
    xlab("Proportion of NAs") + ylab("Assessement quality (R-squared, or Accuracy)") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 10))+
    facet_wrap(~variable)
  
  list(P1, P2)
}  #END OF FUNCTION estimates_initial_NA
  




estimates_according_NA_sp <- function(order_performance = order_performance,
                                   original_data = species_traits_final){
  # Number of initial NAs per order
  order <- unique(original_data$order)[order(unique(original_data$order))]
  ratio_NA <- list()
  
  for(i in order){
    temp <- dplyr::filter(original_data, order == i) |> 
      dplyr::select(-c(species_name:fishbase_name))
    prop_NA <- mean(is.na(temp))
    nb_species <- nrow(temp)
    
    ratio_NA[[i]] <- c(prop_NA, nb_species)
  }
  
  NA_order <- do.call(rbind, ratio_NA) 
  colnames(NA_order) <- c("NA_proportion", "species_nb")
  NA_order <- tibble::rownames_to_column(as.data.frame(NA_order), var="order")

  
  # mean performance per order
  perf_order <- list()
  for(i in order){
    temp <- dplyr::filter(order_performance, order == i) |> 
      dplyr::group_by(method) |> 
      dplyr::summarise(mean = mean(estimate, na.rm=T),
                       sd = sd(estimate, na.rm =T)) |> 
      # tidyr::pivot_wider(names_from = "method", values_from = c("mean", "sd")) |> 
      dplyr::mutate(order = i)
    
    perf_order[[i]] <- temp
  }
  
  res <- do.call(rbind, perf_order) |> 
    dplyr::left_join(NA_order)
  
  
  # Plot
  P1 <- ggplot(res)+
    geom_point(aes(x=NA_proportion, y = mean, col = method))+
    # geom_smooth(aes(x=NA_proportion, y = mean, col = method))+
    xlab("Na's proportion per order") + ylab(paste("Mean per order (on",
                                                   max(order_performance$model),
                                                   "models)")) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10))
  
  
  P2 <- ggplot(res)+
    geom_point(aes(x=log10(species_nb), y = mean, col = method))+
    xlab("log10(Number of species per order)") + ylab(paste("Mean per order (on",
                                                            max(order_performance$model),
                                                            "models)")) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.background = element_rect(color ="black"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 10))
  
  P1+P2
  } # END of estimates_according_NA_sp (6)






#-------------7) Check NA on map-------------

NA_on_map <- function(data=covariates,
                      variable = "coral_algae_10km",
                      xlim = c(-180,180),
                      ylim = c(-60,60),
                      jitter = 1,
                      lat_line = 30,
                      priority= "no"){
  
  worldmap <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')
  not_priority <- ifelse(priority=="no", "yes", "no") #set the order of the points in ggplot
  
  df <- data |> 
    dplyr::select(longitude, latitude, all_of(variable)) |> 
    dplyr::mutate(presence_of_data = ifelse(!is.na(data[,variable]), "yes", "no"))
  
  ggplot() +
    geom_sf(data = worldmap, color = NA, fill = "grey70") +
    geom_jitter(data=df[df$presence_of_data == not_priority,], 
                aes(x=longitude, y=latitude, color= presence_of_data),
                width = jitter, height = jitter)+
    geom_jitter(data=df[df$presence_of_data == priority,], 
                aes(x=longitude, y=latitude, color= presence_of_data),
                width = jitter, height = jitter)+
    scale_color_brewer(palette = "PuOr")+
    coord_sf(ylim= ylim, xlim =xlim ,expand = FALSE) +
    
    geom_hline(aes(yintercept=c(-lat_line,lat_line)), linetype = "dashed", linewidth=0.3)+
    theme_bw()+
    labs(x="Longitude", y= "Latitude", title = variable) +
    theme(axis.title.x = element_text(face = "bold",
                                      size = 15),
          axis.title.y = element_text(face = "bold",
                                      size = 15),
          axis.text = element_text(size=13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          plot.title = element_text(size=10, face="bold"),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
  
}
