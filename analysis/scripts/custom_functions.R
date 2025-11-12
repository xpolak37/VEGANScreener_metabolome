# Functions
suppressMessages(suppressWarnings({
  library(ggplot2)
  library(cluster)    # For silhouette
  library(factoextra) # For easy plotting (optional)
  library(clusterCrit) # For internal clustering metrics (like CH index)
  library(tidyverse)
  library(mclust)
  library(umap)
  library(ggpubr)
  library(Hmisc)
  library(vegan)
  library(caret)
  library(ranger)
  library(doParallel)
  library(dummy)
  library(patchwork)
  library(pROC)
  library(ggrepel)
}))

## Visualization


umap_plot <- function(df, metadata,color_col="Cluster",
                      shape_col=NULL,
                      neighbors=10,
                      min_dist=0.1){
  
  custom.config <- umap.defaults
  custom.config$random_state <- 123
  custom.config$n_neighbors <- neighbors
  custom.config$min_dist <- min_dist
  
  umap_plot <- umap(df,config=custom.config)
  data_umap <- data.frame(umap_plot$layout)
  data_umap <- merge(data_umap %>% rownames_to_column("study_id"), 
                     metadata,by="study_id")
  
  if (is.null(shape_col)){
    p <- ggplot(data=data_umap, aes(x=X1,y=X2,color=!!sym(color_col))) + 
      geom_point() + 
      theme_bw() +
      xlab("UMAP 1") + ylab("UMAP 2") 
  } else {
    p <- ggplot(data=data_umap, aes(x=X1,y=X2,color=!!sym(color_col),shape=!!sym(shape_col))) + 
      geom_point() + 
      theme_bw() +
      xlab("UMAP 1") + ylab("UMAP 2") 
  }
  
  
  return(p)
}

pca_plot <- function(df,metadata,color_col,perform_pca=TRUE,
                     axis=c("PC1","PC2")){
  axis_num <- c(as.numeric(regmatches(axis[1], regexpr ("\\d+" , axis[1] ))),
                as.numeric(regmatches(axis[2], regexpr ("\\d+" , axis[2] ))))
  colors = switch(  
    color_col,
    "eq_included_group"= 
      c("vegan" = "#1b9e77", "omnivore" = "#7570b3"), 
    "Country"= c(
      BE = "#FDDA24",   # gold
      CZ = "#11457E",   # blue
      DE = "#000000",   # black
      ES = "#630356",   # purple
      CH = "#D52B1E"    # red
    ),
    "age" = c("#deebf7", "darkblue"),
    "sex" = c(
      M = "#1f77b4",     # blue
      F = "#e377c2",   # pink/purple
      OTHER = "#7f7f7f"     # gray
    ),
    "bmi" = c("#ffffb2", "#bd0026"),
    "residence_area"=c(
      "Rural" = "darkgreen",   # green — nature, countryside
      "Urban" = "navy"   # purple — modern, built environment
    ),
    "VS_total" = c("#D73027", "#1A9850"),
    "eq_vegan_dur_total" = c("#deebf7", "darkblue"),
    "tertile_label" = c("#D73027", "#1A9850"),
    NULL
  )  
  
  if (perform_pca){
    data.pca <- prcomp(df,scale. = TRUE,center = TRUE)
    var_explained <- (data.pca$sdev)^2
    
    # proportion of variance explained
    prop_var <- var_explained / sum(var_explained)
    
    # cumulative variance explained
    imp_vec <- cumsum(prop_var) *100
    x_lab = paste(axis[1], "(",round(imp_vec[axis_num[1]],2),"%", ")", sep="")
    y_lab = paste(axis[2], "(",round(imp_vec[axis_num[2]],2),"%", ")", sep="")
    pca_df  <- as.data.frame(data.pca$x)
  } else {
    x_lab <- axis[1]
    y_lab <- axis[2]
    pca_df <- df
  }
  
  pca_df <- merge(pca_df %>% rownames_to_column("study_id"),
                  metadata,
                  by="study_id")
  
  p <- ggplot(pca_df) + 
    geom_point(aes(x=!!sym(axis[1]),y=!!sym(axis[2]),color=!!sym(color_col))) +
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          axis.text = element_text(face = "bold",colour = "black"),
          axis.title = element_text(face = "bold",colour = "black"),
          panel.grid = element_blank(),
          legend.position = "right")  + 
    xlab(x_lab) + 
    ylab(y_lab) 
  
  if (is.character(pca_df[,color_col]) & !is.null(colors)) {
    p <- p + 
      scale_color_manual(values=colors) + 
      stat_ellipse(aes(x=!!sym(axis[1]),y=!!sym(axis[2]),
                       color=!!sym(color_col)),
                   show.legend = FALSE)
  } else if (!is.null(colors)){
    p <- p + scale_color_gradient(low = colors[1],high = colors[2]) 
  } 
  
  return(p)
}

roc_curve <- function(model_object, 
                      group,
                      CI=FALSE){
  # Generates a ROC curve as the performance metric of a model 
  # inputs:
  # model_object - model object
  # group - names of groups that were compared, e.g. c("rPSC","non-rPSC")
  # outputs:
  # p - roc_curve
  
  family <- ifelse(length(unique(model_object$prediction$outcome))==2,
                   "binomial",
                   "multinomial")
  
  if (family=="binomial"){
    roc_c <- make_bi_roc_curve(
      model_object,
      group,
      CI=CI) 
  } else if (family=="multinomial") {
    roc_c <- make_multi_roc_curve(model_object,CI=CI) + 
      ggtitle(paste0(
        "Country",
        ' (AUC = ', round(model_object$model_summary$auc,3), 
        ', P = ',(mean(model_object$valid_performances$auc_validation<0.5)*3),
        ')'
      ))  
  }
  
  roc_c <- roc_c + 
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 inherit.aes = FALSE, linetype = "dashed", color = "grey")
  return(roc_c)
}

make_bi_roc_curve <- function(model_object,group,CI=FALSE){
  if (CI){
    p <- ggplot()
    color <- "#FF6347"
    
    # load the object 
    auc_val <- model_object$model_summary$auc
    
    aucs <- sapply(model_object$kfold_rocobjs,function(df) df$auc)
    q1_auc <- quantile(aucs, probs = 0.025)
    q3_auc <- quantile(aucs, probs = 0.975)
    
    differences <- abs(aucs - q1_auc)
    min_auc <- which.min(differences)
    
    differences <- abs(aucs - q3_auc)
    max_auc <- which.min(differences)
    
    differences <- abs(aucs - auc_val)
    true_auc <- which.min(differences)
    
    c <- model_object$kfold_rocobjs[[true_auc]]
    ggroc_data_c <- ggroc(c)$data
    
    # Create a common set of x-axis points
    x_common <- seq(0, 1, length.out = 100)
    my_df <- data.frame(x=x_common)
    
    for (j in 1:length(aucs)){
      # get ggrocdata
      ggroc_data <- ggroc(model_object$kfold_rocobjs[[j]])$data
      
      # adjust the lowest and highest numbers
      ggroc_data$`1-specificity` <- 
        ggroc_data$`1-specificity` + (1:nrow(ggroc_data))/(1e+20)
      ggroc_data$`1-specificity`[nrow(ggroc_data)] <- 1
      
      # approximation (for equal axis coordinates)
      y_aprox <- approx(
        ggroc_data$`1-specificity`, 
        ggroc_data$sensitivity, xout = x_common)$y
      
      my_df <- cbind(my_df,y_aprox)
    }
    
    my_df <- my_df %>% column_to_rownames("x")
    q1_y <- c()
    q2_y <- c()
    q3_y <- c()
    for (j in 1:length(x_common)){
      q1_y <- c(q1_y,quantile(as.numeric(my_df[j,]), 
                              probs = 0.025,na.rm=TRUE))
      q2_y <- c(q2_y,mean(as.numeric(my_df[j,])))
      q3_y <- c(q3_y,quantile(as.numeric(my_df[j,]), 
                              probs = 0.975,na.rm=TRUE))
    }
    df <- data.frame(x = x_common, 
                     y1 = q1_y, 
                     y2 = q3_y, 
                     y3 = q2_y,
                     Comparison=paste(group[1],"vs",group[2]))
    df[is.na(df)] <- 0
    
    my_color <- color
    
    p_value <- (mean(model_object$valid_performances$auc_validation<0.5)*2)
    auc_mean <- round(model_object$model_summary$auc,2)
    auc_cil <- round(model_object$model_summary$auc_CIL,2)
    auc_ciu <- round(model_object$model_summary$auc_CIU,2)
    if ((auc_ciu) == 1) auc_ciu <- 0.99
    if (p_value<0.001) p_value <- "< 0.001"
    roc_c <- ggplot() + 
      geom_ribbon(data=df,
                  aes(x =x,ymin = pmin(y1, y2), ymax = pmax(y1, y2),
                      fill=Comparison), 
                  alpha = 0.5,
                  show.legend = TRUE) +
      geom_line(data=df,aes(x=x, y = y3, color=Comparison),
                size=1.5) + 
      theme_minimal() + 
      theme(panel.border = 
              element_rect(color = "black", fill = NA, size = 0),
            panel.grid = element_blank(),
            axis.ticks.x = 
              element_line(size=0.3,color = "black"),
            axis.ticks.y = 
              element_line(size=0.3,color="black"),
            axis.ticks.length = unit(4,"pt"),
            axis.text=element_text(face="bold"),
            axis.title = element_text(face="bold",size=15)) + 
      ylab("Sensitivity") + xlab("1-specificity") + 
      scale_color_manual(values=color) + 
      scale_fill_manual(values=color) + 
      ggtitle(paste0(
        group[1],' vs ',group[2],
        ' (AUC = ', auc_mean, ", CI: ", auc_cil,"-",auc_ciu, 
        ', p ',p_value,
        ')'
      ))  
  } else {
    ggroc_data <- ggroc(model_object$kfold_rocobjs)$data
    roc_c <- ggplot(data=ggroc_data) + 
      geom_line(aes(x=`1-specificity`, y=sensitivity, 
                    by=name, color="red",alpha=0.9)) +
      theme_minimal() + 
      theme(legend.position = "none",
            axis.text=element_text(face="bold"),
            axis.title = element_text(face="bold",size=15)) + 
      ggtitle(paste0(
        group[1],' vs ',group[2],
        ' (AUC = ', round(model_object$model_summary$auc,3), 
        ', P = ',(mean(model_object$valid_performances$auc_validation<0.5)*3),
        ')'
      ))  
  }
  return(roc_c)
}

make_multi_roc_curve <- function(model_object,CI=FALSE){
  groups <- model_object$mapping
  colors <-  c(
    BE = "#FDDA24",   # gold
    CZ = "#11457E",   # blue
    DE = "#000000",   # black
    ES = "#630356",   # purple
    CH = "#D52B1E"    # red
  )
  
  if (CI){
    roc_c <- ggplot()
    for (group in groups){
      rocs <- lapply(model_object$kfold_rocobjs,function(l){
        true_outcome <- groups[l$response]
        true_outcome <- dummy(as.data.frame(true_outcome))
        colnames(true_outcome) <- paste0(groups,"_true")
        predicted_outcome = l$predictor
        colnames(predicted_outcome) <- paste0(groups,"_pred")
        subdf <- cbind(true_outcome,predicted_outcome) %>%
          dplyr::select(starts_with(group))
        colnames(subdf) <- c("true_outcome","predicted_outcome")
        roc(true_outcome  ~ predicted_outcome , 
            data = subdf,
            direction = '<',
            levels = c(0, 1))
        
      })
      
      aucs <- sapply(rocs,function(x)x$auc)
      auc_val <- mean(aucs) #????? TO DO
      q1_auc <- quantile(aucs, probs = 0.025)
      q3_auc <- quantile(aucs, probs = 0.975)
      
      differences <- abs(aucs - q1_auc)
      min_auc <- which.min(differences)
      
      differences <- abs(aucs - q3_auc)
      max_auc <- which.min(differences)
      
      differences <- abs(aucs - auc_val)
      true_auc <- which.min(differences)
      
      true_roc <- rocs[[true_auc]]
      ggroc_data_true <- ggroc(true_roc)$data
      ggroc_data_true$Group <- "True"
      
      # Create a common set of x-axis points
      x_common <- seq(0, 1, length.out = 100)
      my_df <- data.frame(x=x_common)
      
      for (j in 1:length(aucs)){
        # get ggrocdata
        ggroc_data <- ggroc(rocs[[j]])$data
        
        # adjust the lowest and highest numbers
        ggroc_data$`1-specificity` <- 
          ggroc_data$`1-specificity` + (1:nrow(ggroc_data))/(1e+20)
        ggroc_data$`1-specificity`[nrow(ggroc_data)] <- 1
        
        # approximation (for equal axis coordinates)
        y_aprox <- approx(
          ggroc_data$`1-specificity`, 
          ggroc_data$sensitivity, xout = x_common)$y
        
        my_df <- cbind(my_df,y_aprox)
      }
      
      my_df <- my_df %>% column_to_rownames("x")
      q1_y <- c()
      q2_y <- c()
      q3_y <- c()
      for (j in 1:length(x_common)){
        q1_y <- c(q1_y,quantile(as.numeric(my_df[j,]), 
                                probs = 0.025,na.rm=TRUE))
        q2_y <- c(q2_y,mean(as.numeric(my_df[j,])))
        q3_y <- c(q3_y,quantile(as.numeric(my_df[j,]), 
                                probs = 0.975,na.rm=TRUE))
      }
      df <- data.frame(x = x_common, 
                       y1 = q1_y, 
                       y2 = q3_y, 
                       y3 = q2_y,
                       Comparison=group)
      df[is.na(df)] <- 0
      
      my_color <- colors[group]
      
      roc_c <- roc_c + 
        geom_ribbon(data=df,
                    aes(x=x,ymin = pmin(y1, y2), ymax = pmax(y1, y2),
                        fill=Comparison), 
                    alpha = 0.5,
                    show.legend = TRUE) +
        geom_line(data=df,aes(x=x, y = y3, color=Comparison),
                  size=1.5) + 
        theme_minimal() + 
        theme(panel.border = 
                element_rect(color = "black", fill = NA, size = 0),
              panel.grid = element_blank(),
              axis.ticks.x = 
                element_line(size=0.3,color = "black"),
              axis.ticks.y = 
                element_line(size=0.3,color="black"),
              axis.ticks.length = unit(4,"pt")) + 
        ylab("Sensitivity") + xlab("1-specificity") 
    }
    roc_c <- roc_c + 
      scale_color_manual(values=colors) + 
      scale_fill_manual(values=colors)
  } else {
    p <- ggplot()
    for (i in 1:length(model_object$kfold_rocobjs)){
      true_outcome <- groups[model_object$kfold_rocobjs[[i]]$response]
      true_outcome <- dummy(as.data.frame(true_outcome))
      colnames(true_outcome) <- paste0(groups,"_true")
      predicted_outcome = model_object$kfold_rocobjs[[i]]$predictor
      colnames(predicted_outcome) <- paste0(groups,"_pred")
      subdf <- cbind(true_outcome,predicted_outcome)
      for (group in groups){
        group_subdf <- subdf %>%
          dplyr::select(starts_with(group))
        colnames(group_subdf) <- c("outcome","predicted_num")
        rocobj <- roc(outcome ~ predicted_num, 
                      data = group_subdf,
                      direction = '<',
                      levels = c(0, 1))
        
        ggroc_data <- ggroc(rocobj)$data
        ggroc_data$Country <- group
        p <- p + 
          geom_line(data=ggroc_data,
                    aes(x=`1-specificity`, y=sensitivity,
                        color=Country),alpha=0.9) +
          theme_minimal() + 
          theme(legend.position = "right")  
        
      }
    }
    roc_c <- p + 
      scale_color_manual(name = 'Country',
                         breaks = c('BE', 'CZ', 'DE',"ES","CH"),
                         values = c(
                           'BE' = '#FDDA24', 
                           'CZ' = '#11457E', 
                           'DE' = '#000000',
                           "ES" = '#630356',
                           "CH" = '#D52B1E')
      )
  }
  
  return(roc_c)
}

volcano_plot <- function(data,comparison="diet"){
  if (comparison=="diet"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in vegans",
                             "Enriched in omnivores")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in vegans"="#1b9e77",
                             "Enriched in omnivores"="#7570b3")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"eq_included_group"))
  } else if (comparison=="sex"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in males",
                             "Enriched in females")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in males"="#1f77b4",
                             "Enriched in females"="#e377c2")
    
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"sex"))
  } else if (comparison=="VS_total"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in higher VS_score",
                             "Enriched in lower VS_score")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher VS_score"="#1A9850",
                             "Enriched in lower VS_score"="#D73027")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"VS_total"))
  } else if (comparison=="tertile_label"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in 3rd tertile",
                             "Enriched in 1st tertile")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in 3rd tertile"="#1A9850",
                             "Enriched in 1st tertile"="#D73027")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"tertile_label"))
  }
  volcano_df <- volcano_df %>%
    dplyr::select(metabolite,Estimate,p_adjusted) %>%
    dplyr::mutate(p_to_plot=-log10(p_adjusted)) %>%
    dplyr::mutate(Color_coding=ifelse(p_adjusted>0.05,color_coding_labels[1],
                                      ifelse(Estimate>0,color_coding_labels[2],
                                             color_coding_labels[3]))) %>%
    dplyr::mutate(Label=ifelse(p_adjusted<0.05,metabolite,""))
  
  maximum <- max(c(abs(min(volcano_df$Estimate)), abs(max(volcano_df$Estimate))))
  
  p <- ggplot(data=volcano_df,aes(x=Estimate,y=p_to_plot,color=Color_coding,
                                  label=Label)) + 
    geom_point(shape=19, size=3) + 
    theme_bw() + 
    xlab("Effect size") + ylab("- log"[10]~" p-adjusted value") + 
    scale_color_manual(values = color_coding_colors) + 
    coord_cartesian(xlim = c(-maximum,maximum), clip="on") + 
    geom_text_repel(show.legend = FALSE)
  
  return(p)
}

rf_regression_plot <- function(model_object) {
  agg_preds <- model_object$predictions %>%
    group_by(id, outcome) %>%
    summarise(
      pred_median = median(predicted_num),
      pred_lower = quantile(predicted_num, probs = 0.05),
      pred_upper = quantile(predicted_num, probs = 0.95),
      .groups = "drop"
    ) %>%
    mutate(interval_width = pred_upper - pred_lower)
  
  
  range_all <- range(c(agg_preds$outcome, agg_preds$pred_lower, agg_preds$pred_upper), na.rm = TRUE)
  
  corr_value <- cor(model_object$predictions$predicted_num,
                    model_object$predictions$outcome)
  
  p <- ggplot(agg_preds, 
              aes(x = outcome, y = pred_median, 
                  color = interval_width)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = pred_lower, ymax = pred_upper), width = 0.3) +
    scale_color_viridis_c() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = paste("RF regression, rho =",round(corr_value,2)),
      x = "True value",
      y = "Predicted value",
      color = "Interval width"
    ) +
    theme_minimal() + 
    coord_fixed(ratio = 1, xlim = range_all, ylim = range_all) 
  
  return(p)
}


## RF



calc_rmse <- function(prediction){
  actual <- prediction$outcome
  pred <- prediction$predicted_num
  rmse <- sqrt(mean((pred - actual)^2))
  return(rmse)
}

calc_r2 <- function(prediction){
  actual <- prediction$outcome
  pred <- prediction$predicted_num
  rss <- sum((actual - pred)^2) # Residual sum of squares
  tss <- sum((actual - mean(actual))^2) # Total sum of squares
  r2 <- 1 - (rss / tss)
}
model_preparation <- function(metabo_df,metadata,predicting_col="Group"){
  
  if(predicting_col=="sex"){
    metadata <- metadata %>%
      dplyr::filter(sex %in% c("F","M"))
  }
  
  if (!("study_id" %in% colnames(metabo_df))){
    metabo_df <- metabo_df %>% rownames_to_column("study_id")
  }
  
  metabo_df <- metabo_df %>%
    dplyr::filter(study_id %in% metadata$study_id)
  
  if (nrow(metadata)!=nrow(metabo_df)){
    message("Check the dimensions of provided datasets.")
    return(NULL)
  }
  
  
  # keep only study_id column and desired column to be predicted
  metadata <- metadata %>% 
    dplyr::mutate(Group=!!sym(predicting_col)) %>%
    dplyr::select(study_id,Group) 
  
  
  prepared_df <- merge(metabo_df,
                       metadata,
                       by="study_id")
  
  
  return(prepared_df)
  
}

sampler <- function(data, 
                    outcome,
                    seed = 123, 
                    clust_var=NULL,
                    sample_method = 'atypboot',
                    N = 10,
                    family="binomial"){
  
  set.seed(seed)
  
  data <- data %>% 
    mutate(obs_id = as.character(1:nrow(data)))
  
  if (is.null(clust_var)) clust_var="obs_id" 
  if (clust_var != 'id'){
    data <- data %>% 
      mutate(id = data[[clust_var]]) %>% 
      dplyr::select(-dplyr::all_of(clust_var))
  }
  
  if (colnames(data[outcome]) != 'outcome'){
    data <- data %>% 
      mutate(outcome = data[[outcome]]) %>% 
      dplyr::select(-dplyr::all_of(c(outcome)))
  }
  
  data <- data %>% 
    mutate(obs_id = as.character(1:nrow(data))) %>%
    dplyr::select(id,dplyr::everything())
  
  if (sample_method == 'atypboot') {
    unique_ids <- unique(data$id)
    
    train_data <- list()
    valid_data <- list()
    
    for (i in 1:N) {
      repeat {
        train <- data.frame(
          id = sample(unique(data$id),
                      length(unique(data$id)),
                      replace = TRUE))
        
        temp_train <- train %>% 
          left_join(data, 
                    by = 'id', 
                    relationship = "many-to-many")
        
        temp_valid <- data.frame(
          obs_id = data[!data$obs_id %in% temp_train$obs_id, ]$obs_id) %>%
          left_join(data, by = 'obs_id')
        
        if (family!="regression"){
          train_outcome <- as.numeric(factor(temp_train$outcome)) - 1 
          valid_outcome <- as.numeric(factor(temp_valid$outcome)) - 1
          
          if (!mean(train_outcome) %in% c(0, 1) & 
              !mean(valid_outcome) %in% c(0, 1)) {
            train_data[[i]] <- temp_train
            valid_data[[i]] <- temp_valid
            {break} 
          } 
        } else {
          train_data[[i]] <- temp_train
          valid_data[[i]] <- temp_valid
          {break} 
        }
      }
    }
    return(list(train_data, valid_data))
  }
  
}


make_rf_model <- function(data, 
                          outcome = "Group",
                          sample_method = 'atypboot',
                          clust_var=NULL,
                          N = 10, # number of bootstrap datasets
                          overfitting_check=FALSE,
                          seed = 123,
                          reuse=TRUE,
                          file=NULL) {
  
  # Fits a RF model to dataset 
  # inputs:
  # data - prepared dataframe using binomial_prep(), contains only two groups
  # outcome - vector of outcome (labels)
  # sample_method - atypboot - out-of-sample boostrap
  # clust_var, name of clustering variable (here, it can be Patient), default NULL
  # N, number of bootstrap samples, default 10 (minimum 100 needed for reportable results)
  # family - A string specifying the family for the GLMNET model, default 'binomial'
  # overfitting_check - boolean, if random labels reshuffling should be performed, default FALSE
  # seed - random seed for reproducibility, default 123
  # reuse - boolean, should model just be reloaded?, default FALSE
  # file - name of the file for loading the pre-trained model
  # Q - analysis question - Q1/Q2/Q3..
  # outputs: 
  # rf_model - list(model_summary,valid_performances,
  #                   valid_performance,
  #                   predictions, 
  #                   roc_curve,
  #                   kfold_roc_curves,
  #                   trained_model)
  
  # if we are just reusing precomputed models
  if (reuse){
    # model for overfitting check
    if (overfitting_check) {
      load(
        file.path(
          "../intermediate_files/models_overfitting_check/",
          file,
          "rf_model.RData")
      )
    }
    # normal model
    else load(
      file.path("../intermediate_files/models/",
                file,
                "rf_model.RData")
    )
  } else {
    # if we are computing the model from scratch
    set.seed(seed)
    
    # if we are doing overfitting check, lets sample the Groups to be random
    if (overfitting_check) {
      data$Group <- sample(data$Group)
    } else {
      data <- data[sample(1:nrow(data)),] %>%
        `rownames<-`(NULL)
    }
    
    data <- data.frame(data,check.names=TRUE)
    
    ## where to save relevant information
    auc_validation <- vector('double', N)
    accuracy_validation <- vector('double', N)
    predictions <- vector("list", N)
    conf_matrices <- vector("list", N)
    kfold_rocobjs <- vector("list", N)
    rmse_validation <- vector('double', N)
    r2_validation <- vector('double', N)
    
    ## original data in a matrix form
    original_outcome <- base::as.matrix(data[[outcome]])
    
    if(!is.numeric(original_outcome)){
      original_outcome <- factor(original_outcome)
      mapping <- setNames(
        levels(original_outcome), 
        0:(length(levels(original_outcome))-1))
      original_outcome <- as.numeric(factor(original_outcome)) - 1 
    }
    
    # set important variables - classification
    metric="Accuracy"
    split_rule = "gini"
    importance = "impurity"
    probability=TRUE
    family=NULL
    if (length(unique(original_outcome))==2){
      family="binomial"
    } else if (length(unique(original_outcome))>2 &
               length(unique(original_outcome))<10) {
      family="multinomial"
    } else {
      # REGRESSION
      family="regression"
      metric="RMSE"
      split_rule="variance"
      importance = "permutation"
      probability=FALSE
    }
    # prepare the predictors and outcome
    original_predictors <- data %>% 
      dplyr::select(-dplyr::all_of(c("study_id",outcome,clust_var))) %>% 
      as.matrix()
    
    if (family!="regression") {
      original_outcome <- as.factor(original_outcome)
    } else original_outcome <- as.numeric(as.character(original_outcome))
    
    tune_grid <- expand.grid(
      mtry = seq(1, ncol(original_predictors),2), # Number of predictors 
      min.node.size = c(2, 5), # Minimum node size
      splitrule = split_rule
    )
    
    fitControl <- trainControl(method = "repeatedcv",
                               ## 5-fold CV...
                               number = 5,
                               ## repeated 5 times
                               repeats = 1)
    
    # OPTIMALIZATION Train the random forest model using ranger
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)
    
    rf_model <- caret::train(
      x = original_predictors, 
      y = original_outcome, 
      method = "ranger",  # Use kknn method for more flexibility
      tuneGrid = tune_grid,  # Manhattan or Euclidean distance
      metric = metric, # Specify the evaluation metric
      trControl = fitControl,
      importance = importance
    )
    
    stopCluster(cl)
    
    # store the optimized hyperparameters
    optim_mtry = rf_model$bestTune$mtry
    optim_splitrule = rf_model$bestTune$splitrule
    optim_min.node.size = rf_model$bestTune$min.node.size
    
    # Combine predictors and outcome into a single data frame
    train_data <- data.frame(original_outcome, 
                             original_predictors, 
                             check.names = TRUE)
    
    # RF formula
    fit <- ranger(
      formula = original_outcome ~ ., # Explicitly specify the outcome 
      data = train_data,        # Use the combined data frame
      mtry = optim_mtry,  # Use the best value of mtry 
      splitrule = optim_splitrule, # splitrule
      min.node.size = optim_min.node.size, # minimum node size
      probability = probability, # store the probability,
      importance = importance
    )
    
    # predictions
    if (family=="binomial"){
      predicted_num = stats::predict(
        fit, 
        data = as.data.frame(original_predictors))$predictions[,2]
      
      predicted_orig <- ifelse(predicted_num > 0.5, 1, 0)
      
      prediction <- data.frame(
        predicted = predicted_orig,
        predicted_num=predicted_num,
        outcome = original_outcome)
      
      rocobj <- roc(outcome ~ predicted_num, 
                    data = prediction,
                    direction = '<',
                    levels = c(0, 1))
      
    } else if (family=="multinomial"){
      predicted_num = stats::predict(
        fit,
        data = as.data.frame(original_predictors))$predictions
      
      predicted_orig <- apply(predicted_num,1,function(x) which.max(x)) - 1
      
      prediction <- list(
        predicted = predicted_orig,
        predicted_num=predicted_num,
        outcome = original_outcome)
      
      rocobj <- multiclass.roc(
        original_outcome,
        predicted_num,
        levels = colnames(predicted_num))
    } else if (family=="regression"){
      predicted_num = stats::predict(
        fit,
        data = as.data.frame(original_predictors))$predictions
      
      prediction <- data.frame(
        predicted_num=predicted_num,
        outcome = original_outcome)
      
      rmse <- calc_rmse(prediction)
      r2 <- calc_r2(prediction)
      
    }
    
    if (family!="regression"){
      fitted <- data.frame(
        mtry = optim_mtry,  
        splitrule = optim_splitrule,
        min.node.size = optim_min.node.size,
        auc = rocobj$auc,
        accuracy = mean(prediction$predicted==prediction$outcome))
    } else {
      fitted <- data.frame(
        mtry = optim_mtry,  
        splitrule = optim_splitrule,
        min.node.size = optim_min.node.size,
        rmse = rmse,
        r2 = r2)
    }
    
    
    # BOOTSTRAP
    sampled_data <- sampler(data, 
                            outcome = outcome,
                            clust_var = clust_var,
                            sample_method = sample_method,
                            N = N, 
                            seed = seed,family=family)
    
    if(!is.null(clust_var)) clust_var <- "id"
    
    for (i in 1:N){
      
      ## sampled data in a matrix form
      sampled_outcome <- as.matrix(sampled_data[[1]][[i]]$outcome)
      
      if (!is.numeric(sampled_outcome)){
        sampled_outcome <- factor(sampled_outcome)
        mapping <- setNames(levels(sampled_outcome), 0:(length(levels(sampled_outcome))-1))
        sampled_outcome <- as.numeric(factor(sampled_outcome)) - 1 
      }
      
      if (family!="regression") {
        sampled_outcome <- as.factor(sampled_outcome)
      } else {
        sampled_outcome <- as.numeric(as.character(sampled_outcome))
      }
      
      sampled_predictors <- sampled_data[[1]][[i]] %>% 
        dplyr::select(-dplyr::all_of(
          c("study_id", "outcome","id", "obs_id", clust_var))) %>%
        as.matrix()
      
      ## re-optimize parameters
      # OPTIMALIZATION
      
      cl <- makePSOCKcluster(5)
      registerDoParallel(cl)
      
      rf_model_sampled <- caret::train(
        x = sampled_predictors, 
        y = sampled_outcome, 
        method = "ranger",  # rf
        tuneGrid = tune_grid,  # gini
        metric = metric, # Specify the evaluation metric
        trControl = fitControl,
        importance = importance
      )
      
      stopCluster(cl)
      
      optim_mtry_sampled = rf_model_sampled$bestTune$mtry
      optim_splitrule_sampled = rf_model_sampled$bestTune$splitrule
      optim_min.node.size_sampled = rf_model_sampled$bestTune$min.node.size
      
      # Combine predictors and outcome into a single data frame
      train_data_sampled <- data.frame(sampled_outcome, 
                                       sampled_predictors, 
                                       check.names = TRUE)
      # TRAINING
      sampled_fit <- ranger(
        formula = sampled_outcome ~ .,    # outcome variable 
        data = train_data_sampled,        # Use the combined data frame
        mtry = optim_mtry_sampled,  # Use the best value of k from tuning
        splitrule = optim_splitrule_sampled,
        min.node.size = optim_min.node.size_sampled,
        probability = probability
      )
      
      # VALIDATION
      valid_outcome <- as.matrix(
        sampled_data[[2]][[i]]$outcome)
      
      valid_predictors <- sampled_data[[2]][[i]] %>% 
        dplyr::select(-study_id, -outcome, -id, -obs_id) %>% 
        as.matrix()
      
      
      if (family=="binomial"){
        valid_outcome <- as.numeric(factor(valid_outcome)) - 1 
        valid_outcome <- as.factor(valid_outcome)
        
        valid_predicted_num = stats::predict(
          sampled_fit, 
          data = as.data.frame(valid_predictors)
        )$predictions[,2]
        
        valid_predicted_orig <- ifelse(valid_predicted_num > 0.5, 1, 0)
        
        prediction_onValidation <- data.frame(
          predicted = valid_predicted_orig,
          predicted_num=valid_predicted_num,
          outcome = valid_outcome)
        
        kfold_rocobjs[[i]] <- roc(outcome ~ predicted_num, 
                                  data = prediction_onValidation,
                                  direction = '<',
                                  levels = c(0, 1))
        
        auc_validation[i] <- kfold_rocobjs[[i]]$auc
        accuracy_validation[i] <- mean(
          prediction_onValidation$predicted == prediction_onValidation$outcome)
        
      } else if (family=="multinomial"){
        valid_outcome <- as.numeric(factor(valid_outcome)) - 1 
        valid_outcome <- as.factor(valid_outcome)
        
        valid_predicted_num = stats::predict(
          sampled_fit,
          data = as.data.frame(valid_predictors))$predictions
        
        valid_predicted_orig <- apply(valid_predicted_num,1,function(x) which.max(x)) - 1
        
        prediction_onValidation <- list(
          predicted = valid_predicted_orig,
          predicted_num=valid_predicted_num,
          outcome = valid_outcome,
          iteration = i)
        
        ## record performances
        kfold_rocobjs[[i]] <- multiclass.roc(
          valid_outcome,
          valid_predicted_num,
          levels = colnames(valid_predicted_num))
        
        auc_validation[i] <- kfold_rocobjs[[i]]$auc
        accuracy_validation[i] <- mean(
          prediction_onValidation$predicted == prediction_onValidation$outcome)
        
      } else {
        # regression
        valid_predicted_num = stats::predict(
          sampled_fit,
          data = as.data.frame(valid_predictors))$predictions
        
        prediction_onValidation <- data.frame(
          predicted_num=valid_predicted_num,
          outcome = valid_outcome,
          iteration = i)
        
        ## record performances
        rmse_validation[i] <- calc_rmse(prediction_onValidation)
        r2_validation[i]<- calc_r2(prediction_onValidation)
      }
      
      
      prediction_onValidation[["id"]] <- sampled_data[[2]][[i]]$id
      predictions[[i]] <- prediction_onValidation
      
    }
    
    ## connect predictions
    predictions <- bind_rows(predictions)
    
    ## aggregate obtained information
    valid_performances <- data.frame(
      auc_validation,
      accuracy_validation,
      rmse_validation,
      r2_validation
    )
    
    if (family!="regression"){
      model_summary <- fitted %>% 
        mutate(
          auc = mean(valid_performances$auc_validation),
          auc_CIL = quantile(valid_performances$auc_validation, 
                             probs = 0.025),
          auc_CIU = quantile(valid_performances$auc_validation, 
                             probs = 0.975),
          accuracy = mean(valid_performances$accuracy_validation),
          accuracy_CIL = quantile(valid_performances$accuracy_validation, 
                                  probs = 0.025),
          accuracy_CIU = quantile(valid_performances$accuracy_validation, 
                                  probs = 0.975)
        )
    } else {
      model_summary <- fitted %>% 
        mutate(
          rmse = mean(valid_performances$rmse_validation),
          rmse_CIL = quantile(valid_performances$rmse_validation, 
                              probs = 0.025),
          rmse_CIU = quantile(valid_performances$rmse_validation, 
                              probs = 0.975),
          r2 = mean(valid_performances$r2_validation),
          r2_CIL = quantile(valid_performances$r2_validation, 
                            probs = 0.025),
          r2_CIU = quantile(valid_performances$r2_validation, 
                            probs = 0.975)
        )
      rocobj <- NULL
      mapping <- NULL
    }
    
    
    ## define outputs
    rf_model <- list(model_summary = model_summary, 
                     valid_performances = valid_performances, 
                     prediction = prediction,
                     predictions = predictions,
                     rocobj=rocobj,
                     kfold_rocobjs=kfold_rocobjs,
                     trained_model=fit,
                     mapping=mapping)
    
    # save results
    if (overfitting_check){
      if (!dir.exists(file.path("../intermediate_files/models_overfitting_check/",file))){
        dir.create(file.path("../intermediate_files/models_overfitting_check/",file))
      } 
      
      save(rf_model,
           file=file.path("../intermediate_files/models_overfitting_check/",
                          file,
                          "rf_model.RData")
      )
    } else {
      if (!dir.exists(file.path("../intermediate_files/models/",file))){
        dir.create(file.path("../intermediate_files/models/",file))
      }
      save(rf_model,
           file=file.path("../intermediate_files/models/",file,"rf_model.RData"))
    }
  }
  
  return(rf_model)
}


## Tests


fit_glm <- function(metab,data,
                    tested_variable = NULL,
                    intcheck=FALSE) {
  if (tested_variable=="eq_included_group"){
    if (intcheck) {
      formula <- as.formula(paste0(
        "`", metab, "` ~ age + sex + bmi + eq_included_group"))
    }
    else {
      formula <- as.formula(paste0(
        "`", metab, "` ~ age + sex + bmi + Country * eq_included_group"))
      
    } 
  } else if (tested_variable=="Country"){
    formula <- as.formula(paste0(
      "`", metab, "` ~  age + sex + bmi + Country"))
  } else if (tested_variable=="sex"){
    if (intcheck) {
      formula <- as.formula(paste0(
        "`", metab, "` ~ age + bmi + eq_included_group + sex"))
    } else {
      formula <- as.formula(paste0(
        "`", metab, "` ~ age + bmi + eq_included_group * sex"))
    }
  } else if (tested_variable=="VS_total"){
    formula <- as.formula(paste0(
      "`", metab, "` ~  age + sex + bmi + Country + eq_included_group* VS_total"))
  } else if (tested_variable=="tertile_label"){
    if (intcheck) {
      formula <- as.formula(paste0(
        "`", metab, "` ~ age + bmi + sex + eq_included_group + tertile_label"))
    } else {
      formula <- as.formula(paste0(
        "`", metab, "` ~  age + sex + bmi + Country + eq_included_group* tertile_label"))
    }
  }
  
  
  # Try to fit model, handle errors if any
  model <- glm(formula, data = data)
  
  # Extract coefficients summary, p-values, etc.
  coef_summary <- as.data.frame(summary(model)$coefficients) %>% 
    rownames_to_column("variable")
  
  # Return a data.frame with metabolite name and coefficients
  out <- data.frame(metabolite = metab, coef_summary)
  return(out)
  
}

pairwise_glm <- function(spca_for_uni,tested_variable=NULL){
  co <- data.frame("1"=c("BE","BE","BE","BE","CZ","CZ","CZ","DE","DE","ES"),
                   "2"=c("CZ","DE","ES","CH","DE","ES","CH","ES","CH","CH")) %>% t()
  
  final_results_df <- data.frame()
  for (i in 1:ncol(co)){
    comp <- co[,i]
    subdata <- spca_for_uni %>% 
      dplyr::filter(Country %in% comp)
    
    results_list <- lapply(
      colnames(spca_integrated_df),
      function(x) fit_glm(x,subdata,
                          tested_variable = tested_variable))
    
    results_df <- do.call(rbind, results_list)
    results_df$Comparison <- paste0(comp[1],"_vs_",comp[2])
    results_df$p_adjusted <- p.adjust(results_df$Pr...t..,method = "BH")
    final_results_df <- rbind(final_results_df,results_df)
  }
  
  df_listed <- final_results_df %>%
    # Filter only Country rows (ignore intercept, age, sex, etc.)
    filter(grepl("^Country", variable)) %>%
    # Only keep significant p-values (adjusted)
    filter(p_adjusted < 0.05) %>%
    # Group by Comparison
    group_by(Comparison) %>%
    # Collect the PC values where significant
    summarise(Significant_PCs = paste(unique(metabolite), collapse = ", ")) %>%
    arrange(Comparison)
  
  df_summary <- final_results_df %>%
    # Filter rows for country variable only
    filter(grepl("^Country", variable)) %>%
    # Create a flag for significance
    mutate(significant = ifelse(p_adjusted < 0.05, "YES", "NO")) %>%
    # Select only relevant columns
    select(metabolite, Comparison, significant) %>%
    # Pivot wider to get PCs as columns
    pivot_wider(
      names_from = metabolite,
      values_from = significant,
      values_fill = "NO"
    ) %>%
    arrange(Comparison)
  
  return(list(df_listed,df_summary))
}

rf_regression_test <- function(model_object){
  errors_model <- abs(
    model_object$prediction$predicted_num - model_object$prediction$outcome)
  errors_null <- abs(
    mean(model_object$prediction$outcome) - model_object$prediction$outcome)
  
  # Paired test on absolute errors (or squared errors)
  return(t.test(errors_null, errors_model, paired = TRUE, alternative = "greater"))
}

ggtheme <- function(p){
  p <- p + theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, 
                                      size = 1,linewidth = 1),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          axis.text = element_text(face = "bold",colour = "black"),
          axis.title = element_text(face = "bold",colour = "black"),
          panel.grid = element_blank(),
          legend.position = "right")
  return(p)
}

