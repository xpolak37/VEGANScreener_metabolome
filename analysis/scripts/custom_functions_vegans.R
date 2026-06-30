suppressMessages(suppressWarnings({
  library(ggplot2)
  library(cluster)   
  library(factoextra)
  library(clusterCrit) 
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
  library(emmeans)
  library(igraph)
  library(ggraph)
  library(pheatmap)
  library(ggvenn)
  library(glmnet)
  library(data.table)
  library(skimr)
  library(openxlsx)
  library(tableone)  
  
}))


# Plots PCA/sPCA scores colored by a chosen metadata variable, with custom palettes per variable
pca_plot <- function(df,metadata,color_col,perform_pca=TRUE,
                     axis=c("PC1","PC2")){
  load("../results/preprocessing/spca_result_vegans_INT_v3.RData")
  axis_num <- c(as.numeric(regmatches(axis[1], regexpr ("\\d+" , axis[1] ))),
                as.numeric(regmatches(axis[2], regexpr ("\\d+" , axis[2] ))))
  
  expl_var <- spca_result$prop_expl_var$X[axis_num] *100
  
  colors = switch(  
    color_col,
    "group"= 
      c("VG" = "#1b9e77", "OM" = "#7570b3"), 
    "country"= c(
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
    "VS_total_nosup" = c("#D73027", "#1A9850"),
    "VS_total_plus" = c("#D73027", "#1A9850"),
    "VS_total_minus" = c("#D73027", "#1A9850"),
    "vegan_duration" = c("#deebf7", "darkblue"),
    "duration_vegan_cat"= c("lightgreen","darkgreen"),
    "tertile_label" = c("#D73027", "#1A9850"),
    "hpdi_an_avg" = c("#D73027", "#1A9850"),
    "updi_an_avg" = c("#1A9850", "#D73027"),
    "pdi_an_avg" = c("#D73027", "#1A9850"),
    "UPF_perc_avg" = c("#1A9850", "#D73027"),
    NULL
  )  
  
  
  if (perform_pca){
    data.pca <- prcomp(df,scale. = TRUE,center = TRUE)
    imp_vec <- data.pca$sdev^2
    x_lab = paste(axis[1], "(",round(imp_vec[axis_num[1]],2),"%", ")", sep="")
    y_lab = paste(axis[2], "(",round(imp_vec[axis_num[2]],2),"%", ")", sep="")
    pca_df  <- as.data.frame(data.pca$x)
  } else {
    x_lab = paste(axis[1], " (",round(expl_var[1],2),"%", ")", sep="")
    y_lab = paste(axis[2], " (",round(expl_var[2],2),"%", ")", sep="")
    pca_df <- df
  }
  
  pca_df <- merge(pca_df %>% rownames_to_column("study_id"),
                  metadata,
                  by="study_id")
  
  p <- ggplot(pca_df) + 
    geom_point(aes(x=!!sym(axis[1]),y=!!sym(axis[2]),color=!!sym(color_col)),size=3) +
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
                      CI=FALSE, color="#FF6347"){
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
      CI=CI,
      color=color) 
  } else if (family=="multinomial") {
    p_value <- calc_model_p_value(model_object)
    if (p_value == 0) p_value <- "< 0.001"
    else if (p_value > 1) p_value <- "= 0.99"
    else p_value <- paste0("= ",p_value)
    roc_c <- make_multi_roc_curve(model_object,CI=CI) + 
      ggtitle(paste0(
        "Country",
        ' (AUC = ', round(model_object$model_summary$auc,3), 
        ', P ',p_value,
        ')'
      ))  
  }
  
  roc_c <- roc_c + 
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 inherit.aes = FALSE, linetype = "dashed", color = "grey")
  return(roc_c)
}

# Builds a binary-classification ROC curve, optionally with bootstrap CI ribbon across k-fold ROC objects
make_bi_roc_curve <- function(model_object,group,CI=FALSE,
                              color="#FF6347"){
  if (CI){
    p <- ggplot()
    
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
    
    p_value <- calc_model_p_value(model_object)
    auc_mean <- round(model_object$model_summary$auc,2)
    auc_cil <- round(model_object$model_summary$auc_CIL,2)
    auc_ciu <- round(model_object$model_summary$auc_CIU,2)
    
    p_value <- ifelse(p_value<0.001,"< 0.001",ifelse(
      p_value < 0.01,"< 0.01",ifelse(
        p_value < 0.05, "< 0.05",ifelse(
          p_value > 1, "= 0.99", paste("=",round(p_value,2))
        )
      )
    ))
    
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
            axis.title = element_text(face="bold",size=12)) + 
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

# Builds a multiclass (per-country) ROC curve, optionally with bootstrap CI ribbons per class
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
        true_outcome <- dummy::dummy(as.data.frame(true_outcome))
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
                    alpha = 0.2,
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
                        color=Country),alpha=1) +
          theme_minimal() + 
          theme(legend.position = "right",
                axis.text=element_text(face="bold"),
                axis.title=element_text(face="bold",size=12))  
        
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

# Creates a volcano plot (effect size vs -log10 adjusted p-value) for a chosen comparison variable, coloring points by enrichment direction
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
                             "Enriched in higher VS_total",
                             "Enriched in lower VS_total")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher VS_total"="#1A9850",
                             "Enriched in lower VS_total"="#D73027")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"VS_total"))
  } else if (comparison=="VS_total_plus"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in higher VS_total(+)",
                             "Enriched in lower VS_total(+)")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher VS_total(+)"="#1A9850",
                             "Enriched in lower VS_total(+)"="#D73027")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"VS_total_plus"))
  } else if (comparison=="VS_total_minus"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in higher VS_total(-)",
                             "Enriched in lower VS_total(-)")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher VS_total(-)"="#1A9850",
                             "Enriched in lower VS_total(-)"="#D73027")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"VS_total_minus"))
  } else if (comparison=="pdi_an_avg"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in higher PDI",
                             "Enriched in lower PDI")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher PDI"="#1A9850",
                             "Enriched in lower PDI"="#D73027")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"pdi_an_avg"))
  } else if (comparison=="hpdi_an_avg"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in higher hPDI",
                             "Enriched in lower hPDI")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher hPDI"="#1A9850",
                             "Enriched in lower hPDI"="#D73027")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"hpdi_an_avg"))
  } else if (comparison=="updi_an_avg"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in higher uPDI",
                             "Enriched in lower uPDI")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher uPDI"="#D73027",
                             "Enriched in lower uPDI"="#1A9850")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"updi_an_avg"))
  } else if (comparison=="UPF_perc_avg"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in higher UPF",
                             "Enriched in lower UPF")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher UPF"="#D73027",
                             "Enriched in lower UPF"="#1A9850")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"UPF_perc_avg"))
  } else if (comparison=="VS_total_nosup"){
    color_coding_labels <- c("Insignificant",
                             "Enriched in higher VS_total",
                             "Enriched in lower VS_total")
    color_coding_colors <- c("Insignificant"="lightgrey",
                             "Enriched in higher VS_total"="#1A9850",
                             "Enriched in lower VS_total"="#D73027")
    volcano_df <- data %>%
      dplyr::filter(str_starts(variable,"VS_total_nosup"))
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
    geom_text_repel(show.legend = FALSE) + 
    theme(
      legend.title = element_blank()
    )
  
  return(p)
}

# Plots predicted vs true values for an RF regression model, with bootstrap prediction intervals per sample
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
                    model_object$predictions$outcome,method = "pearson")
  
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

# Computes root mean squared error between predicted and true outcome
calc_rmse <- function(prediction){
  actual <- prediction$outcome
  pred <- prediction$predicted_num
  rmse <- sqrt(mean((pred - actual)^2))
  return(rmse)
}

# Computes R-squared (1 - RSS/TSS) between predicted and true outcome
calc_r2 <- function(prediction){
  actual <- prediction$outcome
  pred <- prediction$predicted_num
  rss <- sum((actual - pred)^2) # Residual sum of squares
  tss <- sum((actual - mean(actual))^2) # Total sum of squares
  r2 <- 1 - (rss / tss)
}
# Merges a metabolite dataframe with metadata, keeping only samples present in both and renaming the chosen predicting column to "Group"
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

# Generates N out-of-sample bootstrap train/validation splits for multivariate (mgaussian) outcomes
sampler_multi <- function(dataX, 
                          dataY,
                          seed = 123, 
                          sample_method = 'atypboot',
                          N = 10,
                          family="mgaussian"){
  
  set.seed(seed)
  
  dataX <- as.data.frame(dataX) %>% 
    mutate(obs_id = as.character(1:nrow(dataX)))
  
  clust_var="obs_id" 
  
  dataX <- dataX %>% 
    mutate(id = dataX[[clust_var]]) %>% 
    dplyr::select(-dplyr::all_of(clust_var))
  
  if (sample_method == 'atypboot') {
    unique_ids <- unique(dataX$id)
    
    train_data <- list()
    valid_data <- list()
    
    for (i in 1:N) {
      train <- data.frame(
        id = sample(unique(dataX$id),
                    length(unique(dataX$id)),
                    replace = TRUE))
      
      temp_train <- train %>% 
        left_join(dataX, 
                  by = 'id', 
                  relationship = "many-to-many")
      
      temp_valid <- data.frame(
        id = dataX[!dataX$id %in% temp_train$id, ]$id) %>%
        left_join(dataX, by = 'id')
      
      train_outcome <- dataY[as.numeric(temp_train$id),]
      valid_outcome <- dataY[as.numeric(temp_valid$id),]
      
      temp_train <- temp_train %>% dplyr::select(-id) %>% as.matrix()
      temp_valid <- temp_valid %>% dplyr::select(-id) %>% as.matrix()
      
      train_data[[i]] <- list(temp_train,train_outcome)
      valid_data[[i]] <- list(temp_valid,valid_outcome)
    }
    return(list(train_data, valid_data))
  }
  
}

# Generates N out-of-sample bootstrap train/validation splits for classification/regression, ensuring both classes present in train/valid for binomial/multinomial outcomes
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

# Creates k-fold cross-validation train/validation index splits
cv_folds <- function(data,folds=10){
  set.seed(123)  # for reproducibility
  n <- nrow(data) 
  k <- folds        # number of folds
  
  # Create folds
  folds <- sample(rep(1:k, length.out = n))
  
  cv_splits <- lapply(1:k, function(i) {
    val_idx <- which(folds == i)
    train_idx <- setdiff(1:n, val_idx)
    list(train = train_idx, val = val_idx)
  })
  
  return(cv_splits)
}

make_rf_model <- function(data, 
                          outcome = "Group",
                          sample_method = 'atypboot',
                          clust_var=NULL,
                          N = 10, # number of bootstrap datasets
                          overfitting_check=FALSE,
                          seed = 123,
                          reuse=FALSE,
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
# Fits a GLM for a metabolite against a tested variable (country, tertile, PDI subscores, etc.) and returns emmeans-based summary/coefficients
fit_netglm <- function(metab,data,
                       tested_variable = NULL,
                       intcheck=FALSE) {
  if (tested_variable=="country"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + VS_total + vegan_duration + country"))
    
    # Try to fit model, handle errors if any
    model <- glm(formula, data = data)
    emm_tests <- emmeans(model, ~ country) %>%
      summary(infer = TRUE) %>%  # adds p-values and CI
      as.data.frame()
    
    #model.emm.s <- emmeans(model, "Country")
    #pairs(model.emm.s)
    #pwpm(model.emm.s)
    #pwpm(model.emm.s, means = FALSE, flip = TRUE,     # args for pwpm()
    #     reverse = TRUE,                             # args for pairs()
    #     side = ">", delta = 0.05, adjust = "none")  # args for test()
    
    #eff_size(model.emm.s, sigma = sigma(model), edf = 349)
    
    # Extract coefficients summary, p-values, etc.
    coef_summary <- as.data.frame(summary(emm_tests)[,c("country","emmean","p.value","lower.CL","upper.CL")])
    
    # Return a data.frame with metabolite name and coefficients
    out <- data.frame(metabolite = metab, coef_summary)
    
    pairwise_comp <- as.data.frame(emmeans(model, "country") |> 
                                     pairs(adjust = "tukey"))
    
  } else if (tested_variable=="tertile_label"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country + tertile_label"))
    
    # Try to fit model, handle errors if any
    model <- glm(formula, data = data)
    emm_tests <- emmeans(model, ~ tertile_label) %>%
      summary(infer = TRUE) %>%  # adds p-values and CI
      as.data.frame()
    
    # Extract coefficients summary, p-values, etc.
    coef_summary <- as.data.frame(summary(emm_tests)[,c("tertile_label","emmean","p.value","lower.CL","upper.CL")])
    
    # Return a data.frame with metabolite name and coefficients
    out <- data.frame(metabolite = metab, coef_summary)
  } else if (tested_variable=="pdi_groups"){
    nutri_names <- colnames(data %>% dplyr::select(PDI_fruit_juices_avg:PDI_whole_grains_avg))
    results_list <- lapply(nutri_names, 
                           function(x) fit_netglm(metab,concat_for_uni,
                                                  tested_variable = x))
    
    results_df <- do.call(rbind, results_list)
    #nutri_names <- paste(nutri_names,collapse = " + ")
    #formula <- paste0(
    #  "`", metab, "` ~ sex + bmi + age + Country + ")
    #formula <- as.formula(paste0(formula,nutri_names))
    #model <- glm(formula, data = data)
    # Tidy output
    #out <- broom::tidy(model)
    #out$Metabolite <- metab
    return(results_df)
  } else if (tested_variable %in%  colnames(data %>% dplyr::select(PDI_fruit_juices_avg:PDI_whole_grains_avg))) {
    formula <- paste0(
      "`", metab, "` ~ sex + bmi + age + country + ",tested_variable)
    formula <- as.formula(formula)
    model <- glm(formula, data = data)
    # Tidy output
    out <- broom::tidy(model)
    out$Metabolite <- metab
    out <- out %>% dplyr::filter(term==tested_variable)
    return(out)
  } else if (tested_variable=="main_groups"){
    nutri_names <- colnames(data %>% dplyr::select(alcoholic_beverages:whole_grains))
    nutri_names <- paste(nutri_names,collapse = " + ")
    formula <- paste0(
      "`", metab, "` ~ sex + bmi + age + country + ")
    formula <- as.formula(paste0(formula,nutri_names))
    model <- glm(formula, data = data)
    # Tidy output
    out <- broom::tidy(model)
    out$Metabolite <- metab
  }
  
  return(out)
  
}

# Fits a GLM for a metabolite against a tested variable and returns pairwise (Tukey-adjusted) comparisons
fit_netglm_pairwise <- function(metab,data,
                                tested_variable = NULL,
                                intcheck=FALSE) {
  if (tested_variable=="country"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + VS_total + vegan_duration + country"))
    
    # Try to fit model, handle errors if any
    model <- glm(formula, data = data)
    
    out <- as.data.frame(emmeans(model, "country") |> 
                           pairs(adjust = "tukey"))
    out$Metabolite <- metab
  } else if (tested_variable=="tertile_label"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country + tertile_label"))
    
    # Try to fit model, handle errors if any
    model <- glm(formula, data = data)
    out <- as.data.frame(emmeans(model, "tertile_label") |> 
                           pairs(adjust = "tukey"))
    out$Metabolite <- metab
    
  } else if (tested_variable=="pdi_groups"){
    nutri_names <- colnames(data %>% dplyr::select(dairy:whole_grains))
    nutri_names <- paste(nutri_names,collapse = " + ")
    formula <- paste0(
      "`", metab, "` ~ sex + bmi + age + country + ")
    formula <- as.formula(paste0(formula,nutri_names))
    model <- glm(formula, data = data)
    # Tidy output
    out <- broom::tidy(model)
    out$Metabolite <- metab
  } else if (tested_variable=="main_groups"){
    nutri_names <- colnames(data %>% dplyr::select(alcoholic_beverages:whole_grains))
    nutri_names <- paste(nutri_names,collapse = " + ")
    formula <- paste0(
      "`", metab, "` ~ sex + bmi + age + country + ")
    formula <- as.formula(paste0(formula,nutri_names))
    model <- glm(formula, data = data)
    # Tidy output
    out <- broom::tidy(model)
    out$Metabolite <- metab
  }
  
  return(out)
  
}


# Fits a GLM for a metabolite using a formula chosen based on the tested_variable (diet, country, sex, VS_total, PDI variants, tertile, etc.) and returns coefficient table
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
  } else if (tested_variable=="country"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + VS_total + vegan_duration + country"))
  } else if (tested_variable=="sex"){
    if (intcheck) {
      formula <- as.formula(paste0(
        "`", metab, "` ~ sex + bmi + age + VS_total + vegan_duration + country"))
    } else {
      formula <- as.formula(paste0(
        "`", metab, "` ~ sex + bmi + age + VS_total + vegan_duration + country"))
    }
  } else if (tested_variable=="VS_total"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country*VS_total"))
  } else if (tested_variable=="VS_total_plus"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country*VS_total_plus"))
  } else if (tested_variable=="VS_total_minus"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country*VS_total_minus"))
  } else if (tested_variable=="pdi_an_avg"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country*pdi_an_avg"))
  } else if (tested_variable=="hpdi_an_avg"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country*hpdi_an_avg"))
  } else if (tested_variable=="updi_an_avg"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country*updi_an_avg"))
  } else if (tested_variable=="UPF_perc_avg"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country*UPF_perc_avg"))
  } else if (tested_variable=="VS_total_nosup"){
    formula <- as.formula(paste0(
      "`", metab, "` ~ sex + bmi + age + vegan_duration + country*VS_total_nosup"))
  } else if (tested_variable=="tertile_label"){
    if (intcheck) {
      formula <- as.formula(paste0(
        "`", metab, "` ~ age + bmi + sex  + tertile_label"))
    } else {
      formula <- as.formula(paste0(
        "`", metab, "` ~  age + sex + bmi + country*tertile_label"))
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

# Runs fit_glm across all pairwise country comparisons and all metabolites, returning significant PCs per comparison plus a wide significance summary table
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

roc_curve_all_custom <- function(objects,model_name,legend=TRUE){
  # Generates a ROC curve as the performance metric of a model
  # This functions will put all comparisons in one plot
  # inputs:
  # objects - roc_cs objects, where plain ROC information is stored
  # Q - question type, will be used for the path where model is stored (Q1/Q2/Q3)
  # model_name - name of the model
  # legend - boolean, if legend should be generated, default TRUE
  # outputs:
  # p - roc_curves plot
  
  #print(names(objects))
  p <- ggplot()
  
  colors <- c("#999900","#2ca02c")
  
  if (length(objects)==2) colors <- c("#FDDA24","#11457E")
  else if (length(objects)==4) colors <- c(
    "orange",
    "#11457E",
    "#FDDA24",
    "lightblue"
  )
  names(colors) <- names(objects)
  
  for (i in 1:length(objects)){
    loaded_name <- load(file.path("../intermediate_files","models",names(objects)[i],paste0(model_name,".RData")))
    my_model <- get(loaded_name)
    auc_optimism_corrected <- my_model$model_summary$auc
    
    my_rocobjs <- objects[[i]]$kfold_rocobjs
    aucs <- objects[[i]]$valid_performances$auc_validation
    q1_auc <- quantile(aucs, probs = 0.025)
    q3_auc <- quantile(aucs, probs = 0.975)
    
    differences <- abs(aucs - q1_auc)
    min_auc <- which.min(differences)
    
    differences <- abs(aucs - q3_auc)
    max_auc <- which.min(differences)
    
    differences <- abs(aucs - auc_optimism_corrected)
    true_auc <- which.min(differences)
    
    c <- my_rocobjs[[true_auc]]
    ggroc_data_c <- ggroc(c)$data
    
    # Create a common set of x-axis points
    x_common <- seq(0, 1, length.out = 100)
    
    my_df <- data.frame(x=x_common)
    for (j in 1:length(aucs)){
      ggroc_data <- ggroc(my_rocobjs[[j]])$data
      ggroc_data$`1-specificity` <- ggroc_data$`1-specificity` + (1:nrow(ggroc_data))/100000000000000000000
      ggroc_data$`1-specificity`[nrow(ggroc_data)] <- 1
      y_aprox <- approx(ggroc_data$`1-specificity`, ggroc_data$sensitivity, xout = x_common)$y
      my_df <- cbind(my_df,y_aprox)
    }
    my_df <- my_df %>% column_to_rownames("x")
    q1_y <- c()
    q2_y <- c()
    q3_y <- c()
    for (j in 1:length(x_common)){
      q1_y <- c(q1_y,quantile(as.numeric(my_df[j,]), probs = 0.025,na.rm=TRUE))
      q2_y <- c(q2_y,mean(as.numeric(my_df[j,])))
      q3_y <- c(q3_y,quantile(as.numeric(my_df[j,]), probs = 0.975,na.rm=TRUE))
    }
    df <- data.frame(x = x_common, 
                     y1 = q1_y, 
                     y2 = q3_y, 
                     y3 = q2_y)
    df[is.na(df)] <- 0
    
    my_color <- colors[i]
    
    p <- p + 
      #  geom_line(data=df,aes(x=x,y = y1,color=!!my_color)) +
      geom_ribbon(data=df,aes(x =x,ymin = pmin(y1, y2), ymax = pmax(y1, y2)), fill=my_color, 
                  alpha = 0.5,show.legend = TRUE) +
      #geom_line(data=df,aes(x=x, y = y3),color=my_color,size=1.5) +
      theme_minimal() + 
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
            panel.grid = element_blank(),
            axis.ticks.x = element_line(size=0.3,color = "black"),
            axis.ticks.y = element_line(size=0.3,color="black"),
            axis.ticks.length = unit(4,"pt"),
            axis.text=element_text(face="bold"),
            axis.title = element_text(face="bold",size=15)) + 
      ylab("Sensitivity") + xlab("1-specificity") + 
      ggtitle("low VS_total vs high VS_total")
  }
  for (i in 1:length(objects)){
    loaded_name <- load(file.path("../intermediate_files","models",names(objects)[i],paste0(model_name,".RData")))
    my_model <- get(loaded_name)
    my_rocobjs <- objects[[i]]$kfold_rocobjs
    aucs <- objects[[i]]$valid_performances$auc_validation
    
    q1_auc <- quantile(aucs, probs = 0.025)
    q3_auc <- quantile(aucs, probs = 0.975)
    
    differences <- abs(aucs - q1_auc)
    min_auc <- which.min(differences)
    
    differences <- abs(aucs - q3_auc)
    max_auc <- which.min(differences)
    
    differences <- abs(aucs - auc_optimism_corrected)
    true_auc <- which.min(differences)
    
    #a <- objects[[i]][[max_auc]]
    #b <- objects[[i]][[min_auc]]
    c <- my_rocobjs[[true_auc]]
    ggroc_data_c <- ggroc(c)$data
    y3_interp <- approx(ggroc_data_c$`1-specificity`, ggroc_data_c$sensitivity, xout = x_common)$y
    
    # Create a common set of x-axis points
    x_common <- seq(0, 1, length.out = 100)
    
    my_df <- data.frame(x=x_common)
    for (j in 1:length(aucs)){
      ggroc_data <- ggroc(my_rocobjs[[j]])$data
      ggroc_data$`1-specificity` <- ggroc_data$`1-specificity` + (1:nrow(ggroc_data))/100000000000000000000
      ggroc_data$`1-specificity`[nrow(ggroc_data)] <- 1
      y_aprox <- approx(ggroc_data$`1-specificity`, ggroc_data$sensitivity, xout = x_common)$y
      my_df <- cbind(my_df,y_aprox)
    }
    my_df <- my_df %>% column_to_rownames("x")
    q1_y <- c()
    q2_y <- c()
    q3_y <- c()
    for (j in 1:length(x_common)){
      q1_y <- c(q1_y,quantile(as.numeric(my_df[j,]), probs = 0.025,na.rm=TRUE))
      q2_y <- c(q2_y,mean(as.numeric(my_df[j,])))
      q3_y <- c(q3_y,quantile(as.numeric(my_df[j,]), probs = 0.975,na.rm=TRUE))
    }
    df <- data.frame(x = x_common, 
                     y1 = q1_y, 
                     y2 = q3_y, 
                     y3 = q2_y)
    df[is.na(df)] <- 0
    
    my_color <- colors[i]
    
    p <- p + 
      geom_line(data=df,aes(x=x, y = y3),color=my_color,size=1.5) +
      theme_minimal() + 
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
            axis.ticks.x = element_line(size=0.3,color = "black"),
            axis.ticks.y = element_line(size=0.3,color="black"),
            axis.ticks.length = unit(4,"pt"), 
            panel.grid = element_blank(),
            axis.text=element_text(face="bold"),
            axis.title = element_text(face="bold",size=15)) + 
      #scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)) + 
      #scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) + 
      ylab("Sensitivity") + xlab("1-specificity") + 
      ggtitle("low VS_total vs high VS_total")
  }
  
  if (legend){
    legend_data <- data.frame(
      color = c("Short-term vegans; AUC = 0.65; CI 0.47-0.82; p = 0.1",
                "Long-term vegans; AUC = 0.81; CI 0.66-0.94; p < 0.01"),
      value = c(1,2),
      color_hex <- colors
    )
    
    # Custom legend plot
    legend_plot <- ggplot(legend_data, aes(x = 0.5, y = value, fill = color_hex)) +
      geom_tile(width = 0.3, height = 0.8) +
      geom_text(size=4,aes(label = color), hjust = 0, nudge_x = 0.2, color = "black") +
      scale_fill_identity() +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) + xlim(0, 5)  # Adjust spacing between tiles and text
    
    p <- ggarrange(p, legend_plot, ncol=1,heights = c(1,0.2)) 
  }
  
  return(p)
}

# Applies a shared minimal black-border ggplot theme to a plot
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

# Regresses each column of a matrix on covariates and returns the residuals (covariate-adjusted values)
adjust_for_covars <- function(mat, covars, meta) {
  residuals_mat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
  colnames(residuals_mat) <- colnames(mat)
  rownames(residuals_mat) <- rownames(mat)
  
  for (i in seq_len(ncol(mat))) {
    df <- data.frame(y = mat[, i], meta[, covars, drop = FALSE])
    model <- lm(y ~ ., data = df)
    residuals_mat[, i] <- resid(model)
  }
  return(as.data.frame(residuals_mat))
  
}

# Computes Spearman correlations between each PC and each metadata variable (numeric, binary, and dummy-encoded nominal)
corr_with_pcs <- function(pca_df,metadata_df){
  ### Metadata columns
  meta_vars <- colnames(metadata_df)
  meta_df <- metadata_df
  
  ### Detect variable types
  # Numeric variables
  numeric_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
  
  # Factors/characters
  cat_vars <- names(meta_df)[sapply(meta_df, function(x) is.factor(x) || is.character(x))]
  
  # Binary categorical
  binary_vars <- cat_vars[sapply(meta_df[cat_vars], function(x) length(unique(x[!is.na(x)])) == 2)]
  
  # Nominal >2 levels (cannot be ordered)
  nominal_vars <- setdiff(cat_vars, c(binary_vars))
  
  ### Process metadata
  
  # Convert binary factors to 0/1
  binary_df <- meta_df[, binary_vars, drop = FALSE] %>%
    mutate(across(everything(), ~ as.numeric(factor(.)) - 1))
  
  # Expand nominal categorical variables into dummy variables
  dummy_list <- list()
  for (v in nominal_vars) {
    vec <- factor(meta_df[[v]], exclude = NULL)  # preserve NAs
    temp_df <- data.frame(temp = vec)
    
    dummies <- model.matrix(~ temp - 1, data = temp_df,
                            contrasts.arg = list(temp = contrasts(vec, contrasts = FALSE)))
    
    # Clean column names
    colnames(dummies) <- paste0(v, "_", levels(vec))
    
    # Convert to data frame
    dummy_list[[v]] <- as.data.frame(dummies)
  }
  
  dummy_df <- bind_cols(dummy_list) %>% dplyr::select(!ends_with("_NA"))
  
  # Combine numeric, binary, and dummy variables
  clean_meta <- bind_cols(
    meta_df[, numeric_vars, drop = FALSE],
    binary_df,
    dummy_df
  )
  
  ### 4. Spearman correlation function
  cor_test_spearman <- function(x, y) {
    suppressWarnings({
      test <- cor.test(x, y, method = "spearman", use = "pairwise.complete.obs")
    })
    tibble(
      rho = as.numeric(test$estimate),
      pvalue = test$p.value
    )
  }
  
  ### Compute correlations: PCs × metadata
  cor_results <- map_df(
    names(pca_df),
    function(pc) {
      map_df(
        names(clean_meta),
        function(meta) {
          res <- cor_test_spearman(pca_df[[pc]], clean_meta[[meta]])
          tibble(
            var1 = pc,
            var2 = meta,
            rho = res$rho,
            pvalue = res$pvalue
          )
        }
      )
    }
  )
  
  return(cor_results)
}

# Saves/updates RF model performance details (AUC, accuracy, hyperparameters, p-value) into an Excel summary sheet
save_rf_details <- function(path_to_file,model_name,rf_model){
  # save results to excel
  rf_sheet <- read.xlsx(path_to_file)
  # mtry
  mtry <- rf_model$model_summary$mtry
  
  # splitrule
  splitrule <- rf_model$model_summary$splitrule
  
  # min.node.size
  minnodesize <- rf_model$model_summary$min.node.size
  
  auc_res <- paste0(round(rf_model$model_summary$auc,3),
                    "(",round(rf_model$model_summary$auc_CIL,3),
                    "-",round(rf_model$model_summary$auc_CIU,3),
                    ")")
  acc_res <- paste0(round(rf_model$model_summary$accuracy,3),
                    "(",round(rf_model$model_summary$accuracy_CIL,3),
                    "-",round(rf_model$model_summary$accuracy_CIU,3),
                    ")")
  
  p_value <- calc_model_p_value(rf_model)
  
  if (model_name %in% rf_sheet$Model_name){
    rf_sheet$mtry[rf_sheet$Model_name==model_name] <- mtry
    rf_sheet$splitrule[rf_sheet$Model_name==model_name] <- splitrule
    rf_sheet$min.node.size[rf_sheet$Model_name==model_name] <- minnodesize 
    rf_sheet$`AUC.(CI)`[rf_sheet$Model_name==model_name] <- auc_res 
    rf_sheet$`ACC.(CI)`[rf_sheet$Model_name==model_name] <- acc_res 
    rf_sheet$p_value[rf_sheet$Model_name==model_name] <- p_value 
  } else {
    column_names <- colnames(rf_sheet)
    new_row <- c(model_name,auc_res,acc_res,p_value,mtry,splitrule,minnodesize)
    rf_sheet <- rbind(rf_sheet,new_row)
    colnames(rf_sheet) <- column_names
  }
  write.xlsx(rf_sheet,path_to_file,overwrite = TRUE)
  print(paste0("Results saved to ",path_to_file))
}

# Computes a two-sided z-test p-value testing whether mean bootstrap AUC differs from 0.5 (chance level)
calc_model_p_value <- function(model_object){
  Abar <- mean(model_object$valid_performances$auc_validation)
  SE <- sd(model_object$valid_performances$auc_validation)
  z <- (Abar - 0.5) / SE
  p_z_two <- 2 * (1 - pnorm(abs(z)))
  return(p_z_two)
}


pairwise.adonis <- function(x, metadata, formula, 
                            main_variable=NULL,
                            cluster_variable=NULL, 
                            sim.method = 'robust.aitchison', 
                            by="terms",
                            perm=999)
{
  # Runs adonis() function (PERMANOVA) for each pairwise comparison of groups in input
  # inputs:
  # x - asv table with SeqID as rownames and samples as colnames
  # metadata - metadata with SampleIDs and columns that you want to include
  # formula - formula for the adonis() function
  # main_variable - what is the name variable you want to analyze (e.g. group)
  # cluster_variable - name of the variable you want to account for (in custom permutations)
  # sim.method - distance metric, default robust.aitchison
  # perm - number of permutations (default 999), if cluster_variable is provided, custom permutations will be generated,
  # outputs:
  # ad_list: results for each pairwise comparison
  
  set.seed(123)
  
  # create df of pairwise comparison of main variable
  metadata[,main_variable] <- as.character(metadata[,main_variable])
  groups <- as.character(unique(metadata[,main_variable]))
  groups <- groups[!is.na(groups)]
  co <- t(as.data.frame(t(combn(groups, 2))))
  
  # adjust formula
  formula <- gsub("(.+) ~","x1 ~",formula)
  formula <- as.formula(formula)
  
  # create list for results
  ad_list <- list()
  
  # do comparison for each column of predefined df
  for(elem in 1:ncol(co)){
    # which samples are in the groups
    sub_samples <- metadata$SampleID[metadata[[main_variable]] %in% c(co[1,elem],co[2,elem])]
    x_sub <- x[,as.character(sub_samples)]
    # calculate distance matrix
    x1 = vegdist(x_sub %>% t(),method=sim.method)
    # subset of metadata for this comparison
    if ("SampleID" %in% colnames(metadata)) x2 = metadata %>% dplyr::filter(SampleID %in% sub_samples)
    else  x2 = metadata[sub_samples,]
    
    # create custm permutations when needed
    if (!(is.null(cluster_variable))){
      x2$Fac = as.character(x2[,main_variable])
      x2$cluster_variable <- as.character(x2[,cluster_variable])
      perm <- custom_permutations(x_sub %>% t(),
                                  factors = x2$Fac,
                                  cluster_labels = x2$cluster_variable)
    }
    
    # compute adonis2()
    ad <- adonis2(formula, data = x2,
                  permutations = perm, by=by);
    
    # create dataframe with results
    ad <- as.data.frame(ad)
    ad$sig <- NA
    ad$sig <- ifelse(ad$`Pr(>F)`<=0.001,"***",
                     ifelse(ad$`Pr(>F)`<=0.01,"**",
                            ifelse(ad$`Pr(>F)`<=0.05,"*",ad$`Pr(>F)`)))
    
    # save results to predefined list
    comp_name <- paste(c(co[1,elem],co[2,elem]),collapse=" vs ")
    ad_list[[comp_name]] <- ad
  }
  return(ad_list)
  
}

# Generates custom permutation matrices that respect clustering (e.g. same patient gets the same permutation), used inside pairwise.adonis()
custom_permutations <- function(x,factors,cluster_labels, perm=999){
  # Generates custom permutations, used in adonis()
  
  # x - data with SampleIDs as rownames
  # factors - vector of the labels
  # cluster_labels - vector of the variable labels you want to account for
  # perm - number of permutations
  
  set.seed(123)
  # create data frame
  cluster_labels_df <- data.frame(patient=cluster_labels,
                                  group=factors)
  # get unique patients
  cluster_labels_df <- cluster_labels_df[!duplicated(cluster_labels_df$patient),]
  
  # permutations across patients
  perm <- how(nperm = perm)
  perm <- shuffleSet(nrow(cluster_labels_df), control = perm)
  colnames(perm) <- cluster_labels_df$patient
  
  perm_final <- matrix(nrow=nrow(perm),ncol = length(factors)) 
  colnames(perm_final) <- rownames(x)
  
  # one patient = same permutation index
  for (j in 1:ncol(perm_final)){
    which_cluster <- cluster_labels[j]
    which_ind_cluster <- perm[,as.character(which_cluster)]
    perm_final[,j] <- which_ind_cluster
  }
  return(perm_final)
}

# Plots a heatmap-style confusion matrix (true vs predicted country) with percentage labels and highlighted diagonal
confusion_matrix_plot <- function(conf_matrix, groups){# Melt to long format
  conf_long <- reshape2::melt(t(conf_matrix))
  conf_df <- data.frame(
    True      = rep(groups, each = 5),
    Predicted = rep(groups, times = 5),
    Value     = conf_long$value
  )
  colnames(conf_df) <- c("True", "Predicted", "Value")
  
  # Add label column (percentage)
  conf_df$label <- sprintf("%.1f%%", conf_df$Value * 100)
  
  # Plot
  ggplot(conf_df, aes(x = Predicted, y = True, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label,
                  color = Value > 0.5),
              size = 4.5, fontface = "bold") +
    scale_fill_gradient2(
      low  = "#f7f7f7",
      mid  = "#d6604d",
      high = "#4d1a1a",
      midpoint = 0.35,
      limits = c(0, 1),
      name = "Proportion"
    ) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"),
                       guide = "none") +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(limits = rev(countries)) +  # <-- this fixes the diagonal
    labs(
      x = "Predicted country",
      y = "True country"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.text     = element_text(face = "bold", size = 12),
      axis.title    = element_text(face = "bold"),
      panel.grid    = element_blank(),
      legend.position = "right"
    ) +
    geom_tile(data = subset(conf_df, True == Predicted),
              aes(x = Predicted, y = True),
              fill = NA, color = "black", linewidth = 1.2)
  
  
}

direction_heatmap <- function(country_effects){
  # CLASSIFY EACH COUNTRY-METABOLITE AS HIGH / LOW / NEUTRAL 
  # Based on whether the CI is clearly above or below zero
  # (p_adj < 0.05 and CI does not cross zero)
  
  direction_df <- country_effects %>%
    mutate(
      direction = case_when(
        p_adj < 0.05 & lower.CL > 0  ~  1L,   # clearly HIGH
        p_adj < 0.05 & upper.CL < 0  ~ -1L,   # clearly LOW
        TRUE                          ~  0L    # neutral / uncertain
      )
    )
  
  # CI OVERLAP FUNCTION 
  # Returns the proportion of the shorter CI that overlaps with the other CI
  ci_overlap_prop <- function(lo1, hi1, lo2, hi2) {
    overlap <- max(0, min(hi1, hi2) - max(lo1, lo2))
    len1    <- hi1 - lo1
    len2    <- hi2 - lo2
    shorter <- min(len1, len2)
    if (shorter == 0) return(0)
    overlap / shorter
  }
  
  # FOR EACH METABOLITE: FIND GROUPS OF COUNTRIES WITH ≥50% CI OVERLAP ────
  countries <- c("BE", "CH", "CZ", "DE", "ES")
  
  group_assignments <- direction_df %>%
    group_by(metabolite) %>%
    group_modify(~ {
      df <- .x
      n  <- nrow(df)
      
      # Build overlap matrix
      overlap_mat <- matrix(0, n, n,
                            dimnames = list(df$country, df$country))
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          if (i != j) {
            overlap_mat[i, j] <- ci_overlap_prop(
              df$lower.CL[i], df$upper.CL[i],
              df$lower.CL[j], df$upper.CL[j]
            )
          }
        }
      }
      
      assigned <- rep(0L, n)
      
      for (i in seq_len(n)) {
        if (df$direction[i] == 0L) {
          assigned[i] <- 0L
          next
        }
        
        overlapping_idx <- which(overlap_mat[i, ] >= 0.50 & df$direction != 0L)
        
        if (length(overlapping_idx) == 0) {
          assigned[i] <- df$direction[i]
        } else {
          group_mean <- mean(df$emmean[c(i, overlapping_idx)])
          overall_sd <- sd(df$emmean)
          
          if (overall_sd == 0) {
            assigned[i] <- 0L
          } else {
            z <- group_mean / overall_sd
            assigned[i] <- case_when(
              z >  0.3 ~  1L,
              z < -0.3 ~ -1L,
              TRUE     ~  0L
            )
          }
        }
      }
      
      df$direction_auto <- assigned
      
      # Addition rule
      # If a non-significant country has ≥50% CI overlap with a significant one,
      # set the significant country back to neutral too
      changed <- TRUE
      while (changed) {  # repeat until no more changes (cascading neutralization)
        changed <- FALSE
        for (i in seq_len(n)) {
          if (df$direction_auto[i] == 0L) {  # neutral country (by direction_auto)
            for (j in seq_len(n)) {
              if (i != j && df$direction_auto[j] != 0L) {  # significant country
                if (overlap_mat[i, j] >= 0.50) {
                  df$direction_auto[j] <- 0L
                  changed <- TRUE  # something changed, loop again
                }
              }
            }
          }
        }
      }
      
      
      df
    }) %>%
    ungroup()
  
  
  group_assignments <- group_assignments %>%
    group_by(metabolite) %>%
    mutate(
      direction_auto = {
        vals <- direction_auto
        if (length(unique(vals)) == 1) rep(0L, n()) else vals
      }
    ) %>%
    ungroup()
  
  # KEEP ONLY METABOLITES WHERE AT LEAST ONE COUNTRY IS SIGNIFICANT 
  sig_metabolites <- group_assignments %>%
    group_by(metabolite) %>%
    filter(any(direction != 0) & any(direction_auto!=0)) %>%
    ungroup()
  
  # COMPUTE DISPERSION AND ORDER METABOLITES 
  heatmap_matrix <- sig_metabolites %>%
    select(metabolite, country, direction_auto) %>%
    mutate(direction_auto = as.numeric(direction_auto)) %>%
    pivot_wider(names_from = country, values_from = direction_auto) %>%
    tibble::column_to_rownames("metabolite")
  
  hclust_res <- hclust(dist(heatmap_matrix), method = "ward.D2")
  
  metabolite_order <- hclust_res$labels[hclust_res$order]
  
  # Country order (columns) - transpose the matrix
  hclust_ctr <- hclust(dist(t(heatmap_matrix)), method = "ward.D2")
  country_order <- hclust_ctr$labels[hclust_ctr$order]
  
  # Apply both to plot_df
  plot_df <- sig_metabolites %>%
    mutate(
      metabolite = factor(metabolite, levels = metabolite_order),
      country    = factor(country,    levels = country_order),
      dir_factor = factor(direction_auto, levels = c(-1L, 0L, 1L))
    )
  
  # plot
  p <- ggplot(plot_df, aes(x = country, y = metabolite, fill = dir_factor)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_point(
      data  = subset(plot_df, direction_auto != 0),
      aes(shape = dir_factor),
      color = "black", fill = "black", size = 2.5,
      inherit.aes = TRUE, show.legend = FALSE
    ) +
    scale_shape_manual(values = c("-1" = 25, "0" = NA, "1" = 24)) +
    scale_fill_manual(
      values = c("-1" = "#4575b4", "0" = "grey92", "1" = "#d73027"),
      labels = c("Low", "Neutral", "High"),
      name   = "Relative level"
    ) +
    labs(
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.text.x   = element_text(face = "bold", size = 12),
      axis.text.y   = element_text(size = 9),
      panel.grid    = element_blank(),
      legend.position = "right"
    )
  
  return(p)
}

# Maps cleaned/abbreviated metabolite names to nicely formatted display labels (lipid notation, chain lengths, etc.); applies them to a plot's y-axis or returns the lookup vector
fix_labels <- function(p,get_labels=FALSE){
  metabolite_labels <- c(
    # clean names (kept identical, included so nothing is dropped) 
    "Acetylornithin"                = "Acetylornithine",      
    "Ergothionenine"                = "Ergothioneine",    
    "Betaine"                       = "Betaine",
    "Asparagine"                    = "Asparagine",
    "Formate"                       = "Formate",
    "Glutamylglutamine"             = "Glutamylglutamine",
    "Hydroxydecanoylcarnitine"      = "Hydroxydecanoylcarnitine",
    "Hydroxyoctanoylcarnitine"      = "Hydroxyoctanoylcarnitine",
    "Glycine"                       = "Glycine",
    "3-Methyl-2-oxovalerate"        = "3-Methyl-2-oxovalerate",
    "Acetate"                       = "Acetate",
    "Citrate"                       = "Citrate",
    "Octenoylcarnitine"             = "Octenoylcarnitine",
    "Dodecenoylcarnitine"           = "Dodecenoylcarnitine",
    "Dimethylamine"                 = "Dimethylamine",
    "Uric.acid"                     = "Uric acid",
    "Alanine"                       = "Alanine",
    "3-Hydroxybutyrate"             = "3-Hydroxybutyrate",
    "Glucose"                       = "Glucose",
    "Hydroxybutyrylcarnitine"       = "Hydroxybutyrylcarnitine",
    "Lactate"                       = "Lactate",
    "Tyrosine"                      = "Tyrosine",
    "Inosine"                       = "Inosine",
    "Histidine"                     = "Histidine",
    "Phenylalanine"                 = "Phenylalanine",
    "Threonine"                     = "Threonine",
    "Phenylalanylphenylalanine"     = "Phenylalanylphenylalanine",
    "Methylguanidine"               = "Methylguanidine",
    "Glycerophosphocholine"         = "Glycerophosphocholine",
    "Acetone"                       = "Acetone",
    "Caffeine"                      = "Caffeine",
    "Hexanoylcarnitine"             = "Hexanoylcarnitine",
    "Phenylbutyrylglutamine"        = "Phenylbutyrylglutamine",
    "Trimethyllysine"               = "Trimethyllysine",
    "Isoleucine"                    = "Isoleucine",
    "Carnitine"                     = "Carnitine",
    "2-Oxoisocaproate"              = "2-Oxoisocaproate",
    "Propionylcarnitine"            = "Propionylcarnitine",
    "Glutamine"                     = "Glutamine",
    "Creatinine"                    = "Creatinine",
    "Valine"                        = "Valine",
    "Butyrylcarnitine"              = "Butyrylcarnitine",
    "Glutarylcarnitine"             = "Glutarylcarnitine",
    "2-Hydroxybutyrate"             = "2-Hydroxybutyrate",
    "Pyruvate"                      = "Pyruvate",
    "Lysine"                        = "Lysine",
    "Tryptophan"                    = "Tryptophan",
    "Diaminonaphthalene"            = "Diaminonaphthalene",
    "Ascorbate"                     = "Ascorbate",
    
    #  lipids: dots restored to ( ) and : 
    "PC.36.2."                      = "PC(36:2)",
    "LysoPC.18.2."                  = "LysoPC(18:2)",
    "SM.d42.2."                     = "SM(d42:2)",
    "PC.P.38.4."                    = "PC(P-38:4)",
    "PC.38.4."                      = "PC(38:4)",
    "SM.d38.2."                     = "SM(d38:2)",
    "PC.36.4."                      = "PC(36:4)",
    "SM.d34.1."                     = "SM(d34:1)",
    "LysoPC.14.0."                  = "LysoPC(14:0)",
    "PC.34.1."                      = "PC(34:1)",
    "PC.40.6."                      = "PC(40:6)",
    "LysoPC.22.6."                  = "LysoPC(22:6)",
    "PA.44.5."                      = "PA(44:5)",
    "PC.32.1."                      = "PC(32:1)",
    "LysoPC.P.16.0."                = "LysoPC(P-16:0)",
    "LysoPC.17.0."                  = "LysoPC(17:0)",
    "SM.d32.1."                     = "SM(d32:1)",
    "SM.d33.1."                     = "SM(d33:1)",
    
    # acylcarnitine chain-length annotations
    "Acylcarnitine..C8."            = "Acylcarnitine (C8)",
    "Acylcarnitine..C5."            = "Acylcarnitine (C5)",
    
    # other
    "Citrulline..M.Na."             = "Citrulline (M+Na)",
    "CE.18.2..M.NH4.."              = "CE(18:2)[M+NH4]",
    "PC.20.3.OH.P.18.1."            = "PC(20:3-OH/P-18:1)",
    "PC.P.38.4..PC.O.38.5."         = "PC(P-38:4)/PC(O-38:5)")
  
  if (get_labels) return(metabolite_labels)
  p <- p + scale_y_discrete(labels = metabolite_labels)
  return(p)
}