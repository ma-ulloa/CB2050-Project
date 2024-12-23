
# Utils functions 

# Load necessary packages
# Load the tidyverse package
suppressMessages(library(tidyverse))
suppressMessages(library(tidymodels))
suppressMessages(library(vip)) # Check importance of features in the model
suppressMessages(library(patchwork))
suppressMessages(library(scales))
suppressMessages(library(paletteer))



plot_roc_auc <- function(final_training_predictions, tittle) {
  pred <- final_training_predictions |> collect_predictions()
  roc_data <- roc_curve(
    pred,
    truth = Disease,
    `.pred_Myositis`, 
    `.pred_Rheumatoid arthritis`, 
    `.pred_Sjögrens syndrome`, 
    `.pred_Systemic sclerosis`, 
    `.pred_Systemisk lupus erythematosus`
  )
  
  autoplot(roc_data) +
    ggtitle(tittle) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal()
}


conf_mat_plot <- function(final_model_pred, title) {
  pred <- final_model_pred |> 
    collect_predictions() |>
    count(Disease, .pred_class) |>
    complete(Disease, .pred_class, fill = list(n = 0))
  
  # Define new labels
  labels <- c(
    "Myositis" = "IMM",
    "Sjögrens syndrome" = "SjD",
    "Systemic sclerosis" = "SSc",
    "Systemisk lupus erythematosus" = "SLE",
    "Rheumatoid arthritis" = "RA"
  )
  
  # Plot confusion matrix
  ggplot(data = pred, mapping = aes(x = Disease, y = .pred_class)) +
    geom_tile(aes(fill = n), colour = "white") +  
    geom_text(aes(label = sprintf("%1.0f", n)), vjust = 0.5, size = 5) +  
    scale_fill_gradient(low = "white", high = "#4685A0FF") +
    scale_x_discrete(labels = labels, expand = expansion(mult = 0)) +  # Apply new labels
    scale_y_discrete(labels = labels, expand = expansion(mult = 0)) +  # Apply new labels
    labs(x = "Actual Disease", y = "Predicted Disease") +
    theme_minimal() +
    theme(
      text = element_text(size = 14),                 
      axis.text = element_text(size = 12),            
      axis.title = element_text(size = 16, face = "bold"),  
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  
      panel.grid = element_blank(),                  
      legend.position = "right",                     
      legend.title = element_text(size = 14),        
      legend.text = element_text(size = 12)          
    ) +
    ggtitle(title)
}



barplot_vip <- function(vi_table, variable_y, feature_num=20, tittle, color) {
  vi_table |> 
    slice_max(order_by = Importance, n = feature_num) |>  
    ggplot(aes(x = reorder(Variable, !!sym(variable_y)), y = !!sym(variable_y))) +  
    geom_bar(stat = "identity", fill = color) +  
    coord_flip() +
    theme_classic() +
    theme(
      axis.text = element_text(size = 16),         
      axis.title = element_text(size = 16)         
    ) +
    labs(x = "Features", y = "Importance") + ggtitle(tittle)
}


lollipop_vip <- function(vi_table, variable_y, feature_num=20, tittle, color) {
  vi_table |> 
    slice_max(order_by = Importance, n = feature_num) |>  
    ggplot(aes(x = reorder(Variable, !!sym(variable_y)), y = !!sym(variable_y))) +
    geom_segment(aes(x = reorder(Variable, !!sym(variable_y)), 
                     xend = reorder(Variable, !!sym(variable_y)), 
                     y = 0, yend = !!sym(variable_y)), color = color) +
    geom_point(color = color, size = 3, alpha = 1) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text = element_text(size = 16),         
      axis.title = element_text(size = 16)         
    ) +
    labs(x = "Features", y = "Importance") + 
    ggtitle(tittle)
}


plot_vi_by_class <- function(final_model, palette, top_number, title = "Relevant Features per Class") {
  # Extract model fit and calculate feature importance
  importance_df <- final_model %>% 
    extract_fit_parsnip() %>% 
    tidy() %>% 
    filter(term != "(Intercept)") %>% 
    mutate(
      abs_estimate = abs(estimate),
      rank = row_number(-abs_estimate)
    ) %>% 
    filter(abs_estimate > 0) %>% 
    mutate(
      Scaled_estimate = rescale(abs_estimate, to = c(0, 100))
    ) %>% 
    group_by(class) %>% 
    slice_max(order_by = abs_estimate, n = top_number) %>%
    arrange(rank) 
  
  # Create the plot
  importance_plot <- importance_df %>% 
    ggplot(aes(x = reorder(term, Scaled_estimate), y = Scaled_estimate, color = class)) +
    geom_segment(aes(xend = term, yend = 0), size = 1) +
    geom_point(size = 3, alpha = 1) +
    facet_wrap(~class, scales = "free_y") +
    coord_flip() +
    scale_color_manual(values = palette) + 
    theme_bw() +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.position = "none"
    ) +
    labs(x = "Features", y = "Importance") + 
    ggtitle(title)
  
  return(importance_plot)
}

# Expression visualization with boxplots

boxplot_expression <- function(protein_name, countmatrix, class, color_palete) {
  countmatrix |> 
    select(all_of(class), all_of(protein_name)) |>
    ggplot(aes(x = !!sym(class), y = !!sym(protein_name), fill = !!sym(class))) +
    geom_violin(trim = FALSE, alpha = 0.6) +  
    geom_boxplot(width = 0.1, color = "black", alpha = 0.7, outlier.shape = NA) +  
    geom_jitter(width = 0.2, size = 1.5, color = "black", alpha = 0.7) + 
    scale_fill_manual(values = color_palete) + 
    theme_bw() +
    theme(
      axis.text = element_text(size = 16),         
      axis.title = element_text(size = 16),
      legend.position = "none"
    ) +
    labs(x = "Disease", y = "NPX") + 
    ggtitle(protein_name)
}

# GO enrichment and process results
protein_functions <- function(protein_list, protein_names) {
  # Convert proteins to ENTREZ IDs
  protein_conversion <- clusterProfiler::bitr(
    unique(protein_list),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(
    gene          = protein_conversion$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",  # Biological Process
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  # Process results to extract the top 20 most common processes
  go_enrichment@result %>%
    separate_rows(geneID, sep = "/") |>  
    left_join(protein_conversion, by = c("geneID" = "ENTREZID")) |>  
    group_by(SYMBOL) |>  
    slice_min(qvalue, with_ties = FALSE) 
  
}




# UMAP function

perform_umap <- function(data, protein_names, title, palette, n_neighbors = 15, min_dist = 0.1, metric = "euclidean") {
  
  # Filter the data for the selected protein names and disease class
  filtered_data_overlap <- data |> dplyr::select(Disease, all_of(protein_names))
  
  # Extract numeric data (exclude Disease column)
  numeric_data <- filtered_data_overlap[, -1]
  disease_labels <- filtered_data_overlap$Disease  # Extract disease labels
  
  # Perform UMAP
  set.seed(502)
  umap_result <- umap(
    numeric_data,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    init = "random"
  )
  
  # Convert UMAP results to a data frame
  umap_data <- as.data.frame(umap_result)
  colnames(umap_data) <- c("Dim1", "Dim2")  
  umap_data$Class <- factor(disease_labels)  
  
  # Plot UMAP
  umap_plot <- ggplot(umap_data, aes(x = Dim1, y = Dim2, color = Class)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = title,
      x = "Dimension 1",
      y = "Dimension 2",
      color = "Disease Class"
    ) +
    scale_color_manual(values = palette) +
    theme(legend.position = "right")
  
  return(umap_plot)
}

# PCA
perform_pca <- function(data, protein_names, tittle, palette) {
  
  filtered_data_overlap <- data |> dplyr::select(Disease, all_of(protein_names))
  
  # Perform PCA
  pca_result <- prcomp(filtered_data_overlap[, -1], center = TRUE, scale. = TRUE)
  pca_data <- as.data.frame(pca_result$x[, 1:2])  # Keep only the first two components
  colnames(pca_data) <- c("PC1", "PC2")
  
  # Add disease class to PCA data
  pca_data$disease <- factor(wide_data_auto$Disease)
  
  # Plot PCA
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = disease)) +
    geom_point(size = 3) +
    scale_color_manual(values = palette) +
    labs(title = tittle, x = "PC1", y = "PC2") +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(pca_plot)
}