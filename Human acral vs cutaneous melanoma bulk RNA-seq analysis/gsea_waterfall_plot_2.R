gsea_waterfall_plot2 <- function(gsea_results, pathways_to_highlight_1, pathways_to_highlight_2 = NULL, labels_for_plot = NULL, n_terms = 250) {
  # gsea_results = your GSEA results table, should have 1 column called "pathway" and another called "NES" (additional columns are fine too)
  # pathways_label = a character vector with the exact names of the pathways you want highlighted on your waterfall plot in red circles
  # labels_for_plot (optional) = a character vector with the exact names of the pathways you want labelled as text on your plot
  
  # this will only plot the **top** pathways since with chip seq the enrichment scores are all positive
  require(tidyverse)
  
  # find the top 250 pathways by NES
  gsea.up <- gsea_results %>% 
    arrange(desc(NES)) %>%   
    head(., n_terms) %>% 
    dplyr::select(pathway, NES) %>% 
    rownames_to_column(var = "index")
  # 
  # # find the bottom 250 pathways by NES
  # gsea.down <- gsea_results %>% 
  #   arrange(desc(NES)) %>% 
  #   tail(., n_terms) %>% 
  #   dplyr::select(pathway, NES) %>% 
  #   rownames_to_column(var = "index") 
  
  # organize data frame for plotting
  gsea.up$index = 1:n_terms
  # gsea.down$index = n_terms:1
  # gsea.waterfall <- rbind(gsea.up, gsea.down)
  gsea.waterfall <- gsea.up
  label_pathways_1 <- intersect(pathways_to_highlight_1$pathway, gsea.waterfall$pathway)
  if (length(label_pathways_1) < 1) {
    stop("None of your pathways to label are within the top or bottom terms")
  }
  label_pathways_1 <- data.frame(pathway = label_pathways_1, label_pathway = "Limb Pathways")
 
  label_pathways_2 <- intersect(pathways_to_highlight_2$pathway, gsea.waterfall$pathway)
  label_pathways_2 <- data.frame(pathway = label_pathways_2, label_pathway = "insulin/IGF Pathways")
  
  labelled_pathways <- rbind(label_pathways_1, label_pathways_2)
  
  unlabelled_pathways <- data.frame(pathway = gsea.waterfall$pathway[!gsea.waterfall$pathway %in% labelled_pathways$pathway],
                                    label_pathway = "Other")
  pathways_all <- rbind(labelled_pathways, unlabelled_pathways)
  gsea.waterfall <- merge(gsea.waterfall,
                          pathways_all,
                          by = "pathway") 
  
  # create plot
  plot <- ggplot(gsea.waterfall, aes(x = index, y = log2(NES))) + 
    geom_point(aes(color = label_pathway), size = 8) +
    # scale_color_manual(values = c("#E01111", "#B2B2B220")) +
    scale_color_manual(values = c("#E01111", "navyblue", "#B2B2B220")) +
    geom_hline(yintercept = 0, size = 1) +
    theme_minimal() +
    # scale_y_continuous(limits = c(-3.1, 3.1), breaks = seq(-3,3,1)) + # change this to change the y axis limits
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 32, face = "bold"),
          axis.text.y = element_text(size = 24, color = "black"),
          panel.grid.major.x = element_blank()) +
    xlim(0, n_terms+1) +
    ylab("log2(normalized enrichment score)")
  
  if (is.null(labels_for_plot)) {
    return(plot)
  } else {
    # optional: add text labels to the plot highlighting certain pathways
    require(ggrepel)
    plot <- plot + geom_text_repel(data = subset(gsea.waterfall, pathway %in% labels_for_plot),
                                   aes(label = pathway),
                                   box.padding = unit(1, "lines"),
                                   point.padding = unit(1, "lines"))
    return(plot)
  }
}
