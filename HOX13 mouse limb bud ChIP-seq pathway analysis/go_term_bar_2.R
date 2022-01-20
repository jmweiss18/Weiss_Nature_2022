go_term_bar <- function(gsea_results, n_terms = 10, metric = "pval", fill = T, change_labels = T, pval_fill = T) {
  # metric must be either "NES" or "pval"
  # set fill = T to fill with gradient colours, set fill = F to fill with grey only
  # set change_labels = T to change pathway names to lowercase and remove "GO:" from the pathway names, set change_labels = F to keep them as original
  # set pval_fill = T to color the bars based on p value
  
  require(tidyverse)
  require(pals)
  
  if (change_labels == T) {
    gsea_results$pathway <- gsea_results$pathway %>%
      gsub(pattern = "GO_", replacement = "", x = ., ignore.case = F) %>%
      gsub(pattern = "_", replacement = " ", x = .) %>%
      gsub(pattern = "plus", replacement = "+", x = .) %>%
      tolower()
  }
  
  if (fill == T) {
    cols_pval <- viridis(n_terms)
    cols_NES <- magma(n_terms)
  } else if (fill == F) {
    cols_pval <- "grey"
    cols_NES <- "grey" 
  } 
  
  if (metric == "pval") {
    
    data.bar <- gsea_results %>% 
      filter(NES > 0) %>% 
      arrange(padj) %>%
      slice(., 1:n_terms, .preserve = TRUE)
    
    # bar <- ggplot(data.bar, aes(x = log10pval, y = reorder(pathway, log10pval))) + 
    #   geom_bar(stat = "identity", fill = cols_pval) + 
    #   theme_minimal() + 
    #   labs(x = "-log10 p-value",
    #        y = "GO term") +
    #   theme(axis.title.x = element_text(size = 13, color = "black"),
    #         axis.title.y = element_blank(),
    #         axis.text.y = element_text(size = 13, color = "black"),
    #         axis.text.x = element_text(size = 12, color = "black"))
    
    if (pval_fill) {
      
      # updated to fix issue with ordering of Y axis
      data.bar$pathway <- factor(data.bar$pathway, levels = rev(data.bar$pathway))
      bar <- ggplot(data.bar, aes(x = padj, y = pathway, fill = padj)) +
        geom_bar(stat = "identity") + 
        scale_fill_gradientn(colours = ocean.matter(n = 1000) %>% rev()) +
        theme_minimal() + 
        labs(x = "padj",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar) 
      
    } else {
      
      data.bar$pathway <- factor(data.bar$pathway, levels = rev(data.bar$pathway))
      bar <- ggplot(data.bar, aes(x = padj, y = pathway)) +
        geom_bar(stat = "identity", fill = cols_pval) + 
        theme_minimal() + 
        labs(x = "padj",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar) 
      
    }
    
  } else if (metric == "NES") {
    
    data.bar <- gsea_results %>% 
      arrange(-NES) %>% 
      slice(., 1:n_terms, .preserve = TRUE)
    
    # bar <- ggplot(data.bar, aes(x = NES, y = reorder(pathway, NES))) + 
    #   geom_bar(stat = "identity", fill = cols_NES) + 
    #   theme_minimal() + 
    #   labs(x = "normalized enrichment score",
    #        y = "GO term") +
    #   theme(axis.title.x = element_text(size = 13, color = "black"),
    #         axis.title.y = element_blank(),
    #         axis.text.y = element_text(size = 13, color = "black"),
    #         axis.text.x = element_text(size = 12, color = "black"))
    
    if (pval_fill) {
      
      data.bar$pathway <- factor(data.bar$pathway, levels = rev(data.bar$pathway))
      bar <- ggplot(data.bar, aes(x = log2(NES), y = pathway, fill = padj)) +
        geom_bar(stat = "identity") + 
        scale_fill_gradientn(colours = ocean.matter(n = 1000)) +
        theme_minimal() + 
        labs(x = "log2(normalized enrichment score)",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar) 
      
    } else {
      
      data.bar$pathway <- factor(data.bar$pathway, levels = rev(data.bar$pathway))
      bar <- ggplot(data.bar, aes(x = NES, y = pathway)) + 
        geom_bar(stat = "identity", fill = cols_NES) + 
        theme_minimal() + 
        labs(x = "normalized enrichment score",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar) 
    }
    
  } else {
    stop("metric must be NES or pval")
  }
}
