top_bottom_goterms_plot <- function(gsea_results, nterms = 20, metric = "NES", split = T) {
  
  if (split == T) {
  
  if (metric == "pval") {
    gsea_results <- gsea_results %>% filter(NES > 0)
    gsea_results <- gsea_results %>% arrange(desc(log10pval))
    data.bar <- data.frame(pathway = gsea_results$pathway[1:nterms], log10pval = gsea_results$log10pval[1:nterms]) 
    
    bar <- ggplot(data.bar, aes(x = log10pval, y = reorder(pathway, log10pval))) + 
      geom_bar(stat = "identity") + 
      theme_minimal() + 
      labs(x = "-log10 p-value",
           y = "GO term") +
      theme(axis.title.x = element_text(size = 13, color = "black"),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 13, color = "black"),
            axis.text.x = element_text(size = 12, color = "black"))
    return(bar) 
    
  } else if (metric == "NES") {
    
    gsea_results <- gsea_results %>% arrange(desc(NES))
    
    data.bar <- rbind(head(gsea_results, nterms), tail(gsea_results, nterms)) %>% 
      dplyr::select(pathway, NES)
    
    bar <- ggplot(data.bar, aes(x = NES, y = reorder(pathway, NES))) + 
      geom_bar(stat = "identity", fill = brewer.rdbu(nterms*2)) + 
      theme_minimal() + 
      labs(x = "normalized enrichment score",
           y = "GO term") +
      theme(axis.title.x = element_text(size = 13, color = "black"),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.text.x = element_text(size = 12, color = "black")) 
    return(bar) 
    
  } else {
    stop("metric must be NES or pval")
  }
    
  } else if (split == F) {
    
    if (metric == "pval") {
      gsea_results <- gsea_results %>% filter(NES > 0)
      gsea_results <- gsea_results %>% arrange(desc(log10pval))
      data.bar <- data.frame(pathway = gsea_results$pathway[1:nterms], log10pval = gsea_results$log10pval[1:nterms]) 
      
      bar <- ggplot(data.bar, aes(x = log10pval, y = reorder(pathway, log10pval))) + 
        geom_bar(stat = "identity") + 
        theme_minimal() + 
        labs(x = "-log10 p-value",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar) 
      
    } else if (metric == "NES") {
      
      gsea_results <- gsea_results %>% arrange(desc(NES))
      
      data.bar <- head(gsea_results, nterms) %>% 
        dplyr::select(pathway, NES)
      
      bar <- ggplot(data.bar, aes(x = NES, y = reorder(pathway, NES))) + 
        geom_bar(stat = "identity", fill = (brewer.purd(nterms)) %>% rev()) + 
        theme_minimal() + 
        labs(x = "normalized enrichment score",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 11, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar) 
      
    } else {
      stop("metric must be NES or pval")
    }
    
  } else {
    stop("split must be T or F")
  }
}