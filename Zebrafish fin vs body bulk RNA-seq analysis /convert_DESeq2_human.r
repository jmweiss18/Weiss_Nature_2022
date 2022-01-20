convert_DESeq2_human <- function(results_deseq) {
  # convert output of DESeq2 package to human genes for pathway analysis
  
  fish.human.convert <- readxl::read_excel("/Users/weiss/Documents/R/Learning RNA-seq/ZMEL_ANC_vs_EVO/Zebrafish_to_Human_new_R_XYZ.xlsx")
  
  results_deseq <- data.frame(results_deseq) %>% rownames_to_column(var = "fish_gene")
  
  genelist <- results_deseq
  merge <- merge(x = genelist,
                 y = fish.human.convert,
                 by.x = "Ensembl",
                 by.y = "Ensembl") 
  
  h.genelist.c <- merge$Human_Symbol_Unique %>% as.character()
  
  if (unique(h.genelist.c[duplicated(h.genelist.c)]) == "NA") {
    message("No duplicated human genes remain, returning gene list")
    genelist_final <- merge %>% dplyr::select(Human_Symbol_Unique, baseMean, log2FoldChange, pvalue, padj, everything())
    genelist_final <- genelist_final[!genelist_final$Human_Symbol_Unique == "NA", ] # remove genes with no human ortholog
    return(genelist_final)
  } else {
    stop("Some duplicated human genes remain")
  }
}
  
  
#   h.genelist.c <- merge$Human_Symbol %>% as.character()
#   h.genelist <- list(h.genelist.c)
#   if (length(unique(duplicated(h.genelist.c))) == 1) {
#     message("No duplicated genes, conversion complete")
#     return(h.genelist)
#   } else {
#     message("Some genes are duplicated, removing...")
#     gene.dups <- as.character(unique(merge$Human_Symbol_Unique[duplicated(merge$Human_Symbol_Unique)]))
#     dups <- vector("list", length(gene.dups))
#     for (ii in 1:length(gene.dups)) {
#       genes <- merge[merge$Human_Symbol == gene.dups[ii], ] %>% 
#         arrange(desc(avg_logFC))
#       dups[[ii]] <- genes[1,]
#     }
#     dups <- data.table::rbindlist(dups)
#     nondups <- merge[grep(paste(pattern = gene.dups,
#                                 collapse = "|"),
#                           x = merge$Human_Symbol,
#                           invert = TRUE), ]
#     test1 <- unique(duplicated(nondups$Human_Symbol)) 
#     test2 <- unique(duplicated(dups$fish_gene)) 
#     if (test1 == FALSE & test2 == FALSE) {
#       h.genelist <- rbind(dups, nondups) 
#       test3 <- unique(duplicated(h.genelist$fish_gene))
#       test4 <- unique(duplicated(h.genelist$Human_Symbol))
#       if (test3 == FALSE & test4 == FALSE) {
#         message("Duplicated genes removed, conversion complete")
#         return(h.genelist)
#       } else {
#         stop("Some duplicated genes weren't removed properly") }
#     } else {
#       stop("Some duplicated genes weren't removed properly")
#     }
#   }
# }
