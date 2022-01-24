library (dplyr)
library (ggplot2)
library (tidyr)
library (stringr)
library (ggrepel)
library (gplots)


setwd ("D:\\white_IP_for_submission\\IP3\\result")

allProtQuant <- read.csv ("2020_11_02_protein_quantitation_exp3.csv", stringsAsFactors = FALSE) %>%
                mutate (interactor = adj.pvalue_OEIPtoOEC < 0.05 & log2FC_OEIPtoOEC > 1) 

ggplot(allProtQuant , aes (log2FC_OEIPtoOEC, -log10(adj.pvalue_OEIPtoOEC), label = ifelse (interactor, SYMBOL, ""), colour = interactor)) + 
       geom_point() + 
       geom_text_repel (size = 2, show.legend = FALSE) + 
       geom_hline(yintercept= -log10 (0.05), linetype="dashed") + 
       geom_vline(xintercept= 1, linetype="dashed") + 
       xlab ("Log2 fold change of IP to control") + 
       ylab ("-Log10 of adjusted p-value") + theme_classic() + 
       xlim (-4.6, 4.6)


ggsave ("2021_10_10_volcano_plot.pdf", width = 10, height = 7)  

intQuant <- allProtQuant %>% filter (adj.pvalue_OEIPtoOEC < 0.05, log2FC_OEIPtoOEC > 1)

pdf ("2021_10_10_heatmap_plot.pdf", width = 8, height = 8)
heatmap.2 (as.matrix (intQuant[6:16]), trace = "none", labRow = xx$SYMBOL, labCol = c(rep ("IP", 6), rep ("control", 5)))
dev.off()           
