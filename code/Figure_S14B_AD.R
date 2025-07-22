library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggrepel)
library(ggpubr)

rm(list = ls())


# set basic info.
project_dir = "SVPro" # modify this to the work dir
input_dir = paste0(project_dir,"\\", "input")
output_dir = paste0(project_dir,"\\", "output")

obj_list = readRDS(paste0(input_dir, "\\AD_obj_list.rds"))
meta_df = obj_list$meta_df        
data_df = obj_list$data_df

enrich_df = obj_list$enrich_df

figure_no = "Figure_S14B_"

# enrichment analysis of overlap protein lists
{

  enrich_df = na.omit(enrich_df)
  enrich_df$pid = gsub(row.names(enrich_df), pattern=";.*", replacement = "")
  
  
  
  tran_df = bitr(enrich_df$pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  formula_res <- clusterProfiler::enrichGO(gene = tran_df$ENTREZID,
                                           OrgDb = org.Mm.eg.db,minGSSize=10,maxGSSize=2000,
                                           ont = "BP", pAdjustMethod = "BH",#universe = bg_pid_df$ENTREZID,
                                           pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                           
  )
  formula_res_cutoff = formula_res
  formula_res_cutoff@result = formula_res@result[formula_res@result$p.adjust<=0.05,]
  

  # selected_items = c("vesicle-mediated transport", "neuron projection development",
  #                    "neuron differentiation", "cell morphogenesis involved in differentiation",
  #                    "cell morphogenesis involved in differentiation","chemical synaptic transmission",
  #                    "cytoskeleton organization", "cell-cell signaling", "synaptic vesicle cycle",
  #                    "export from cell", "protein localization to membrane"
  #                    )
  
  selected_ID = c("GO:0016192", "GO:0031175", "GO:0030182", "GO:0000904", "GO:0007268",
                  "GO:0007010", "GO:0007267", "GO:0099504", "GO:0140352", "GO:0072657")
  selected_GO = formula_res_cutoff@result
  formula_res_cutoff@result = selected_GO %>% dplyr::filter(ID %in% selected_ID)
  
  p2=formula_res_cutoff %>% barplot(label_format=50,showCategory=10,font.size=14,fill="p.adjust",size="count") +
    theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
    scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30))
  p2$data$p.adjust <- -log10(p2$data$p.adjust)
  p2 <- p2 + aes(fill = p.adjust) +
    scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30), name = "-log10(p.adjust)")
  
  # save dot plot
  pdf(file=paste0(output_dir,"\\",figure_no,"_AD_GOBP_barplot_",Sys.Date(),".pdf"),width = 8,height = 5)
  print(p2)
  dev.off() 
  
  
}
