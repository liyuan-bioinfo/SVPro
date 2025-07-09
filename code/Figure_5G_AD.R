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
dep_df = obj_list$dep_ct4_df_FC1_2
ct_order = obj_list$ct4_order

figure_no = "Figure_5G_"

# cluster enrichment analysis using enrich and sig. enrich
{
  table(dep_df$region)
  # Amy  CC  Hi  Th 
  # 52 174  16 160
  
  dep_df$pid = gsub(dep_df$pid, pattern=";.*", replacement = "")
  tran_df = bitr(dep_df$pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  plot_df = dep_df %>% merge(tran_df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)#482
  plot_df$region = factor(plot_df$region, levels=ct_order)
  plot_df = plot_df %>% dplyr::arrange(region)
  

  formula_res <- compareCluster(data = plot_df,ENTREZID~region,fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,minGSSize=10,maxGSSize=2000,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                
  )
  formula_res_cutoff = formula_res
  formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<=0.05,]
  
  
  selected_items_id = c("GO:1902284", "GO:0035066", "GO:0016239", "GO:0010498",
                        "GO:0043405", "GO:1902959", "GO:0150104", "GO:0031646")
  
  temp_df = formula_res_cutoff@compareClusterResult
  temp_df = temp_df %>% dplyr::filter(ID %in% selected_items_id)

  formula_res_cutoff@compareClusterResult = temp_df
  formula_res_cutoff@compareClusterResult$`log10P` = -log10(formula_res_cutoff@compareClusterResult$`p.adjust`)
  
  p2=dotplot(formula_res_cutoff, label_format=50,showCategory=2,font.size=14,color="log10P",size="count") + 
    theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
    scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30))
  
  # save dot plot
  pdf(file=paste0(output_dir,"\\",figure_no,"_AD_GOBP_compare.pdf"),width = 7,height = 6)
  print(p2)
  dev.off()   
  
}
