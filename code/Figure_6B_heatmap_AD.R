library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readxl)

rm(list = ls())

# set basic info.
project_dir = "SVPro" # modify this to the work dir

input_dir = paste0(project_dir,"\\", "input")
output_dir = paste0(project_dir,"\\", "output")

input_meta_path = paste0(input_dir,"\\", "AD_data_20250226.xlsx")

obj_list = readRDS(file=paste0(input_dir,"\\","AD_obj_list.rds"))
meta_df = obj_list$meta_df        
data_df = obj_list$data_df

figure_no = "Figure_6B"

# Heat map of selected 30 proteins
{
  
  selected_df = readxl::read_xlsx(input_meta_path, sheet="selected") %>% as.data.frame() %>% 
    dplyr::select(pid, gene)
  

  plot_df = na.omit(data_df[selected_df$pid,])
  row.names(plot_df) = selected_df$gene
  
  scale_plot_df = t(apply(plot_df, 1, scale)) %>% as.data.frame()
  names(scale_plot_df) =  names(plot_df)
  
  range_all <- range(c(-2, 2))    
  my_palette <- colorRampPalette(c("navy", "grey","white","#ffb74d","red"))(n=100)
  breaks = seq(range_all[1], range_all[2], length.out = 101)
  p1 = pheatmap::pheatmap(scale_plot_df,scale="none", color=my_palette,breaks=breaks,cellwidth = 13,cellheight = 10,
                          cluster_cols=F,cluster_rows=F,silence=T,show_rownames=T,fontsize=10)  
  pdf(file=paste0(output_dir, "\\", figure_no,"_heatmap_enrich_selected_30_",Sys.Date(),".pdf"),width=7,height = 7)
  print(p1)
  dev.off()
  
}

