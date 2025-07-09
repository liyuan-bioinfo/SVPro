library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggrepel)
library(ggpubr)

rm(list = ls())


# set basic info.
project_dir = "D:\\01_科研项目\\Project3_zhengzhendong\\submit_version"
input_dir = paste0(project_dir,"\\", "input")
output_dir = paste0(project_dir,"\\", "output")

# region_df = read.csv(file=paste0(input_dir,"\\enrich_proteins_WT_CC.csv"),header = T,row.names = 1, check.names = F)
obj_list = readRDS(paste0(input_dir, "\\AD_obj_list.rds"))
meta_df = obj_list$meta_df        
data_df = obj_list$data_df


figure_no = "Figure_6_boxplots_"


{
  
  
  anno_df = na.omit(obj_list$anno_df[row.names(data_df),])
  
  i_list = which(anno_df$gene %in% c("Apoe","Clu","Slc6a11","Gpc4"))
  
  plot_list = list()
  for (i in i_list){
    
    
    pid = row.names(data_df)[i]
    
    gene = obj_list$anno_df[pid,"gene"]
    
    subtitle = paste0(pid,"_",gene)
    
    temp_df = t(data_df[i,]) %>% as.data.frame()
    names(temp_df) = "abundance"
    
    temp_df$celltype = meta_df$Region
    
    comparisons <- combn(levels(temp_df$celltype), 2, simplify = FALSE)
    
    p1=temp_df %>% 
      ggplot(aes(x = celltype, y = abundance, color = celltype)) +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8
        
      ) +
      geom_jitter(
        width = 0.15,
        alpha = 1,   
        size = 2
      
      ) +
      
      ggpubr::stat_compare_means(
        comparisons = comparisons,
        method = "t.test",
        label = "p.signif",
        step.increase = 0.15,
        tip.length = 0.01    
      ) +
      labs(
        x = "",
        y = "Log2(Intensity)",
        subtitle = subtitle
      ) +
      theme_bw(base_size = 16) +
      scale_y_continuous(expand = c(0, 0), limits = c(12, 22))+
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none", 
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank()
      )
    plot_list[[pid]] = p1
  }
  
  # arrange and save plot
  p1 = ggpubr::ggarrange(plotlist = plot_list, ncol = 2,nrow = 2)
  pdf(file=paste0(output_dir,"\\",figure_no,"_t-test.pdf"),width = 8,height = 8)
  print(p1)
  dev.off()
  
  
  
}  
