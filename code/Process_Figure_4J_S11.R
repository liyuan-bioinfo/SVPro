library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readxl)

rm(list = ls())

# set basic info.
project_dir = "D:\\01_科研项目\\Project3_zhengzhendong\\submit_version"

input_dir = paste0(project_dir,"\\", "input")
output_dir = paste0(project_dir,"\\", "output")

input_path = paste0(input_dir,"\\", "SVPro_resource.xlsx")

{

  # Venn plot, obtain union proteins for three cell types
  {
    ct6_order = c("WT_CC-GFAP","WT_Hi-GFAP","WT_CC-NeuN","WT_Hi-NeuN","WT_CC-NF200","WT_Hi-NF200")
    union_pids <- c()
    protein_lists <- list()
    
    for (ct in ct6_order) {
      temp_df = readxl::read_xlsx(input_path, sheet=ct) %>% as.data.frame()
      temp_df = temp_df[which(temp_df[10] >= -log10(0.05)  & #-Log Student's T-test p-value 
                                temp_df[12] >= log2(2)),] #Student's T-test Difference CC-NeuN_CC-Ctr
      
      temp_pid <- temp_df$PG.ProteinGroups
    
      protein_lists[[ct]] <- temp_pid
      
      union_pids <- unique(c(union_pids, temp_pid))
    }
    
    all_pids <- unique(unlist(protein_lists))
    
    result_df <- data.frame(ProteinID = all_pids)
    
    for (ct in ct6_order) {
      current_pids <- protein_lists[[ct]]
      
      result_df[[ct]] <- ifelse(result_df$ProteinID %in% current_pids, result_df$ProteinID, NA)
    }
    
    # save file
    write.csv(result_df, file = paste0(input_dir,"\\enrich_proteins_WT.csv"), row.names = FALSE)  
    dim(result_df) # 1349
  }
  
  
  # Venn plot, obtain union proteins for two regions
  {
    region_order = c("CC","Hi")
    for (region in region_order){
      if(region == "CC"){
        ct6_order = c("WT_CC-GFAP","WT_CC-NeuN",  "WT_CC-NF200")    
      }else if(region == "Hi"){
        ct6_order = c("WT_Hi-GFAP",  "WT_Hi-NeuN",  "WT_Hi-NF200")    
      }
      
      
      
      union_pids <- c()
      protein_lists <- list()
      
      for (ct in ct6_order) {

        temp_df = readxl::read_xlsx(input_path, sheet=ct) %>% as.data.frame()
        temp_df = temp_df[which(temp_df[10] >= -log10(0.05)  &
                                  temp_df[12] >= log2(2)),]
        
        temp_pid <- temp_df$PG.ProteinGroups
        
        protein_lists[[ct]] <- temp_pid
        
        union_pids <- unique(c(union_pids, temp_pid))
      }
      
      all_pids <- unique(unlist(protein_lists))
      
      result_df <- data.frame(ProteinID = all_pids)
      
      for (ct in ct6_order) {

        current_pids <- protein_lists[[ct]]
        
        result_df[[ct]] <- ifelse(result_df$ProteinID %in% current_pids, result_df$ProteinID, NA)
      }
      
      # save file
      write.csv(result_df, file = paste0(input_dir,"\\enrich_proteins_",region,".csv"), row.names = FALSE)  
    }
  }
  
  
}

