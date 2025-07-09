library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readxl)

rm(list = ls())


# Basic pre-process
{
  project_dir = "D:\\01_科研项目\\Project3_zhengzhendong\\submit_version"
  
  input_dir = paste0(project_dir,"\\", "input")
  output_dir = paste0(project_dir,"\\", "output")
  
  input_path = paste0(input_dir,"\\", "SVPro_resource.xlsx")
  input_meta_path = paste0(input_dir,"\\", "AD_data_20250226.xlsx")


  # load meta files and order meta
  meta_df = readxl::read_xlsx(input_meta_path, sheet="meta") %>% as.data.frame() # 12 * 3
  
  ct4_order = c("AD_Amy", "AD_CC", "AD_Hi", "AD_Th")
  meta_df$Region = factor(meta_df$Region, levels = ct4_order)
  meta_df = meta_df %>% arrange(Region)
  
  # load quantified data
  data_df = readxl::read_xlsx(input_path, sheet="Normalized data_AD_Exp") %>% as.data.frame() # 2658 * 15
  dim(data_df)
  row.names(data_df) = data_df$`PG.ProteinGroups`
  
  # load sig. enriched proteins for each region
  union_pids = c()
  protein_lists <- list()
  for (ct in ct4_order){
    # print(ct)
    temp_df = readxl::read_xlsx(input_path, sheet=ct) %>% as.data.frame() # 2658
    temp_df = temp_df[which(temp_df$`-Log Student's T-test p-value NAB228_Ab.omitted` >= -log10(0.05)  &
                              temp_df$`Student's T-test Difference NAB228_Ab.omitted` >= log2(2)),]
    temp_pid = temp_df$PG.ProteinGroups
    union_pids = unique(c(union_pids, temp_pid))
    
    protein_lists[[ct]] <- temp_pid
  }
  

  all_pids <- unique(unlist(protein_lists))
  
  enrich_df <- data.frame(ProteinID = all_pids)
  
  for (ct in ct4_order) {
    current_pids <- protein_lists[[ct]]
    
    enrich_df[[ct]] <- ifelse(enrich_df$ProteinID %in% current_pids, enrich_df$ProteinID, NA)
  }
  row.names(enrich_df) = enrich_df$ProteinID
  enrich_df$ProteinID = NULL
  
  # retain annotation files
  anno_df = data_df[,c("PG.ProteinGroups","PG.Genes")]
  names(anno_df) = c("pid", "genes")
  anno_df$gene = gsub(anno_df$genes, pattern=";.*", replacement = "")
  
  # Select samples and re-name sampleID
  data_df = data_df[union_pids,meta_df$RawSampleID] %>% unique() # 2195
  names(data_df) = meta_df$SampleID # 2195 * 12
  
  # replace NaN with 0
  data_df = data_df %>%
    mutate(across(everything(), ~ as.numeric(.) %>% replace_na(0)))
  
  # obtain the mean value of each cell type
  median_data_df = data.frame(row.names = row.names(data_df))
  
  for (ct in ct4_order){
    temp_id = meta_df[which(ct == meta_df$Region),"SampleID"]
    temp_df = data_df[,temp_id]
    temp_median = apply(log2(temp_df+1),1,median)
    median_data_df = cbind(median_data_df, temp_median)
  }
  names(median_data_df) = ct4_order  
  
  
  # save RDATA
  obj_list = list()
  obj_list$data_df = log2(data_df+1)
  obj_list$median_data_df = median_data_df
  obj_list$enrich_df = enrich_df
  obj_list$meta_df = meta_df
  obj_list$anno_df = anno_df
  obj_list$ct4_order = ct4_order
  obj_list$ct4_color = c("#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF")
  
  saveRDS(obj_list, file=paste0(input_dir,"\\AD_obj_list.rds"))
  
}

# DEP analysis
{
  rm(list=ls())
  project_dir = "D:\\01_科研项目\\Project3_zhengzhendong\\submit_version"
  input_dir = paste0(project_dir,"\\", "input")
  output_dir = paste0(project_dir,"\\", "output")
  
  obj_list = readRDS(file = paste0(input_dir,"\\AD_obj_list.rds"))
  
  meta_df = obj_list$meta_df        
  data_df = obj_list$data_df
  median_data_df = obj_list$median_data_df
  enrich_df = obj_list$enrich_df
  protein_num = dim(data_df)[1]
  ct4_order = obj_list$ct4_order
  ct_num = 4
  enrich_fc = 1.2
  

  # to find enriched cell-types with mean abundance
  dep_df = data.frame()
  for(i in 1:ct_num){
    temp_df = median_data_df[,i] - median_data_df[,1:ct_num]
    temp_pid = names(which(rowSums(temp_df>log2(1.2)) == (ct_num-1)))
    temp_df2 = data.frame(pid=temp_pid)
    temp_df2$region = names(median_data_df)[i]
    temp_df2$log2FC = apply(temp_df[temp_pid,-i],1,median)        
    dep_df = rbind(dep_df,temp_df2)            
  }
  
  ## pvalue, with each sample
  row.names(dep_df) = dep_df$pid
  enrich_data_df = na.omit(data_df[dep_df$pid,]) # 1281
  protein_num = dim(enrich_data_df)[1] #update for enriched proteins
  
  Pvalue = c()
  for (i in 1:protein_num) {
    
    pid_df = enrich_data_df[i,] %>% t() %>% as.data.frame()
    names(pid_df) = "pid"
    pid_df$SampleId = row.names(pid_df)
    pid_df$region = meta_df$Region
    
    enriched_ct = dep_df[i, "region"]
    
    target_data = pid_df[pid_df$region == enriched_ct, "pid"]
    
    other_regions = unique(pid_df$region[pid_df$region != enriched_ct])
    
    p_values = c()
    
    for (region in other_regions) {
      
      other_data = pid_df[pid_df$region == region, "pid"]
      
      t_test_result = t.test(target_data, other_data)
      
      p_values = c(p_values, t_test_result$p.value)
    }
    
    # combine p-value（using Fisher's method
    combined_pvalue = RecordTest::fisher.method(p_values)$"p.value"[[1]]
    
    Pvalue = c(Pvalue, combined_pvalue)
  }  
  
  
  dep_df$pvalue = Pvalue
  dep_df$fdr = p.adjust(Pvalue,method = "BH")                        
  dep_df$gene = obj_list$anno_df[dep_df$pid,"gene"]
  length(which(dep_df$pvalue < 0.05))# 447
  dep_df = dep_df %>% dplyr::filter(pvalue < 0.05)
  table(dep_df$region)
  
  dep_df = dep_df %>% arrange(region, log2FC)
  
  obj_list$dep_df = dep_df
  saveRDS(obj_list, file=paste0(input_dir,"\\AD_obj_list.rds"))
  
}
  

# obtain 402 enrich and sig. proteins
{
 
  rm(list=ls())
  
  project_dir = "D:\\01_科研项目\\Project3_zhengzhendong\\submit_version"
  
  input_dir = paste0(project_dir,"\\", "input")
  output_dir = paste0(project_dir,"\\", "output")
  
  obj_list = readRDS(file = paste0(input_dir,"\\AD_obj_list.rds"))
  meta_df = obj_list$meta_df        
  data_df = obj_list$data_df

  enrich_df = obj_list$enrich_df
  dep_df = obj_list$dep_df

  ct4_order = obj_list$ct4_order

  enrich_df = enrich_df[dep_df$pid,]
  dep_df = cbind(dep_df, enrich_df)
  
  # retain enriched for each region
  output_df = NA
  for (ct in ct4_order){
    temp_df = dep_df[dep_df$region == ct,]
    print(dim(temp_df))
    temp_df = temp_df[!is.na(temp_df[ct]),]
    print(dim(temp_df))
    if(length(output_df)==1){
      output_df = temp_df
    }else{
      output_df = rbind(output_df, temp_df)  
    }
    
  }
  dim(output_df) # 402 * 10
  # write.csv(output_df,file=paste0("analysis/impute_enrichment_celltype4_FC1.2_402_", Sys.Date(),".csv"))
  

  obj_list$dep_ct4_df_FC1_2 = output_df
  saveRDS(obj_list, file=paste0(input_dir,"\\AD_obj_list.rds"))
  
}

  

