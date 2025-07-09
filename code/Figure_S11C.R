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

figure_no = "Figure_S11C"
total_enrich_df = read.csv(paste0(input_dir,"\\enrich_proteins_WT.csv"),header = T,row.names = 1)

# NeuN [Soma]
{
  
  enrich_df = total_enrich_df[,c(3:4)]# NeuN_CC NeuN_Hi
  head(enrich_df)
  
  enrich_df[!is.na(enrich_df)] <- 1
  enrich_df[is.na(enrich_df)] <- 0
  
  temp_sum = t(apply(enrich_df, 1, as.numeric))
  temp_sum = apply(temp_sum, 1, sum)
  
  enrich_df$pids = row.names(enrich_df)
  enrich_df$pid = gsub(enrich_df$pids, pattern=";.*", replacement = "")
  enrich_df$type = NA
  enrich_df[which(temp_sum == 1),"type"] = "specific"
  
  specific_df = enrich_df %>% dplyr::filter(type == "specific")
  specific_df$type = NULL
  specific_df = specific_df %>% tidyr::gather(key="group",value="value",-c(pids,pid)) %>% 
    dplyr::filter(value == 1) %>% dplyr::select(pids,pid, group)
  table(specific_df$group)
  # NeuN_CC NeuN_Hi 
  # 233     229
  
  tran_df = bitr(specific_df$pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)  
  plot_df = specific_df %>% merge(tran_df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)
  
  formula_res <- compareCluster(data = plot_df,ENTREZID~group,fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,minGSSize=10,maxGSSize=2000,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                
  )
  formula_res_cutoff = formula_res
  formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<=0.05,]
  
  # saveRDS(formula_res, file=paste0(save_path,"_formula_res_",Sys.Date(),".rds"))
  # write.csv(x = formula_res_cutoff@compareClusterResult, file=paste0("xxx.csv"))  

  # load labeled items
  labeled_df_CC = readxl::read_xlsx(path=paste0(input_dir,"\\WT_NeuN_Compare_GOBP__formula_res_2025-03-17.xlsx"),sheet="CC") %>% 
    as.data.frame()%>% dplyr::filter(label == 1)
  labeled_df_Hi = readxl::read_xlsx(path=paste0(input_dir,"\\WT_NeuN_Compare_GOBP__formula_res_2025-03-17.xlsx"),sheet="Hi") %>% 
    as.data.frame()%>% dplyr::filter(label == 1)
  
  
  selected_ID = c(labeled_df_CC$ID, labeled_df_Hi$ID)
  
  selected_GO = formula_res_cutoff@compareClusterResult
  formula_res_cutoff@compareClusterResult = selected_GO %>% dplyr::filter(ID %in% selected_ID)
  
  formula_res_cutoff@compareClusterResult$`log10Adj.P` = -log10(formula_res_cutoff@compareClusterResult$`p.adjust`)
  p1=dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=10,color="log10Adj.P",size="count") + 
    theme(panel.grid = element_blank(),axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    
    scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30), name = "-log10(p.adjust)")    
  
  
  # save dot plot
  pdf(file=paste0(output_dir,"\\",figure_no,"_WT_NeuN_GOBP.pdf"),width = 6,height = 4)
  print(p1)
  dev.off()     
} 
