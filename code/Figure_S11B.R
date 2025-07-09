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

region_df = read.csv(file=paste0(input_dir,"\\enrich_proteins_Hi.csv"),header = T,row.names = 1)
figure_no = "Figure_S11B"

# Hi [Hippocampus]
{
  
  region_df[!is.na(region_df)] <- 1
  region_df[is.na(region_df)] <- 0
  
  temp_sum = t(apply(region_df, 1, as.numeric))
  temp_sum = apply(temp_sum, 1, sum)
  region_df = region_df[names(which(temp_sum == 1)),]
  
  region_df$pids = row.names(region_df)
  
  temp_df = region_df %>% tidyr::gather(key="region",value="value",-pids)
  venn_df = temp_df %>% dplyr::filter(value == 1)
  table(venn_df$region)
  # GFAP_Hi  NeuN_Hi NF200_Hi 
  # 131      226       81
  
  venn_df$pid = gsub(venn_df$pids, pattern=";.*", replacement = "")
  row.names(venn_df) = venn_df$pids
  
  venn_df$gene = anno_df[row.names(venn_df),"genes"]
  
  tran_df = bitr(venn_df$pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  plot_df = venn_df %>% merge(tran_df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)

  formula_res <- compareCluster(data = plot_df,ENTREZID~region,fun="enrichGO",
                                OrgDb = org.Mm.eg.db,minGSSize=10,maxGSSize=2000,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                
  )
  
  formula_res_cutoff = formula_res
  formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<=0.05,]
  
  # save this file to select meaningful items
  # write.csv(x = formula_res_cutoff@compareClusterResult, file=paste0("xxxx.csv"))
  
  # load labeled items  
  
  labeled_df_NF200 = readxl::read_xlsx(path=paste0(input_dir,"\\Venn_protein438_Hi_GOBP_formula_res_2025-03-17.xlsx"),sheet="NF200") %>% 
    as.data.frame()%>% dplyr::filter(label == 1)
  labeled_df_NeuN = readxl::read_xlsx(path=paste0(input_dir,"\\Venn_protein438_Hi_GOBP_formula_res_2025-03-17.xlsx"),sheet="NeuN") %>% 
    as.data.frame()%>% dplyr::filter(label == 1)
  labeled_df_GFAP = readxl::read_xlsx(path=paste0(input_dir,"\\Venn_protein438_Hi_GOBP_formula_res_2025-03-17.xlsx"),sheet="GFAP") %>% 
    as.data.frame()%>% dplyr::filter(label == 1)
  
  
  selected_ID = c(labeled_df_NF200$ID, labeled_df_NeuN$ID, labeled_df_GFAP$ID)
  
  
  selected_GO = formula_res_cutoff@compareClusterResult
  formula_res_cutoff@compareClusterResult = selected_GO %>% dplyr::filter(ID %in% selected_ID)
  
  formula_res_cutoff@compareClusterResult$`log10Adj.P` = -log10(formula_res_cutoff@compareClusterResult$`p.adjust`)
  p1=dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=10,color="log10Adj.P",size="count") + 
    theme(panel.grid = element_blank(),axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))  +
    
    scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30), name = "-log10(p.adjust)")    
  
  
  # save dot plot
  pdf(file=paste0(output_dir,"\\",figure_no,"_Hi_GOBP.pdf",Sys.Date(),".pdf"),width = 6,height = 6)
  print(p1)
  dev.off() 
  
  
}  
