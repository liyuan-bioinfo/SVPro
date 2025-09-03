library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readxl)

rm(list = ls())

project_dir = "/home/path/to/SVPro" # set your own project directory
source(paste(project_dir,"code","tools.R",sep="/"))
input_dir = paste0(project_dir,"/", "input")
output_dir = paste0(project_dir,"/", "output")

input_path = paste0(input_dir,"/", "SVPro_resource.xlsx")

figure_no = "Figure_5E"

# main
{
  
  target_gene = c("Serpina3", "C1qa", "C1qc", "Eno1", "Cfh", "Pgk1", "Bcan", "Apod", "Apoe",
                  "Clu", "Ncam1", "Vtn", "Actb", "Gsn", "Olfml3", "Pcsk1n", "Cst3", "App",
                  "PSEN1", "Abeta42", "Mapt")
  
  nature_color_top7 = c("#F798B6" , "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF")
  
  #Marker proteins
  NeuN_pro = c("Calm1", "Gria2", "Map1a", "Rps16", "Hnrnpa0",
               "Hnrnpa2b1", "Hnrnpdl", "Hnrnpl", "Hnrnpa3", "Hnrnpf", "Hnrnph2",
               "Rpl26", "Rps11", "Rbfox3")
  
  
  NF200_pro = c("Bsn", "Homer1", "Slc17a7", "Nefm", "Sv2b", "Vamp2", "Syp",
                "Ppp1r9b", "Nefl", "Sptbn1", "Syn1/Syn2", "Ina", "Cltc", "Camk2A/B",
                "Syt1", "Tuba4a (TUBA1)", "Mapt", "Dlg4", "Nefh", "Vamp3", "Syn1", "Syn2", "Camk2a", "Camk2b", "Tuba4a")
  
  
  GFAP_pro = c("Slc1a2 (Glt1)", "Aqp4", "Gfap", "Aldh1l1", "Acsbg1", 
               "Glul", "Acsbg1", "Slc6a11", "Atp1b2", "Atp1a2","Apoe",
               "Gja1", "Slc25a18", "Gstm1", "Slc1a2")
  
  
  # load data
  df_v1_Amy = read_xlsx(input_path, sheet = paste0("Fig 5E_","AD_Amy"))
  df_v2_CC = read_xlsx(input_path, sheet = paste0("Fig 5E_","AD_CC"))
  df_v3_Hi = read_xlsx(input_path, sheet = paste0("Fig 5E_","AD_Hi"))
  df_v4_Th = read_xlsx(input_path, sheet = paste0("Fig 5E_","AD_Th"))
  
  # process data
  df_v1_Amy_vol = process_data_zone_func(dat1 = df_v1_Amy, file_name = "v1_Amy")
  df_v2_CC_vol = process_data_zone_func(dat1 = df_v2_CC, file_name = "v2_CC")
  df_v3_Hi_vol = process_data_zone_func(dat1 = df_v3_Hi, file_name = "v3_Hi")
  df_v4_Th_vol = process_data_zone_func(dat1 = df_v4_Th, file_name = "v4_Th")
  
  df_merge = rbind(df_v1_Amy_vol, df_v2_CC_vol, df_v3_Hi_vol, df_v4_Th_vol)
  
  # plot and save
  pdf(paste0("output","/",figure_no,"_Mutli_volcano_munal_add_gene.pdf"))
  Multi_volcano_label_gene_manual_func(df = df_merge, compare_num = 4, manual_gene = TRUE, BG_length_equal = TRUE, size_num =3)
  dev.off()
}

