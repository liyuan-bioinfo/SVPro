library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readxl)

rm(list = ls())

# set basic info.
project_dir = "D:\\01_科研项目\\Project3_zhengzhendong\\submit_version"
source(paste(project_dir,"code\\tools.R",sep="\\"))
input_dir = paste0(project_dir,"\\", "input")
output_dir = paste0(project_dir,"\\", "output")

input_path = paste0(input_dir,"\\", "SVPro_resource.xlsx")

figure_no = "Figure_4H"

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
  df_v1_CC_GFAP = read_xlsx(input_path, sheet = "WT_CC-GFAP")
  df_v2_CC_NeuN = read_xlsx(input_path, sheet = "WT_CC-NeuN")
  df_v3_CC_NF200 = read_xlsx(input_path, sheet = "WT_CC-NF200")
  df_v4_Hi_GFAP = read_xlsx(input_path, sheet = "WT_Hi-GFAP")
  df_v5_Hi_NeuN = read_xlsx(input_path, sheet = "WT_Hi-NeuN")
  df_v6_Hi_NF200 = read_xlsx(input_path, sheet = "WT_Hi-NF200")
  
  # process
  df_v1_CC_GFAP_vol = process_data_zone_func(dat1 = df_v1_CC_GFAP, file_name = "v1_CC_GFAP")
  df_v2_CC_NeuN_vol = process_data_zone_func(dat1 = df_v2_CC_NeuN, file_name = "v2_CC_NeuN")
  df_v3_CC_NF200_vol = process_data_zone_func(dat1 = df_v3_CC_NF200, file_name = "v3_CC_NF200")
  df_v4_Hi_GFAP_vol = process_data_zone_func(dat1 = df_v4_Hi_GFAP, file_name = "v4_Hi_GFAP")
  df_v5_Hi_NeuN_vol = process_data_zone_func(dat1 = df_v5_Hi_NeuN, file_name = "v5_Hi_NeuN")
  df_v6_Hi_NF200_vol = process_data_zone_func(dat1 = df_v6_Hi_NF200, file_name = "v6_Hi_NF200")
  
  df_merge_six = rbind(df_v1_CC_GFAP_vol, df_v2_CC_NeuN_vol , df_v3_CC_NF200_vol, df_v4_Hi_GFAP_vol, df_v5_Hi_NeuN_vol,df_v6_Hi_NF200_vol)
  
  # plot
  pdf(paste0(output_dir,"\\",figure_no,"_Mutli_volcano_marker_varied_cell_type.pdf"))
  Multi_volcano_vary_marker_func(df = df_merge_six, compare_num = 6, manual_gene = TRUE, BG_length_equal = TRUE, size_num =3)
  dev.off()

}
