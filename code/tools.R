library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readxl)

Multi_volcano_label_gene_manual_func = function(df = df_merge, compare_num = 4, manual_gene = TRUE, BG_length_equal = TRUE, size_num =3 )
{
  
  
  df <- df %>%
    mutate(
      Cluster = Type,
      Log2FC = Difference,
      Neg_log_P = P_value,
      Label = ifelse(Neg_log_P >= -log10(0.05) & Difference >= 1,
                     "1_sig_P_value<=0.05_Log2FC>=1",
                     "2_other_no_sig")
    ) %>%
    mutate(
      Label = ifelse(Gene %in% target_gene & Neg_log_P >= -log10(0.05) & Difference >= 1,
                     "0_sig_label",
                     Label),
      size = 1
    )
  
  
  
  
  manual_gene_index = which(df$Gene %in% target_gene)
  
  df$size[manual_gene_index] = size_num
  
  
  
  
  df <- df %>%
    dplyr::mutate(Cluster_num = as.numeric(factor(Cluster)) - 1)
  
  
  df$Gene_slim = sapply(df$Gene, function(x){unlist(strsplit(x, split = "_"))[1]})
  
  df_top10 = df[df$size==size_num,]
  df_not_top10 = df[df$size == 1,]
  
  df_not_top10_v2_nosig = df_not_top10 %>% dplyr::filter(Label == "2_other_no_sig")
  df_not_top10_v1_sig = df_not_top10 %>% dplyr::filter(Label == "1_sig_P_value<=0.05_Log2FC>=1")
  
  df_sig_num = as.data.frame(table(df$Cluster, df$Label))
  
  colnames(df_sig_num) = c("Cluster", "Label", "Counts")
  
  df_sig_num$Cluster = as.character(df_sig_num$Cluster)
  df_sig_num$Label = as.character(df_sig_num$Label)
  
  
  df_filtered <- df_sig_num %>%
    filter(Label != "2_other_no_sig")
  
  df_filtered$Label = NULL
  
  df_filtered = df_filtered %>%   
    dplyr::group_by(Cluster) %>%                   
    dplyr::summarise(Total_Counts = sum(Counts))   # Sum counts
  
  
  
  
  # Calculate Max and Min Log2FC for each cluster
  min_max_per_cluster <- df %>%
    group_by(Cluster) %>%
    summarise(
      Max_Log2FC = max(Log2FC),
      Min_Log2FC = min(Log2FC),
      .groups = 'drop'
    )
  
  
  
  
  
  
  
  #Check length of rectangle for background
  
  if(!BG_length_equal)
    dfbar<-data.frame(x=seq(0,compare_num - 1,by=1),
                      y=  min_max_per_cluster$Max_Log2FC+0.1)
  
  dfbar1<-data.frame(x=seq(0,compare_num - 1,by=1),
                     y=min_max_per_cluster$Min_Log2FC - 0.1)
  
  if(BG_length_equal)
  {
    max_value = max(min_max_per_cluster$Max_Log2FC) + 0.1
    min_value = min(min_max_per_cluster$Min_Log2FC) - 0.1
    
    dfbar<-data.frame(x=seq(0,compare_num - 1,by=1), y = max_value)
    dfbar1 <- data.frame(x=seq(0,compare_num - 1,by=1), y = min_value)
    
  }
  
  p2 = ggplot()+
    geom_col(data = dfbar,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = dfbar1,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = df_not_top10_v2_nosig,
                aes(x = Cluster_num, y = Log2FC, color = Label),
                size = 0.8,
                #width =0.4
                position = position_jitter(seed = 123))+
    
    geom_jitter(data = df_not_top10_v1_sig,
                aes(x = Cluster_num, y = Log2FC, color = Label),
                size = 0.8,
                #width =0.4
                position = position_jitter(seed = 123))+
    
    geom_jitter(data = df_top10,
                aes(x = Cluster_num, y = Log2FC, color = Label),
                #size = 1.6,
                size = size_num,
                #width =0.4
                position = position_jitter(seed = 123))
  
  #p2
  #########
  #
  #Add cluster-colored blocks as labels on the X-axis.
  #
  ##########
  dfcol<-data.frame(x=c(0:(compare_num - 1)),
                    y=0,
                    label=c(0:(compare_num - 1)),
                    label_name = sort(unique(df$Cluster, decreasing = FALSE)))
  
  
  mycol <- nature_color_top7[1:compare_num]
  
  
  p3 <- p2 + geom_tile(data = dfcol,
                       aes(x=x,y=y),
                       height=0.4,
                       color = "black",
                       fill = mycol,
                       alpha = 0.6,
                       show.legend = F)
  
  p3 <- p3 +
    scale_color_manual(name=NULL,
                       values = c("darkred","red","grey"))
  
  
  p3 <- p3+
    labs(x="Cluster",y="average logFC")+
    geom_text(data=dfcol,
              aes(x=x,y=y,label=label_name),
              size =3,
              color ="black")
  
  
  
  
  
  p3 <- p3+
    theme_minimal()+
    theme(
      axis.title = element_text(size = 13,
                                color = "black",
                                face = "bold"),
      axis.line.y = element_line(color = "black",
                                 size = 1.2),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 15)
    )
  
  #p3
  
  p4 <- p3+
    geom_text(
      data=df_top10,
      aes(x=Cluster_num,y=Log2FC,label=Gene_slim),
      size = 2.5, col = "black",
      position = position_jitter(seed = 123)
    )
  
  df_text = df_filtered
  df_text$x = seq(0.2,3.2,by=1)
  df_text$y = 12
  
  
  ########
  #
  #add text and plot the figure 
  #
  ########
  
  p5 <- p4 + geom_text(data = df_text, aes(x = x, y = y, label = paste("sig_num: ",Total_Counts, sep = "")), size = 3, fontface = "bold") 
  
  print(p5)
  
  
  
}



Multi_volcano_vary_marker_func = function(df = df_merge_six, compare_num = 6, manual_gene = TRUE, BG_length_equal = TRUE, size_num =3 )
{
  
  df <- df %>%
    mutate(
      Gene = sapply(Gene, function(x) unlist(strsplit(x, split = ";"))[1]),
      Cluster = Type,
      Log2FC = Difference,
      Neg_log_P = P_value,
      Label = ifelse(Neg_log_P >= -log10(0.05) & Difference >= 1,
                     "1_sig_P_value<=0.05_Log2FC>=1", "2_other_no_sig"),
      Cell_Type = sapply(Type, function(x) unlist(strsplit(x, split = "_"))[3])
    ) %>%
    mutate(
      Label = ifelse(Cell_Type == "NeuN" & Gene %in% NeuN_pro & Neg_log_P >= -log10(0.05) & Difference >= 1,
                     "0_sig_label", Label),
      Label = ifelse(Cell_Type == "NF200" & Gene %in% NF200_pro & Neg_log_P >= -log10(0.05) & Difference >= 1,
                     "0_sig_label", Label),
      Label = ifelse(Cell_Type == "GFAP" & Gene %in% GFAP_pro & Neg_log_P >= -log10(0.05) & Difference >= 1,
                     "0_sig_label", Label),
      size = 1
    )
  
  
  df <- df %>%
    mutate(
      size = if_else(Label == "0_sig_label", size_num, size)
    )
  
  
  
  
  df <- df %>%
    dplyr::mutate(Cluster_num = as.numeric(factor(Cluster)) - 1)
  
  #df$Gene_slim = sapply(df$GeneID, function(x){unlist(strsplit(x, split = "_"))[2]})
  
  df$Gene_slim = sapply(df$Gene, function(x){unlist(strsplit(x, split = "_"))[1]})
  
  
  
  df_top10 = df[df$size==size_num,]
  df_not_top10 = df[df$size == 1,]
  
  df_not_top10_v2_nosig = df_not_top10 %>% dplyr::filter(Label == "2_other_no_sig")
  df_not_top10_v1_sig = df_not_top10 %>% dplyr::filter(Label == "1_sig_P_value<=0.05_Log2FC>=1")
  
  df_sig_num = as.data.frame(table(df$Cluster, df$Label))
  
  colnames(df_sig_num) = c("Cluster", "Label", "Counts")
  
  df_sig_num$Cluster = as.character(df_sig_num$Cluster)
  df_sig_num$Label = as.character(df_sig_num$Label)
  
  
  df_filtered <- df_sig_num %>%
    dplyr::filter(Label != "2_other_no_sig")
  
  df_filtered$Label = NULL
  
  df_filtered = df_filtered %>%   
    dplyr::group_by(Cluster) %>%                   
    dplyr::summarise(Total_Counts = sum(Counts))  
  
  
  
  
  
  min_max_per_cluster <- df %>%
    group_by(Cluster) %>%
    summarise(
      Max_Log2FC = max(Log2FC),
      Min_Log2FC = min(Log2FC),
      .groups = 'drop'
    )
  
  
  
  
  
  
  if(!BG_length_equal)
    dfbar<-data.frame(x=seq(0,compare_num - 1,by=1),
                      y=  min_max_per_cluster$Max_Log2FC+0.1)
  
  dfbar1<-data.frame(x=seq(0,compare_num - 1,by=1),
                     y=min_max_per_cluster$Min_Log2FC - 0.1)
  
  if(BG_length_equal)
  {
    max_value = max(min_max_per_cluster$Max_Log2FC) + 0.1
    min_value = min(min_max_per_cluster$Min_Log2FC) - 0.1
    
    dfbar<-data.frame(x=seq(0,compare_num - 1,by=1), y = max_value)
    dfbar1 <- data.frame(x=seq(0,compare_num - 1,by=1), y = min_value)
    
  }
  
  p2 = ggplot()+
    geom_col(data = dfbar,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = dfbar1,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = df_not_top10_v2_nosig,
                aes(x = Cluster_num, y = Log2FC, color = Label),
                size = 0.8,
                #width =0.4
                position = position_jitter(seed = 123))+
    
    geom_jitter(data = df_not_top10_v1_sig,
                aes(x = Cluster_num, y = Log2FC, color = Label),
                size = 0.8,
                #width =0.4
                position = position_jitter(seed = 123))+
    
    geom_jitter(data = df_top10,
                aes(x = Cluster_num, y = Log2FC, color = Label),
                #size = 1.6,
                size = size_num,
                #width =0.4
                position = position_jitter(seed = 123))
  
  #p2
  
  
  dfcol<-data.frame(x=c(0:(compare_num - 1)),
                    y=0,
                    label=c(0:(compare_num - 1)),
                    label_name = sort(unique(df$Cluster, decreasing = FALSE)))
  
  
  mycol <- nature_color_top7[1:compare_num]
  
  #"#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF"
  
  p3 <- p2 + geom_tile(data = dfcol,
                       aes(x=x,y=y),
                       height=0.4,
                       color = "black",
                       fill = mycol,
                       alpha = 0.6,
                       show.legend = F)
  
  p3 <- p3 +
    scale_color_manual(name=NULL,
                       values = c("yellow","red","grey"))
  
  
  p3 <- p3+
    labs(x="Cluster",y="average logFC")+
    geom_text(data=dfcol,
              aes(x=x,y=y,label=label_name),
              size =3,
              color ="black")
  
  
  
  
  
  p3 <- p3+
    theme_minimal()+
    theme(
      axis.title = element_text(size = 13,
                                color = "black",
                                face = "bold"),
      axis.line.y = element_line(color = "black",
                                 size = 1.2),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 15)
    )
  
  #p3
  
  p4 <- p3+
    geom_text(
      data=df_top10,
      aes(x=Cluster_num,y=Log2FC,label=Gene_slim),
      size = 2.5, col = "black",
      position = position_jitter(seed = 123)
    )
  
  df_text = df_filtered
  df_text$x = seq(0.2,0.2 + compare_num - 1,by=1)
  df_text$y = 12
  
  
  
  p5 <- p4 + geom_text(data = df_text, aes(x = x, y = y, label = paste("sig_num: ",Total_Counts, sep = "")), size = 3, fontface = "bold") 
  
  print(p5)
  
  
  
}

process_data_zone_func = function(dat1 = df_v1_Amy, file_name = "v1_Amy")
{
  #remove BSA proteins
  #Select columns that calculate P value, Log2FC, etc
  #Rename columns
  #...
  
  dat1 <- dat1 %>%
    filter(PG.Genes != "(Bos") %>%
    dplyr::select(matches("Student|Gene|PG.ProteinGroups")) %>%
    dplyr::select(c(1, 3, 5, 7, 8)) %>%
    rename(
      Significant = 1,
      P_value = 2,
      Difference = 3,
      ID = 4,
      Gene = 5
    ) %>%
    mutate(
      Type = file_name,
      P_value = as.numeric(P_value),
      Difference = as.numeric(Difference)
    )
  
  return(dat1)    
}