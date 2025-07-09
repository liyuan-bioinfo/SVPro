# All plugins can be split into 3 parts
# 1. Reading the command line arguments provided by Perseus and parsing the data.
# 2. Perform the desired functionality on the data.
# 3. Write the results to the expected locations in the Perseus formats.

# 1. Parse command line arguments passed in from Perseus,
# including input file and output file paths.
args = commandArgs(trailingOnly=TRUE)
#if (length(args) != 2) {
#  stop("Do not provide additional arguments!", call.=FALSE)
#}
paramFile = args[1]
inFile <- args[2]
outFile <- args[3]


# Use PerseusR to read and write the data in Perseus text format.
library(PerseusR)
mdata <- read.perseus(inFile)

# The mdata object can be easily deconstructed into a number of different
# data frames. Check reference manual or help() for full list.
mainMatrix <- main(mdata)
aCols <- annotCols(mdata)

# 2. Run any kind of analysis on the extracted data.
#df <- mainMatrix + 1

#test
#pg = read.delim("L:\\Ke_mi\\Protein_annotaion\\From_Peizhong\\z_separate_search\\proteinGroups_BTLA.txt", sep = "\t", header = T, stringsAsFactors = F, encoding = "latin1")
#pg.v1 = dplyr::select(pg, matches("LFQ"))
#x = pg.v1

#add paramFile a random value
if(paramFile == "")
{
  paramFile = "a"
}


Timecourse.median.normalization.proteinGroup.func <- function(x = mainMatrix, bait_index = paramFile)
{
  df.start = x
  #lfq.col = grep("LFQ.intensity", colnames(x))
  lfq.col = 1:ncol(x)
  x = x[,lfq.col]
  sample.num = ncol(x) 
  
  
  
  x_med = apply(x, 2, function(x){median(x)})
  
  if(any(x_med == 0))
  {
    print("median of some columns are equal to zero")
  }
  
  if(!any(x_med == 0))
  {
    x_med_m = mean(x_med)
    
    STS1_cof = x_med_m/x_med
    
    #STS1_cof <- STS1_mean/as.numeric(x[STS1_index,c(1:sample.num)])
    STS1_cof_mat <- diag(x = STS1_cof,nrow = sample.num,ncol = sample.num)
    tpl_mat <- t(as.matrix(x[,c(1:sample.num)]))
    
    tpl_mat_prod <- STS1_cof_mat%*%tpl_mat
    tpl_mat_prod <- t(tpl_mat_prod)
    
    #rownames(tpl_mat_prod) <- x$Majority.protein.IDs
    colnames(tpl_mat_prod) <- colnames(x)[c(1:sample.num)]
    
    df.start[,lfq.col] = tpl_mat_prod
    return(df.start)
    #apply(tpl_mat_prod,2,median)
  }


  
}


#v2: calculate median after removing zero value.

Timecourse.median.normalization.proteinGroup.v2.nozero.func <- function(x = mainMatrix, bait_index = paramFile)
{
  df.start = x
  #lfq.col = grep("LFQ.intensity", colnames(x))
  lfq.col = 1:ncol(x)
  x = x[,lfq.col]
  sample.num = ncol(x) 
  
  x[x == 0] = NA
  
  x_med = apply(x, 2, function(x){median(x,na.rm = T)})
  
  if(any(x_med == 0))
  {
    print("median of some columns are equal to zero")
  }
  
  if(!any(x_med == 0))
  {
    x_med_m = mean(x_med)
    
    STS1_cof = x_med_m/x_med
    STS1_cof = as.vector(STS1_cof)
    
    for(i in 1:ncol(x))
    {
      x[,i] = x[,i] * STS1_cof[i]
    }
    
   
    df.start[,lfq.col] = x
    
    #df.start[is.na(df.start)] = 0
    
    return(df.start)
    
    
  }
  
  
  
}



#df_v2 = Timecourse.median.normalization.proteinGroup.func(x = mainMatrix, bait_index = paramFile)
df_v2 = Timecourse.median.normalization.proteinGroup.v2.nozero.func(x = mainMatrix, bait_index = paramFile)

#df_v2 = mainMatrix

#df_v3 = cbind(df_v2, aCols)
#df_v3 = as.data.frame(df_v3)

# 3. Create a matrixData object which can be conveniently written to file
# in the Perseus txt format.
outMdata <- matrixData(main=df_v2, annotCols = annotCols(mdata), annotRows = annotRows(mdata))
write.perseus(outMdata, outFile)