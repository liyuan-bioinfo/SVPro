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

#Example: paramFile = "Grb2"

gene_col = grep("Gene", colnames(aCols))
#bait_index = which(aCols[,gene_col] == paramFile)

# 2. Run any kind of analysis on the extracted data.
#df <- mainMatrix + 1
Timecourse.bait.normalization.proteinGroup.func <- function(x = mainMatrix, bait_index)
{
  df.start = x
  #lfq.col = grep("LFQ.intensity", colnames(x))
  lfq.col = 1:ncol(x)
  x = x[,lfq.col]
  sample.num = ncol(x) 
  STS1_index <- bait_index
  STS1_mean <- mean(as.numeric(x[STS1_index,c(1:sample.num)]))
  
  STS1_cof <- STS1_mean/as.numeric(x[STS1_index,c(1:sample.num)])
  STS1_cof_mat <- diag(x = STS1_cof,nrow = sample.num,ncol = sample.num)
  tpl_mat <- t(as.matrix(x[,c(1:sample.num)]))
  
  tpl_mat_prod <- STS1_cof_mat%*%tpl_mat
  tpl_mat_prod <- t(tpl_mat_prod)
  
  #rownames(tpl_mat_prod) <- x$Majority.protein.IDs
  colnames(tpl_mat_prod) <- colnames(x)[c(1:sample.num)]
  
  df.start[,lfq.col] = tpl_mat_prod
  return(df.start)
  
}


#df_v3 = cbind(df_v2, aCols)
#df_v3 = as.data.frame(df_v3)

# 3. Create a matrixData object which can be conveniently written to file
# in the Perseus txt format.
#outMdata <- matrixData(main=df_v2)

Timecourse.bait.normalization.proteinGroup.v2.all.biotinylated.proteins.func <- function(x = mainMatrix, biotin_index)
{
  df.start = x
  #lfq.col = grep("LFQ.intensity", colnames(x))
  lfq.col = 1:ncol(x)
  x = x[,lfq.col]
  sample.num = ncol(x) 
  
  df_biotin = x[biotin_index, c(1:sample.num)]
  
  df_biotin_sum = apply(df_biotin, 2, function(x){sum(x, na.rm = TRUE)})
  
  df_biotin_mean = mean(df_biotin_sum)
  
  STS1_cof = df_biotin_mean/df_biotin_sum

  
  
  #STS1_index <- bait_index
  #STS1_mean <- mean(as.numeric(x[STS1_index,c(1:sample.num)]))
  #STS1_cof <- STS1_mean/as.numeric(x[STS1_index,c(1:sample.num)])
  STS1_cof_mat <- diag(x = STS1_cof,nrow = sample.num,ncol = sample.num)
  tpl_mat <- t(as.matrix(x[,c(1:sample.num)]))
  
  tpl_mat[is.na(tpl_mat)] = 0
  
  tpl_mat_prod <- STS1_cof_mat%*%tpl_mat
  tpl_mat_prod <- t(tpl_mat_prod)
  
  #rownames(tpl_mat_prod) <- x$Majority.protein.IDs
  colnames(tpl_mat_prod) <- colnames(x)[c(1:sample.num)]
  
  df.start[,lfq.col] = tpl_mat_prod
  
  #df.start[df.start == 0] = NA
  
  return(df.start)
  
}



if(paramFile != "All")
{
	bait_index = which(aCols[,gene_col] == paramFile)
	df_v2 = Timecourse.bait.normalization.proteinGroup.func(x = mainMatrix, bait_index)
	outMdata <- matrixData(main=df_v2, annotCols = annotCols(mdata), annotRows = annotRows(mdata))
	write.perseus(outMdata, outFile)

}

#Five biotniylation proteins
#Acaca, Pc, Mccc1, Pcca, Acacb

if(paramFile == "All")
{
	  ebp = c("ACACA","PC","PCCA","MCCC1", "ACACB", 
          "Acaca","Pc","Pcca","Mccc1", "Acacb")

	biotin_index = which(aCols[,gene_col] %in% ebp)
	df_v2 = Timecourse.bait.normalization.proteinGroup.v2.all.biotinylated.proteins.func(x = mainMatrix, biotin_index)
	outMdata <- matrixData(main=df_v2, annotCols = annotCols(mdata), annotRows = annotRows(mdata))
	write.perseus(outMdata, outFile)
	
  #z <- df[df$Gene %in% ebp, gene_col]
}

#outMdata <- matrixData(main=df_v2, annotCols = annotCols(mdata))


