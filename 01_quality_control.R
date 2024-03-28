#'@time 202403
#'@author Yuan
#'@desc This is template for performing quality control of proteomics data. 

library(dplyr)
library(RColorBrewer)

setwd("")

# -------------------------------------------------------------------------
#                            I-Basic pre-treat                            #
# -------------------------------------------------------------------------
# Desc: Create .RDS and save files from Perseus output, including Quantified, Filter and Impute
# Output: .RDS
{
    rm(list=ls())

    # prepare meta file of experiment design
    meta_df = read.delim("meta.txt",header=T,check.names=FALSE) 
    ct3_order = c("Acinar","Lymph","Tumor")
    meta_df$CellType = factor(meta_df$CellType, levels=ct3_order)
    meta_df = meta_df %>% arrange(CellType)

    # prepare annotation file of Protein to GeneNames
    # Entry of Protein Groups as Row.Names
    anno_df = read.delim("anno.txt",header=T,check.names=FALSE) 
    row.names(anno_df) = anno_df$pids
    anno_df = anno_df %>% dplyr::select(pids, genes)
    anno_df$gene = gsub(xxx)
    anno_df$pid = gsub(xxx)
    

    # prepare Quantified file
    quantified_df = read.delim("Quantified_2917.txt",header=T,check.names=FALSE) #6406 * 69
    quantified_df = quantified_df[,meta_df$SampleID]#re-order
    names(quantified_df) = meta_df$SampleID #re-name

    # prepare Impute file
    impute_df = read.delim("Impute_2608.txt",header=T,check.names=FALSE) #5900 * 69
    row.names(impute_df) = impute_df$Protein
    impute_df = impute_df[,meta_df$SampleID]
    names(impute_df) = meta_df$SampleID

    # save  
    Obj_list = list() # save RDS
            
    Obj_list$quantified_df = quantified_df
    Obj_list$impute_df = impute_df    
    Obj_list$meta_df = meta_df
    Obj_list$anno_df = anno_df
    Obj_list$ct3_order = ct3_order
    
    saveRDS(Obj_list, file="dataset_01.rds")       
}
