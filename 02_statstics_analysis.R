#'@time 202403
#'@author Yuan
#'@desc This is template for performing statistics analysis of proteomics data. 
#'@function  1-Create RDS after imputation of missing values; 
#            2-PCA plot before data filtering; 
#            3-Corr plot after imputation of missing values;

library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)

    # Desc: ANOVA analysis of celltype groups, including 14 cell types and 4 cell lineages
    # Input: Obj_list$impute_df; Obj_list$impute_mean_df
    # Output: Obj_list$aov_ct3_df
    # for Fig.5e,f;  Fig.6c,d, Extended Data Fig. 10 and Fig. 12
    {
        rm(list=ls())
        Obj_list = readRDS(file="dataset_01.rds") 

        meta_df = Obj_list$meta_df        
        data_df = Obj_list$impute_df
        data_mean_df = Obj_list$impute_mean_df
        
        ct_num = 3        
        
        ## log2 FC, to find enriched cell-types with mean abundance
        dep_df = data.frame()
        for(i in 1:ct_num){
            temp_df = data_mean_df[,i] - data_mean_df[,1:ct_num]
            temp_pid = names(which(rowSums(temp_df>log2(1.2)) == (ct_num-1)))
            temp_df2 = data.frame(pids=temp_pid)
            temp_df2$CellType = names(data_mean_df)[i]
            temp_df2$log2FC = apply(temp_df[temp_pid,-i],1,mean)        
            dep_df = rbind(dep_df,temp_df2)            
        }

        ## pvalue, with each sample
        row.names(dep_df) = dep_df$pids
        # dep_df = dep_df[row.names(data_df),] # re-order
        data_df = data_df[dep_df$pids,]
        protein_num = dim(data_df)[1]

        Pvalue = c()
        for(i in 1:protein_num){
            pid_df = data_df[i,] %>% t() %>% as.data.frame()
            names(pid_df) = "pid"
            pid_df$SampleId = row.names(pid_df)
            pid_df$CellType = meta_df$CellType
            model = aov(data=pid_df,pid~CellType)
            # temp_pvalue = summary(model)[[1]]$`Pr(>F)`[1]
            sig_df = TukeyHSD(model)$CellType #which ct is enriched
            enriched_ct = dep_df[i,"CellType"]
            enriched_ct_sig_df = sig_df[grep(row.names(sig_df),pattern=paste0("-",enriched_ct,"$|","^",enriched_ct,"-")),] #whethor this ct is sig. with other cts
            temp_pvalue = RecordTest::fisher.method(as.numeric(enriched_ct_sig_df[,"p adj"]))$"p.value"[[1]]  
            # temp_pvalue = max(enriched_ct_sig_df[,"p adj"])                                    
            Pvalue = c(Pvalue,temp_pvalue)            
        }        

        dep_df$pvalue = Pvalue
        dep_df$fdr = p.adjust(Pvalue,method = "BH")
