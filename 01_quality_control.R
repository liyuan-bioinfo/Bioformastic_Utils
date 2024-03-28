#'@time 202403
#'@author Yuan
#'@desc This is template for performing quality control of proteomics data. 
#'@function  1-Create RDS after imputation of missing values; 
#            2-PCA plot before data filtering; 
#            3-Corr plot after imputation of missing values;

library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)

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

# -------------------------------------------------------------------------
#                            II-Quality Control                           #
# -------------------------------------------------------------------------
# PCA analysis
# Input: Obj_list$quantified_df; Obj_list$meta_df;
{
    rm(list=ls())
    Obj_list = readRDS("dataset_01.rds")    
    data_df = Obj_list$quantified_df
    meta_df = Obj_list$meta_df
    
    ## PCA analysis
    data_df = log2(data_df+1) %>% t() %>% as.data.frame() 

    pca_result <- prcomp(data_df,scale=T,center = T)
    plot_df = as.data.frame(pca_result$x)
    summ1 <- summary(pca_result)
    xlab1 <- paste0("PC1 (",round(summ1$importance[2,1]*100,2),"%)")
    ylab1 <- paste0("PC2 (",round(summ1$importance[2,2]*100,2),"%)")

    p1=ggplot(data = plot_df,aes(x = PC1,y = PC2,color = meta_df$CellType))+
        stat_ellipse(aes(fill = meta_df$CellType),
                        type = "norm",geom = "polygon",alpha = 0.25, color = NA)+ 
        geom_point(size = 2)+
        labs(x = xlab1,y = ylab1)+ guides(fill = "none")+ theme_bw()+            
        theme(plot.background = element_blank(),legend.background = element_blank(),
            panel.background = element_blank(),panel.grid = element_blank(),
            axis.text = element_text(size = 12),axis.title = element_text(size = 16),
            legend.text = element_text(size = 16))

    ## save plot
    pdf(file=paste0("write/Quantified/Fig_2i_PCA_",Sys.Date(),".pdf"),
        width = 7,height=5)
    print(p1)
    dev.off()    

    ## save table        
    write.csv(plot_df, file = paste0("write/Quantified/Fig_2i_PCA_PCA_",Sys.Date(),".csv"))
}

# Correlation analysis
# Input: Obj_list$quantified_df; Obj_list$meta_df;
{
    rm(list=ls())
    Obj_list = readRDS("dataset_01.rds")
    
    data_df = Obj_list$impute_df
    meta_df = Obj_list$meta_df

    plot_df = cor(data_df[])
    
    # search the max value
    row_max = c()
    for(i in 1:dim(plot_df)[1]){
        temp_row = plot_df[i,]
        row_max = c(row_max,max(temp_row[which(temp_row!=1)]))            
    }

    # replace "1" with max corr
    row_max = max(row_max)
    for(i in 1:dim(plot_df)[1]){
        temp_row = plot_df[i,]            
        plot_df[i,which(temp_row==1)]=row_max            
    }    

    p1 = pheatmap(plot_df,cluster_rows = T,cluster_cols = T,width = 7,height = 6,#cellwidth = 3,cellheight = 3,
                    show_rownames = T,show_colnames = T,border_color = NA,
                    annotation_names_row = F,annotation_names_col = F,                        
                    annotation_legend = F,fontsize_row = 10,fontsize_col = 10,
                    color = colorRampPalette(c("blue","white","red"))(100),silent = TRUE)

    # add legends
    llim = round(min(plot_df),2)
    ulim = round(max(plot_df),2)

    # Fig lengends           
    df <- data.frame(matrix(nrow = 10, ncol = 2))
    df[] <- rnorm(20)
    p2=ggplot(data = df, aes(x = X1, y = X2, fill = X2)) + geom_tile() + theme_void()+
        scale_fill_gradientn(colors=colorRampPalette(c("blue","white","red"))(100),
        limits = c(llim, ulim), breaks = round(seq(llim, ulim, length.out=3),2),guide = guide_colorbar(barwidth = 5, barheight = 1, 
                                                                                        ticks = TRUE, ticks.colour="black",ticks.linewidth=1/.pt,direction="horizontal",#"vertical",                                                                                            
                                                                                        label.theme = element_text(size = 8),
                                                                                        draw.ulim = TRUE, 
                                                                                        draw.llim = TRUE, 
                                                                                        draw.separator = TRUE,
                                                                                        title.position = "top"
                                                                                        )
        )     

    ## save plot
    pdf(file=paste0("write/Quantified/sFig_2b_Corr_",Sys.Date(),".pdf"),
        width = 7,height=5)
    print(p1)
    print(p2)
    dev.off()    
}  
