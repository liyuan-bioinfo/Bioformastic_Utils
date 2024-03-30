library(pheatmap)
library(dplyr)
# Pheatmap with balanced legend breaks
{
    rm(list=ls())
    obj_list = readRDS(file="dataset_demo.rds")
    meta_df = obj_list$meta
    data_df = obj_list$impute[,meta_df$SampleID]
    #data_df$pid = obj_list$anno[row.names(data_df),"pid"]
    
    dep_df = obj_list$PDAC_stage %>% filter(log2FC>log2(2) & Pvalue<0.05) # filter sig. proteins from ANOVA test
    data_df = data_df.df[dep_df$pid,] # use Sig. proteins for plotting heatmap

    # scale data_df
    scale_data_df = t(apply(data_df, 1, scale)) %>% as.data.frame()
    names(scale_data_df) =  names(data_df)

    # set break range for pheatmap 
    range_all <- range(c(-2, 2))    
    my_palette <- colorRampPalette(c("navy", "white","firebrick3"))(n=100)
    breaks = seq(range_all[1], range_all[2], length.out = 101)

    # plot heatmap
    p1=pheatmap::pheatmap(scale_data_df,
                            scale = "none",show_rownames = F,show_colnames = T,
                            legend = T,border_color = NA,                        
                            color = my_palette, 
                            breaks = breaks,
                            cluster_cols = F,cluster_rows=T,
                            cellwidth = 25,cellheight = 0.3,fontsize = 14,fontsize_number = 10,
                            silent = F)
    # save plot
    pdf(file=paste0("write/Figure_4c.pdf"),width = 7,height = 8)
    print(p1)
    dev.off()
}
