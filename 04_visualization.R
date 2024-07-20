library(pheatmap)
library(dplyr)

# set color
{
    ggsci::pal_npg()(10)    
}

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

library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
# dot plot after enrichment analysis
{
    rm(list=ls())
    obj_list = readRDS(file = "sp_dataset_03.rds")
    dep_df = obj_list$PDAC_stage %>% filter(Pvalue<0.05 & log2FC>log2(2))
    dep_df$ENTREZID = obj_list$anno[dep_df$pid,"GeneID"] # set region of dep_df

    formula_res <- compareCluster(data = dep_df,ENTREZID~region,fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
    )
    
    formula_res_cutoff = formula_res
    formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue < 0.05,] # this may set the cutoff of p.adjust
    formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)

    selected_GO = formula_res_cutoff@compareClusterResult
    selected_items = c("digestion","pancreatic juice secretion","intestinal cholesterol absorption","oxidative phosphorylation","lipid digestion",
                "Ras protein signal transduction","regulation of lymphocyte mediated immunity","regulation of stress fiber assembly","negative regulation of natural killer cell mediated cytotoxicity","regulation of actin filament bundle assembly",
                "positive regulation of response to wounding","barbed-end actin filament capping","positive regulation of wound healing","negative regulation of actin filament depolymerization","mitochondrial translation"
                )
    formula_res_cutoff@compareClusterResult = selected_GO %>% filter(Description %in% selected_items)

    # begion dot plot using enrichplot
    p1=dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=14,color="log10P",size="count") + 
      theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
      scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30))
    
    pdf(file=paste0("write/Figure_4d.pdf"),width = 8,height = 6)
    print(p1)
    dev.off()
}

# 设置major celltype and subtype的颜色
library(RColorBrewer)
display.brewer.all() #显示所有调色板
# 获得Set1的颜色
display.brewer.pal(9,"Set1")# "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3"
brewer.pal(9,"Set1")

colorRampPalette(c("white",brewer.pal(9,"Set1")[1]))(100)[c(100)] # PC
colorRampPalette(c("white",brewer.pal(9,"Set1")[4]))(100)[c(100,80,60,40)] #CAF
colorRampPalette(c("white",brewer.pal(9,"Set1")[2]))(100)[c(100,80,60,40)] # Lym
colorRampPalette(c("white",brewer.pal(9,"Set1")[3]))(100)[c(100,80,60,40,20)] # MYE
