#'@time 202403
#'@author Yuan
#'@desc This is template for performing enrichment analysis of proteomics data. 
#'@function  1-compareCluster for enrichment analysis; 
#            2-GSEA for enrichment analysis; 

library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)

# -------------------------------------------------------------------------
#                I - Enrichment analysis with compareCluster              #
# -------------------------------------------------------------------------
# Desc: GO enrichment of cell-type enriched proteins
# Significance was calculated by hypergeometric distribtution followed by Benjamini-Hochberg correction for multiple comparison (adjusted p-value < 0.05) 
# Input: Obj_list$aov_ct3_df
# Output: Dot plot; Related table
{
    rm(list=ls())        
    Obj_list = readRDS(file="dataset_03.rds") 
    
    meta_df = Obj_list$meta_df %>% filter(CellType %in% Obj_list$ct3_order)                
    data_df = Obj_list$impute_df[,meta_df$SampleID]
    dep_df = Obj_list$aov_ct3_df %>% dplyr::filter(pvalue<0.05 & log2FC>log2(1.2))        
    dep_df$pid = gsub(dep_df$pids,pattern=";.*",replacement="")        
    
    tran_df = bitr(dep_df$pids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    plot_df = dep_df %>% merge(tran_df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)    
    plot_df$CellType = factor(plot_df$CellType, levels=Obj_list$ct3_order)    
    plot_df = plot_df %>% dplyr::arrange(CellType)
    
    bg_df = Obj_list$quantified_df
    row.names(bg_df) = gsub(row.names(bg_df),pattern=";.*",replacement="")
    bg_pid_df = bitr(row.names(bg_df), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    
    formula_res <- compareCluster(data = plot_df,ENTREZID~CellType,fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,minGSSize=3,maxGSSize=2000,
                                ont = "BP", pAdjustMethod = "BH",#universe = bg_pid_df$ENTREZID,
                                pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                
    )
    formula_res_cutoff = formula_res
    formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<=0.05,]
    
    write.csv(x = formula_res@compareClusterResult, file=paste0("write/Functional/Fig4_d_GOBP_",Sys.Date(),".csv"))  
    
    # selcted GO items
    selected_GO = formula_res_cutoff@compareClusterResult
    selected_items = c("digestion","pancreatic juice secretion","intestinal cholesterol absorption","oxidative phosphorylation","lipid digestion",
    "Ras protein signal transduction","regulation of lymphocyte mediated immunity","regulation of stress fiber assembly","negative regulation of natural killer cell mediated cytotoxicity","negative regulation of actin filament bundle assembly",
    "regulation of response to wounding","barbed-end actin filament capping","regulation of wound healing","negative regulation of actin filament depolymerization","mitochondrial translation"
    
    )    
    formula_res_cutoff@compareClusterResult = selected_GO %>% filter(Description %in% selected_items)
    formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)
    p2=dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=14,color="log10P",size="count") + 
        theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
        scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30))

    # save dot plot
    pdf(file=paste0("write/Functional/ct3_GOBP_Selected_compareCluster_",Sys.Date(),".pdf"),width = 8,height = 6)
    print(p2)
    dev.off() 
  }

# -------------------------------------------------------------------------
#                    II - Enrichment analysis with GSEA                   #
# -------------------------------------------------------------------------
# Desc: GSEA enrichment of cell-type enriched proteins
# Significance was calculated by hypergeometric distribtution followed by Benjamini-Hochberg correction for multiple comparison (adjusted p-value < 0.05) 
# Input: Obj_list$aov_ct3_df
# Output: GSEA plot; Related table
{
    rm(list=ls())        
    obj_list = readRDS(file = "obj_list_ct8_20240316.rds")      
    meta_df = obj_list$meta_df        

    dep_df = obj_list$dep_df
    dep_df$`-logp` = -log10(dep_df$pvalue)
    dep_df$rank = dep_df$log2fc * (dep_df$`-logp`)
    dep_df = dep_df %>% arrange(desc(rank))
    
    temp_geneList = dep_df$rank
    names(temp_geneList) = dep_df$pid
    re = clusterProfiler::gseGO(geneList =temp_geneList,OrgDb = org.Mm.eg.db,
                                keyType = "UNIPROT",minGSSize = 3,maxGSSize = 2000,
                                by = "fgsea",pvalueCutoff = 1,seed=123)
    
    saveRDS(re, paste0("write/Functional/ct8_DEP_GSEA_",Sys.Date(),".rds"))
    write.csv(x = re@result, file=paste0("write/Functional/Fig6h_ct8_DEP_GSEA_",Sys.Date(),".csv"))  

    # PLOT
    # re = readRDS("output_immune/GSEA.rds")
    selected_items = c("myeloid leukocyte activation","myeloid leukocyte cytokine production")
        
    selected.re = re
    selected.re@result = selected.re@result %>% filter(Description %in% selected_items)
    selected.re@result$Description = factor(selected.re@result$Description,levels = selected_items)
    selected.re@result = selected.re@result %>% arrange(Description)
            
    p1=gseaplot2(selected.re, "GO:0002274", base_size=14,pvalue_table = T,col = "firebrick")# + theme_bw()
    p2=gseaplot2(selected.re, "GO:0061082", base_size=14,pvalue_table = T,col = "firebrick")# + theme_bw()        
        
    pdf(file=paste0("write/Functional/Fig6h_ct8_DEP_GSEA_",Sys.Date(),".pdf"),width = 8,height = 6)
    print(p1)
    print(p2)
    dev.off() 

}
