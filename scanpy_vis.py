# 设置palette以及categories
# Panel-1, umap of single-cell proteomics 
# cell-type proteomics data
os.chdir("...")
sc_adata = sc.read_h5ad("...h5ad")

# basic filter
sc.pp.normalize_total(sc_adata)
sc.pp.log1p(sc_adata)
sc.pp.filter_cells(sc_adata, min_genes=10)
sc.pp.filter_genes(sc_adata, min_cells=3)

# set categories
# sc_adata.obs["celltype"] = sc_adata.obs["celltype"].astype("category")
# sc_adata.obs["celltype"] = sc_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)


# run umap
sc.pp.neighbors(sc_adata, n_pcs=20, n_neighbors=10)
sc.tl.umap(sc_adata)

# cluster or not
# sc.tl.leiden(sc_adata,resolution=0.1)
# sc.pl.umap(sc_adata, color=['celltype'],size=100)

# set color
my_palette = ["#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF","#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"]
sc_adata.uns["celltype_colors"] = my_palette

# add color to obs
# colors = pd.Series(my_palette, index=sc_adata.obs['celltype'].cat.categories)
# sc_adata.obs['color'] = sc_adata.obs['celltype'].map(colors)

sc.pl.umap(sc_adata, color=["celltype"], size=1000)

# plot using outer palette
# sc.pl.umap(sc_adata, color=["celltype"], size=10,palette=my_palette, save="HumanTonsil_panel_1_scp_ref46.pdf")
