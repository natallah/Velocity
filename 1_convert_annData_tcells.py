import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("t_cell_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("t_cells_metadata_renamed.csv")

# load gene names:
with open("t_cell_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()


# set anndata observations and index obs by barcodes (name them cells since the index cannot have the same name as a column name), var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.obs.index.name = "cells"
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("t_cells_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
#resolution for macs was 0.5
#resolution for t cells was 0.3
adata.obs['seurat_clusters'] = list(adata.obs['seurat_clusters'])
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

gyr = ['#FF6F00FF','#C71000FF','#008EA0FF','#8A4198FF','#5A9599FF','#FF6348FF','#84D7E1FF','#FF95A8FF','#3D3B25FF','#ADE2D0FF','#1A5354FF','#3F4041FF']

sc.pl.umap(adata, color=['seurat_clusters'], palette=gyr, frameon=False, save="_t_cells_orig")

# save dataset as anndata format
adata.write('t_cells.h5ad')

# reload dataset
#adata = sc.read_h5ad('t_cells.h5ad')
