import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("macs_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("macs_metadata_renamed.csv")

# load gene names:
with open("macs_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()


# set anndata observations and index obs by barcodes (name them cells since the index cannot have the same name as a column name), var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.obs.index.name = "cells"
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("macs_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
#resolution for macs was 0.5
#resolution for t cells was 0.3

adata.obs['seurat_clusters'] = list(adata.obs['seurat_clusters'])
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

gyr = ['#AE3121','#964D00', '#755F00','#3C6D00','#007700','#007E4A','#007F85','#0077B1']
sc.pl.umap(
        adata,
        color=['seurat_clusters'],
        frameon=False,
        palette= gyr, 
        save="_macs_orig")

# save dataset as anndata format
adata.write('macs.h5ad')

# reload dataset
#adata = sc.read_h5ad('macs.h5ad')
