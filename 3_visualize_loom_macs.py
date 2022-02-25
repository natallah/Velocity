import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

gyr = ['#AE3121','#964D00', '#755F00','#3C6D00','#007700','#007E4A','#007F85','#0077B1']


adata = scv.read('macs_velocity.h5ad')

scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', palette=gyr, save='embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap',color='seurat_clusters', palette=gyr,save='embedding_stream.pdf')

