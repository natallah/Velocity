import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

adata = sc.read_h5ad('macs.h5ad')

#read in cell IDs
sample_obs = pd.read_csv("macs_cellID_obs.csv")

#get cell IDs for each sample
cellID_obs_003 = sample_obs[sample_obs["x"].str.contains("_003_small")]
cellID_obs_004 = sample_obs[sample_obs["x"].str.contains("_004_small")]
cellID_obs_006 = sample_obs[sample_obs["x"].str.contains("_006_small")]
cellID_obs_007 = sample_obs[sample_obs["x"].str.contains("_007_small")]
cellID_obs_008 = sample_obs[sample_obs["x"].str.contains("_006_small")]
cellID_obs_009 = sample_obs[sample_obs["x"].str.contains("_090_small")]
cellID_obs_010 = sample_obs[sample_obs["x"].str.contains("_010_small")]
cellID_obs_013 = sample_obs[sample_obs["x"].str.contains("_013_small")]
cellID_obs_1144 = sample_obs[sample_obs["x"].str.contains("_1144_small")]
cellID_obs_1196 = sample_obs[sample_obs["x"].str.contains("_1196_small")]

# load loom files for spliced/unspliced matrices for each sample:

ldata003 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/003_small_kb/counts_unfiltered/003_small.h5ad')
ldata004 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/004_small_kb/counts_unfiltered/004_small.h5ad')
ldata006 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/006_small_kb/counts_unfiltered/006_small.h5ad')
ldata007 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/007_small_kb/counts_unfiltered/007_small.h5ad')
ldata008 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/008_small_kb/counts_unfiltered/008_small.h5ad')
ldata009 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/009_small_kb/counts_unfiltered/009_small.h5ad')
ldata010 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/010_small_kb/counts_unfiltered/010_small.h5ad')
ldata013 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/013_small_kb/counts_unfiltered/013_small.h5ad')
ldata1144 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/1144_small_kb/counts_unfiltered/1144_small.h5ad')
ldata1196 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/1196_small_kb/counts_unfiltered/1196_small.h5ad')

print(ldata003)
print(ldata003.var)

# rename barcodes in order to merge:
barcodes = [bc[0:len(bc)] + '_003_small' for bc in ldata003.obs.index.tolist()]
ldata003.obs.index = barcodes
ldata003 = ldata003[np.isin(ldata003.obs.index, cellID_obs_003["x"])]

barcodes = [bc[0:len(bc)] + '_004_small' for bc in ldata004.obs.index.tolist()]
ldata004.obs.index = barcodes
ldata004 = ldata004[np.isin(ldata004.obs.index, cellID_obs_004["x"])]

barcodes = [bc[0:len(bc)] + '_006_small' for bc in ldata006.obs.index.tolist()]
ldata006.obs.index = barcodes
ldata006 = ldata006[np.isin(ldata006.obs.index, cellID_obs_006["x"])]

barcodes = [bc[0:len(bc)] + '_007_small' for bc in ldata007.obs.index.tolist()]
ldata007.obs.index = barcodes
ldata007 = ldata007[np.isin(ldata007.obs.index, cellID_obs_007["x"])]

barcodes = [bc[0:len(bc)] + '_008_small' for bc in ldata008.obs.index.tolist()]
ldata008.obs.index = barcodes
ldata008 = ldata008[np.isin(ldata008.obs.index, cellID_obs_008["x"])]

barcodes = [bc[0:len(bc)] + '_090_small' for bc in ldata009.obs.index.tolist()]
ldata009.obs.index = barcodes
ldata009 = ldata009[np.isin(ldata009.obs.index, cellID_obs_009["x"])]

barcodes = [bc[0:len(bc)] + '_010_small' for bc in ldata010.obs.index.tolist()]
ldata010.obs.index = barcodes
ldata010 = ldata010[np.isin(ldata010.obs.index, cellID_obs_010["x"])]

barcodes = [bc[0:len(bc)] + '_003_small' for bc in ldata013.obs.index.tolist()]
ldata013.obs.index = barcodes
ldata013 = ldata013[np.isin(ldata013.obs.index, cellID_obs_013["x"])]

barcodes = [bc[0:len(bc)] + '_1144_small' for bc in ldata1144.obs.index.tolist()]
ldata1144.obs.index = barcodes
ldata1144 = ldata1144[np.isin(ldata1144.obs.index, cellID_obs_1144["x"])]

barcodes = [bc[0:len(bc)] + '_1196_small' for bc in ldata1196.obs.index.tolist()]
ldata1196.obs.index = barcodes
ldata1196 = ldata1196[np.isin(ldata1196.obs.index, cellID_obs_1196["x"])]

# make variable names unique
#ldata003.obs_names_make_unique()
#ldata004.obs_names_make_unique()
#ldata006.obs_names_make_unique()


# concatenate the three loom
ldata = ldata003.concatenate([ldata004, ldata006,ldata007,ldata008,ldata009,ldata010,ldata013,ldata1144,ldata1196])

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

#plot a umap
gyr = ['#AE3121','#964D00', '#755F00','#3C6D00','#007700','#007E4A','#007F85','#0077B1']
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', palette= gyr, title='', save='_macs_celltypes.pdf')

#inspect spliced and unspliced reads
#scv.pl.proportions(adata,save='_macs_celltypes.pdf')

# pre-processing
#scv.pp.filter_and_normalize(adata)


