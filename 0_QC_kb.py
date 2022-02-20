import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib
import matplotlib.pyplot as plt

adata003 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/003_small_kb/counts_unfiltered/003_small.h5ad")
adata004 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/004_small_kb/counts_unfiltered/004_small.h5ad")
adata006 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/006_small_kb/counts_unfiltered/006_small.h5ad")
adata007 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/007_small_kb/counts_unfiltered/007_small.h5ad")
adata008 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/008_small_kb/counts_unfiltered/008_small.h5ad")
adata009 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/009_small_kb/counts_unfiltered/009_small.h5ad")
adata010 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/010_small_kb/counts_unfiltered/010_small.h5ad")
adata013 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/013_small_kb/counts_unfiltered/013_small.h5ad")
adata1144 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/1144_small_kb/counts_unfiltered/1144_small.h5ad")
adata1196 = anndata.read_h5ad("/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/1196_small_kb/counts_unfiltered/1196_small.h5ad")

adata003.var["gene_id"] = adata003.var.index.values
adata004.var["gene_id"] = adata004.var.index.values
adata006.var["gene_id"] = adata006.var.index.values
adata007.var["gene_id"] = adata007.var.index.values
adata008.var["gene_id"] = adata008.var.index.values
adata009.var["gene_id"] = adata009.var.index.values
adata010.var["gene_id"] = adata010.var.index.values
adata013.var["gene_id"] = adata013.var.index.values
adata1144.var["gene_id"] = adata1144.var.index.values
adata1196.var["gene_id"] = adata1196.var.index.values

t2g = pd.read_csv("homo_sapiens/transcripts_to_genes.txt", header=None, names=["tid", "gene_id", "gene_name"], sep="\t")
t2g.index = t2g.gene_id
t2g = t2g.loc[~t2g.index.duplicated(keep='first')]

adata003.var["gene_name"] = adata003.var.gene_id.map(t2g["gene_name"])
adata003.var.index = adata003.var["gene_name"]

adata004.var["gene_name"] = adata004.var.gene_id.map(t2g["gene_name"])
adata004.var.index = adata004.var["gene_name"]

adata006.var["gene_name"] = adata006.var.gene_id.map(t2g["gene_name"])
adata006.var.index = adata006.var["gene_name"]

adata007.var["gene_name"] = adata007.var.gene_id.map(t2g["gene_name"])
adata007.var.index = adata007.var["gene_name"]

adata008.var["gene_name"] = adata008.var.gene_id.map(t2g["gene_name"])
adata008.var.index = adata008.var["gene_name"]

adata009.var["gene_name"] = adata009.var.gene_id.map(t2g["gene_name"])
adata009.var.index = adata009.var["gene_name"]

adata010.var["gene_name"] = adata010.var.gene_id.map(t2g["gene_name"])
adata010.var.index = adata010.var["gene_name"]

adata013.var["gene_name"] = adata013.var.gene_id.map(t2g["gene_name"])
adata013.var.index = adata013.var["gene_name"]

adata1144.var["gene_name"] = adata1144.var.gene_id.map(t2g["gene_name"])
adata1144.var.index = adata1144.var["gene_name"]

adata1196.var["gene_name"] = adata1196.var.gene_id.map(t2g["gene_name"])
adata1196.var.index = adata1196.var["gene_name"]

