import scanpy as sc
import anndata
import MultiMAP
import numpy as np

### data integration with MultiMAP 

sc.settings.set_figure_params(dpi=150)
sc.settings.verbosity = 3
sc.settings.autosave=True

# load PBMC resting and stimulated filtered dataset
pb = sc.read('pbmc.h5ad')
# log normalization
sc.pp.normalize_total(pb, target_sum=10000)
sc.pp.log1p(pb)

#determine highly variable genes
sc.pp.highly_variable_genes(pb, min_mean=0.0125, max_mean=3, min_disp=0.5)

#select highly variable genes
pb = pb[:, pb.var.highly_variable]

# data integration and visualization
MultiMAP.Batch(pb, batch_key='orig.ident')

# visualization
sc.pl.embedding(pb, 'X_multimap', color='orig.ident')

