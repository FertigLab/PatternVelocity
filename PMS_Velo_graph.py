"""

#Construct velocity graph in the Pattern-space: Project expression and sliced/unsliced data onto the 
Pattern space. Construct the velocity graph. Note that this provides both the magnitude and direction of velocity.

Parameters
----------
dat : AnnData object 
    contains spliced/unspliced RNA data
PMS_df: PatternMarkerStatistics dataframe 

Returns
-------
dat: PMS-weighted velocity graph   

"""

import scanpy as sc
import scvelo as scv
import anndata as ad
import pandas as pd
import numpy as np
import csv

dat = ad.read_h5ad('dataonly.h5ad')


dat.obs_names_make_unique()
dat.var_names_make_unique()

#Note: No filtering - keep all genes and cells - only normalize every cell by its initial size & insert log of counts
scv.pp.filter_and_normalize(dat)


# +
###Subset data - only subset data if limited by memory

##subset by cell-type:
#dat_subset = dat[dat.obs['umap2_CellType'].isin(['Horizontal Cells']),:] #.isin(['Early RPCs','Late RPCs']),:]#  
#del dat
#dat = dat_subset
###


##check data counts 
print(dat.obs["age"].value_counts())
print(dat.obs['umap2_CellType'].value_counts())


# -
##subset by age:
dat_subset = dat[dat.obs['age'].isin(['P0', 'P2', 'P5']),:] 
del dat
dat = dat_subset
print(dat.obs["age"].value_counts())
print(dat.obs['umap2_CellType'].value_counts())
###

PMS_df = pd.read_csv('patMarStatScores.csv')

PMS_df = PMS_df.set_index(PMS_df.columns[0])

PMS_mat = PMS_df.as_matrix()

# +

#subset velocity data to those genes that were in the CoGAPS run
common_genes = list(set(PMS_df.index) & set(dat.var.index))
dat_subset_genes = dat[:,dat.var.index.isin(common_genes)] 
PMS_mat_subset = PMS_mat[PMS_df.index.isin(common_genes),:]

# +

#PMS-weighted spliced/sunspliced:
PMS_wgt_spliced = np.transpose(PMS_mat_subset) @ np.transpose(dat_subset_genes.layers['spliced'])
PMS_wgt_unspliced = np.transpose(PMS_mat_subset) @ np.transpose(dat_subset_genes.layers['unspliced'])

# +
#PMS-weighted genes/Patterns:

PMS_wgt_genes = np.transpose(PMS_mat_subset) @ np.transpose(dat_subset_genes.X)


# -
new_var_df = pd.DataFrame(index=PMS_df.columns)

adata=sc.PMS_wgt_genes

dat_new = ad.AnnData(X=np.transpose(PMS_wgt_genes), obs=dat_subset_genes.obs, var=new_var_df)

# +
#replace original spliced/unsplced with PMS-matrix weighted spliced/unspliced 
dat_new.layers['spliced'] = np.transpose(PMS_wgt_spliced)

dat_new.layers['unspliced'] = np.transpose(PMS_wgt_unspliced)

# -

del dat

dat = dat_new

#Compute the first-order moments of spliced/unspliced abundances 
#This first computes neighborhood graph. Computes neighborhoods based on "connectivities" distance-metric (default method = umap)
n_pcs=10
n_neighbors=30
scv.pp.moments(dat, n_pcs=n_pcs, n_neighbors=n_neighbors)

#Compute magnitude of velocities in gene-expression states
scv.tl.velocity(dat)

dat.write('dataonlyPostProcessing_PMS_P0P2P5.h5ad')
Following is to compute the angle of the velocity vector - No need to run if only velocity magnitude is needed
#Construct umap embedding 
#sc.tl.umap(dat,n_components=2,min_dist=0.25)
#Check visualization 
#sc.pl.umap(dat,color="age",projection='2d',save="data_Horizontal.png")
#Compute velocity-graph (cell-cell Transition matrix) based on cosine similarities
#scv.tl.velocity_graph(dat, n_neighbors=n_neighbors) 
#Compute the each cell's velocity projection on the embedding
#scv.tl.velocity_embedding(dat, basis='umap')
#Check visualization 
#scv.pl.velocity_embedding_grid(dat, legend_loc="on data", arrow_size=1,arrow_length=1,basis='umap',color="age")
#scv.pl.velocity_embedding_stream(dat, legend_loc="on data", basis='umap', color="age",save="umap_Horizontal.png")




