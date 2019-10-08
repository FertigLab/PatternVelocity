"""

#Originally written by Loyal, edited by Zeinab so that only necessary functions are called

#Construct velocity graph on retina-only data 
Parameters
----------
dat : AnnData object 
    contains spliced/unspliced RNA data 

Returns
-------
datOnlyPostProcessing.h5ad 
    AnnData object containing the velocity graph  

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





#Compute the first-order moments of spliced/unspliced abundances 
#This first computes neighborhood graph. Computes neighborhoods based on "connectivities" distance-metric (default method = umap)
n_pcs=10
n_neighbors=30
scv.pp.moments(dat, n_pcs=n_pcs, n_neighbors=n_neighbors)

#Compute magnitude of velocities in gene-expression states
scv.tl.velocity(dat)

#Following is to compute the angle of the velocity vector - No need to run if only velocity magnitude is needed
#Construct umap embedding 
sc.tl.umap(dat,n_components=2,min_dist=0.25)
#Check visualization 
sc.pl.umap(dat,color="age",projection='2d',save="data_Horizontal.png")
#Compute velocity-graph (cell-cell Transition matrix) based on cosine similarities
scv.tl.velocity_graph(dat, n_neighbors=n_neighbors) 
#Compute the each cell's velocity projection on the embedding
scv.tl.velocity_embedding(dat, basis='umap')
#Check visualization 
scv.pl.velocity_embedding_grid(dat, legend_loc="on data", arrow_size=1,arrow_length=1,basis='umap',color="age")
scv.pl.velocity_embedding_stream(dat, legend_loc="on data", basis='umap', color="age",save="umap_Horizontal.png")

# +
#dat.write('dataonlyPostProcessing_P0P2P5.h5ad')
# + {}
#dat.write('dataonlyPostProcessing_Horizontal.h5ad')
# -




