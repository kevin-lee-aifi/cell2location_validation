import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import seaborn as sns
import os
import pandas as pd
from functions import *

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

results_folder = '/home/jupyter/EXP-01058/cell2location_results'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# -----------------------------------------------------------------------------------------------
# Loading Visium HD data
# -----------------------------------------------------------------------------------------------

home_dir = "/home/jupyter"
vis_hd_bin_dir_ = "EXP-01058/TIS05740-001-021/outs/binned_outputs/square_008um"
bin_dir = os.path.join(home_dir, vis_hd_bin_dir_)

# Reading in visium HD data
hd_adata = preprocess_vis_adata(sc.read_visium(bin_dir))

# find mitochondria-encoded (MT) genes
hd_adata.var['MT_gene'] = [gene.startswith('MT-') for gene in hd_adata.var['gene_ids']]

# remove MT genes for spatial mapping (keeping their counts in the object)
hd_adata.obsm['MT'] = hd_adata[:, hd_adata.var['MT_gene'].values].X.toarray()
hd_adata = hd_adata[:, ~hd_adata.var['MT_gene'].values]

# -----------------------------------------------------------------------------------------------
# Loading reference model
# -----------------------------------------------------------------------------------------------

results_folder = '/home/jupyter/EXP-01058/cell2location_results'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

adata_file = f"{ref_run_name}/sc.h5ad"
rna_adata = sc.read_h5ad(adata_file)
ref_mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", rna_adata)

# -----------------------------------------------------------------------------------------------
# Extracting reference cell types signatures as a pd.DataFrame
# -----------------------------------------------------------------------------------------------

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in rna_adata.varm.keys():
    inf_aver = rna_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in rna_adata.uns['mod']['factor_names']]].copy()
else:
    inf_aver = rna_adata.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in rna_adata.uns['mod']['factor_names']]].copy()
inf_aver.columns = rna_adata.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# -----------------------------------------------------------------------------------------------
# Cell2location: spatial mapping
# -----------------------------------------------------------------------------------------------

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(hd_adata.var_names, inf_aver.index)
hd_adata = hd_adata[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=hd_adata, batch_key="in_tissue")

# create and train the model
mod = cell2location.models.Cell2location(
    hd_adata, cell_state_df=inf_aver, 
    # the expected average cell abundance: tissue-dependent 
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
) 
mod.view_anndata_setup()

print("Training")

mod.train(max_epochs=30000, 
          # train using full data (batch_size=None)
          batch_size=None, 
          # use all data points in training because 
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=False,
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);

# -----------------------------------------------------------------------------------------------
# Saving cell2location model
# -----------------------------------------------------------------------------------------------

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis,
    sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file
