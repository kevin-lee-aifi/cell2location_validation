#!bin/bash

conda create -y --prefix /home/jupyter/EXP-01058/envs/cell2loc_env python=3.9

conda activate "/home/jupyter/EXP-01058/envs/cell2loc_env"

conda install -c conda-forge mamba -y
conda install -c conda-forge ipykernel -y
conda install -c conda-forge pip -y

conda config --add channels bioconda

# cell2location and relevant depedencies
pip install cell2location[tutorial] jax==0.4.20 jaxlib==0.4.20 scvi-tools==1.0.4 scipy==1.11.1

pip install pyarrow

mamba install -c conda-forge rpy2 -y
mamba install -c conda-forge h5py -y

mamba install -c conda-forge r-base=4.3.1 -y
mamba install -c conda-forge r-seurat=4.1.1 -y
mamba install -c bioconda r-seurat-disk -y
mamba install -c bioconda r-azimuth -y
mamba install -c conda-forge r-hdf5r -y
mamba install -c conda-forge r-devtools -y
mamba install -c conda-forge r-tidyverse -y
mamba install -c conda-forge r-fast -y
mamba install -c bioconda r-anndata -y
