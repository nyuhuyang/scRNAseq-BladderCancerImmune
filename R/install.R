conda create --name r4.0
conda activate r4.0
conda install -c conda-forge r-base rpy2 python matplotlib numpy pandas scipy anndata r-dplyr r-qlcmatrix r-tidyr r-matrix xlrd notebook openssl pyopenssl python-igraph r-cellranger r-data.table r-ggpubr h5py hdf5 ipykernel ipython jpeg jupyterlab r-tidyverse r-matrix r-reshape2 r-devtools leidenalg r-leiden
conda install -c bioconda r-seurat scanpy bbknn anndata2ri bioconductor-dropletutils bioconductor-scater bioconductor-singlecellexperiment

devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
