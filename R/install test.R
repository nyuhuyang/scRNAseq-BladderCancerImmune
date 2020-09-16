conda create --name r4.0.2
conda activate r4.0.2
conda install -c conda-forge r-base=4.0.2

BiocManager::install("SingleR")
BiocManager::install("SingleCellExperiment")


conda install -c bioconda scanpy bbknn anndata2ri bioconductor-singler bioconductor-scrnaseq bioconductor-singlecellexperiment bioconductor-mast bioconductor-monocle gseapy bioconductor-fgsea r-enrichr r-openxlsx r-harmony r-seurat #bioconductor-dropletutils bioconductor-scater


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")


pip install harmonyTS
remotes::install_github("mojaveazure/seurat-disk")

#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
