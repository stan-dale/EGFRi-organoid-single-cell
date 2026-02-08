install.packages("remotes")
remotes::install_version("SeuratObject",  version = "5.0.2", repos = "https://cloud.r-project.org")
remotes::install_version("Seurat",        version = "5.0.3", repos = "https://cloud.r-project.org")


packageVersion("Seurat")


setwd("/Users/stanleydale/user_generated/breault-lab/single-cell/dropbox-data/individualSeuratObjects")


d2 = readRDS("250525_egfDuodenoids_D2_DZ.rds")


library(Seurat)


DefaultAssay(d2) <- "RNA"

# Run log-normalization
d2 <- NormalizeData(
  d2,
  normalization.method = "LogNormalize",  # default
  scale.factor = 10000                    # default library size scale
)






install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

library(Seurat)                  # keep your v5 as-is
library(SingleCellExperiment)
library(zellkonverter)

sce <- as.SingleCellExperiment(d2)
zellkonverter::writeH5AD(
  sce,
  "../../data/egfDuod_D2_DZ.h5ad",
  X_name = "logcounts"          # or "counts" if you want raw in X
)