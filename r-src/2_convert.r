remotes::install_github("theislab/zellkonverter")


packageVersion("Seurat")

packageVersion("zellkonverter")


setwd("/Users/stanleydale/user_generated/breault-lab/single-cell/dropbox-data/individualSeuratObjects")


d4 = readRDS("250525_egfDuodenoids_D4_AS_DZ.rds")


library(Seurat)


DefaultAssay(d4) <- "RNA"

# Run log-normalization
d4 <- NormalizeData(
  d4,
  normalization.method = "LogNormalize",  # default
  scale.factor = 10000                    # default library size scale
)




library(Seurat)                  # keep your v5 as-is
library(SingleCellExperiment)
library(zellkonverter)

sce <- as.SingleCellExperiment(d4)
gc()


class(assay(sce, "logcounts"))



rm(d4)

zellkonverter::writeH5AD(
  sce,
  "../../data/egfDuod_D4_AS_DZ.h5ad",
  X_name = "logcounts"
 # write 1000 cells at a time
)



rm(sce)




# D4_DZ -------------------------------------------------------------------




d4 = readRDS("250525_egfDuodenoids_D4_DZ.rds")


library(Seurat)


DefaultAssay(d4) <- "RNA"

# Run log-normalization
d4 <- NormalizeData(
  d4,
  normalization.method = "LogNormalize",  # default
  scale.factor = 10000                    # default library size scale
)




library(Seurat)                  # keep your v5 as-is
library(SingleCellExperiment)
library(zellkonverter)

sce <- as.SingleCellExperiment(d4)
gc()


class(assay(sce, "logcounts"))



rm(d4)

zellkonverter::writeH5AD(
  sce,
  "../../data/egfDuod_D4_DZ.h5ad",
  X_name = "logcounts"
  # write 1000 cells at a time
)






# D2 lapa -----------------------------------------------------------------




d2 = readRDS("250525_egfDuodenoids_D2_Lapa_DZ.rds")


library(Seurat)


DefaultAssay(d2) <- "RNA"

# Run log-normalization
d2 <- NormalizeData(
  d2,
  normalization.method = "LogNormalize",  # default
  scale.factor = 10000                    # default library size scale
)




library(Seurat)                  # keep your v5 as-is
library(SingleCellExperiment)
library(zellkonverter)

sce <- as.SingleCellExperiment(d2)
gc()


class(assay(sce, "logcounts"))



rm(d2)

zellkonverter::writeH5AD(
  sce,
  "../../data/egfDuod_D2_Lapa_DZ.h5ad",
  X_name = "logcounts"
  # write 1000 cells at a time
)




# D4 lapa dz --------------------------------------------------------------



d4 = readRDS("250525_egfDuodenoids_D4_Lapa_DZ.rds")


library(Seurat)


DefaultAssay(d4) <- "RNA"

# Run log-normalization
d4 <- NormalizeData(
  d4,
  normalization.method = "LogNormalize",  # default
  scale.factor = 10000                    # default library size scale
)




library(Seurat)                  # keep your v5 as-is
library(SingleCellExperiment)
library(zellkonverter)

sce <- as.SingleCellExperiment(d4)
gc()


class(assay(sce, "logcounts"))



rm(d4)

zellkonverter::writeH5AD(
  sce,
  "../../data/egfDuod_D4_Lapa_DZ.h5ad",
  X_name = "logcounts"
  # write 1000 cells at a time
)




# G6 dz -------------------------------------------------------------------




d6 = readRDS("250525_egfDuodenoids_G6_DZ.rds")


library(Seurat)


DefaultAssay(d6) <- "RNA"

# Run log-normalization
d6 <- NormalizeData(
  d6,
  normalization.method = "LogNormalize",  # default
  scale.factor = 10000                    # default library size scale
)




library(Seurat)                  # keep your v5 as-is
library(SingleCellExperiment)
library(zellkonverter)

sce <- as.SingleCellExperiment(d6)
gc()


class(assay(sce, "logcounts"))



rm(d6)

zellkonverter::writeH5AD(
  sce,
  "../../data/egfDuod_G6_DZ.h5ad",
  X_name = "logcounts"
  # write 1000 cells at a time
)




# D10 ---------------------------------------------------------------------




d10 = readRDS("fixedSeuratObjects/250530_egfDuodenoids_D10_Lapa_DZ_fixedV4.rds")


library(Seurat)


DefaultAssay(d10) <- "RNA"

# Run log-normalization
d10 <- NormalizeData(
  d10,
  normalization.method = "LogNormalize",  # default
  scale.factor = 10000                    # default library size scale
)




library(Seurat)                  # keep your v5 as-is
library(SingleCellExperiment)
library(zellkonverter)

sce <- as.SingleCellExperiment(d10)
gc()


class(assay(sce, "logcounts"))



rm(d10)

zellkonverter::writeH5AD(
  sce,
  "../../data/egfDuod_D10_Lapa_DZ.h5ad",
  X_name = "logcounts"
  # write 1000 cells at a time
)




