setwd("/Users/stanleydale/user_generated/breault-lab/single-cell/dropbox-data/individualSeuratObjects")


d2 = readRDS("250525_egfDuodenoids_D2_DZ.rds")

library(Seurat)


# Show what metadata columns exist
colnames(d2@meta.data)

# Peek at first rows
head(d2@meta.data)

# Check the identities (clusters or assigned cell types)
head(Idents(d2))



## I dont think there are any cell type labels for these data




d4 = readRDS("250525_egfDuodenoids_D2_DZ.rds")

library(Seurat)


# Show what metadata columns exist
colnames(d4@meta.data)

# Peek at first rows
head(d4@meta.data)

# Check the identities (clusters or assigned cell types)
head(Idents(d4))
