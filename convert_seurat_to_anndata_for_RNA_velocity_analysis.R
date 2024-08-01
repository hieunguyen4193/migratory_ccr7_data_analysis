##### clean up #####
gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))

outdir <- "/media/hieunguyen/HNSD01/outdir"
orig.dataset <- "GSM5764245"
config.version <- "v0.1"
output.version <- "20240723"
PROJECT <- "FHager_datasets"

input.config <- config.params[[config.version]]
dataset.name <- sprintf("%s_%s", orig.dataset, config.version)

path.to.main.input <- file.path(outdir,
                                PROJECT,
                                output.version, 
                                dataset.name, 
                                "s8a_output",
                                sprintf("%s.output.s8a.rds", dataset.name))

path.to.seurat2anndata <- file.path(outdir, PROJECT, output.version, "seurat2anndata", dataset.name)
dir.create(path.to.seurat2anndata, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(path.to.main.input)

object.name <- dataset.name

s.obj$barcode <- colnames(s.obj)

if (dataset.name = sprintf("integrate_GSE192742_LIVER_%s", config.version)){
  s.obj$UMAP_1 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,1]
  s.obj$UMAP_2 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,2]
} else {
  s.obj$UMAP_1 <- s.obj@reductions$RNA_UMAP@cell.embeddings[,1]
  s.obj$UMAP_2 <- s.obj@reductions$RNA_UMAP@cell.embeddings[,2]
}

# write dimesnionality reduction matrix, in this example s.obj.case pca matrix
write.csv(s.obj@reductions$INTE_UMAP@cell.embeddings, file=file.path(path.to.seurat2anndata, sprintf('pca_%s.csv', object.name)), quote=F, row.names=F)
write.csv(s.obj@meta.data, file=file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', object.name)), quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(s.obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=file.path(path.to.seurat2anndata, sprintf('counts_%s.mtx', object.name)))

# write gene names
write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(path.to.seurat2anndata, sprintf('gene_names_%s.csv', object.name)),
             quote=F,row.names=F,col.names=F)
