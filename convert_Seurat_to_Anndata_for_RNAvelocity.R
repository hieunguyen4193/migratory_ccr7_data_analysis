##### clean up #####
gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))

all.datasets <- c("220907_FH",
                  "GSM5764259", 
                  "integrate_GSE192742_LIVER",
                  "230228_FH",
                  "GSM5764288",            
                  "GSM5764245",
                  "gutcellatlas_myeloid",
                  "220907_FH_cDC1",
                  "220907_FH_cDC2"
)

outdir <- "/media/hieunguyen/HNSD01/outdir"

config.version <- "v0.1"
output.version <- "20240806"
PROJECT <- "FHager_datasets"

input.config <- config.params[[config.version]]

for (orig.dataset in all.datasets){
  print(sprintf("working on dataset %s", orig.dataset))
  dataset.name <- sprintf("%s_%s", orig.dataset, config.version)
  path.to.seurat2anndata <- file.path(outdir, PROJECT, output.version, "seurat2anndata", dataset.name)
  dir.create(path.to.seurat2anndata, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(path.to.seurat2anndata, sprintf('finished_%s.csv', dataset.name))) == FALSE){
    if (orig.dataset == "integrate_GSE192742_LIVER"){
      path.to.main.input <- file.path(outdir,
                                      PROJECT,
                                      output.version, 
                                      dataset.name, 
                                      "s8_output",
                                      sprintf("%s.output.s8.rds", dataset.name))
    } else if (orig.dataset == "220907_FH") {
      path.to.main.input <- "/media/hieunguyen/HNSD01/outdir/FHager_datasets/20240806/220907_FH_v0.1/data_analysis/01_output/220907_FH_v0.1.rds"
    } else {
      path.to.main.input <- file.path(outdir,
                                      PROJECT,
                                      output.version, 
                                      dataset.name, 
                                      "s8a_output",
                                      sprintf("%s.output.s8a.rds", dataset.name))
    }
    
    
    
    
    s.obj <- readRDS(path.to.main.input)
    
    s.obj$barcode <- colnames(s.obj)
    
    if (dataset.name == sprintf("integrate_GSE192742_LIVER_%s", config.version)){
      s.obj$UMAP_1 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,1]
      s.obj$UMAP_2 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,2]
      write.csv(s.obj@reductions$INTE_UMAP@cell.embeddings, file=file.path(path.to.seurat2anndata, sprintf('pca_%s.csv', dataset.name)), quote=F, row.names=F)
      write.csv(s.obj@meta.data, file=file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', dataset.name)), quote=F, row.names=F)
    } else {
      s.obj$UMAP_1 <- s.obj@reductions$RNA_UMAP@cell.embeddings[,1]
      s.obj$UMAP_2 <- s.obj@reductions$RNA_UMAP@cell.embeddings[,2]
      write.csv(s.obj@reductions$RNA_UMAP@cell.embeddings, file=file.path(path.to.seurat2anndata, sprintf('pca_%s.csv', dataset.name)), quote=F, row.names=F)
      write.csv(s.obj@meta.data, file=file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', dataset.name)), quote=F, row.names=F)
    }
    
    # write expression counts matrix
    library(Matrix)
    counts_matrix <- GetAssayData(s.obj, assay='RNA', slot='data')
    writeMM(counts_matrix, file=file.path(path.to.seurat2anndata, sprintf('counts_%s.mtx', dataset.name)))
    
    # write gene names
    write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(path.to.seurat2anndata, sprintf('gene_names_%s.csv', dataset.name)),
                 quote=F,row.names=F,col.names=F)
    write.csv(data.frame(status = c("finished converting for the dataset")), file.path(path.to.seurat2anndata, sprintf('finished_%s.csv', dataset.name)))
  } else {
    print(sprintf("File %s exists", file.path(path.to.seurat2anndata, sprintf('finished_%s.csv', dataset.name))))
  }
}

