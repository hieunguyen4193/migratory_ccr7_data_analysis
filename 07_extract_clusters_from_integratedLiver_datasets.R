##### clean up #####
gc()
rm(list = ls())
my_random_seed <- 42
set.seed(my_random_seed)

#> Monocle2 must be run on the docker tronghieunguyen/r36_monocle2
#> this docker contains the monocle version 2.12 installed in a R 3.6 env
#> which guarantees the reproducibility of the results
#> 
path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))

dataset.name <- "integrate_GSE192742_LIVER"

# outdir <- "/media/hieunguyen/HNSD01/outdir"
outdir <- "/media/hieunguyen/GSHD_HN01/outdir"

config.version <- "v0.1"
output.version <- "20240828"

dataset.name <- sprintf("%s_%s", dataset.name, config.version)
PROJECT <- "FHager_datasets"
path.to.main.output <- file.path(outdir, 
                                 PROJECT, 
                                 output.version, 
                                 dataset.name, 
                                 "data_analysis")
path.to.07.output <- file.path(path.to.main.output, "07_output")
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)


path.to.s.obj <- file.path(path.to.main.output, "01_output", sprintf("%s.rds", dataset.name))
s.obj <- readRDS(path.to.s.obj)
DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE)
selected.clusters <- c(2, 9, 19)

if (file.exists(file.path(path.to.07.output, sprintf("%s.subset_%s.rds", dataset.name, paste(selected.clusters, collapse = "_")))) == FALSE){
  path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline"
  source(file.path(path.to.pipeline.src, "/scRNA_GEX_pipeline/processes_src/s8_integration_and_clustering.R"))
  
  s.obj <- subset(s.obj, seurat_clusters %in% selected.clusters)
  
  num.dim.integration <- 30
  num.PCA <- 30
  num.PC.used.in.UMAP <- 30
  num.PC.used.in.Clustering <- 30
  cluster.resolution <- 0.5
  
  count.cell.in.samples <- table(s.obj$name)
  keep.samples <- names(count.cell.in.samples[count.cell.in.samples >= 30])
  s.obj <- subset(s.obj, name %in% keep.samples)
  s.obj <- NormalizeData(s.obj) # ---> use Log Normalized
  s.obj <- FindVariableFeatures(s.obj, selection.method = "vst")
  s.obj <- ScaleData(s.obj, features = rownames(s.obj))
  
  s.obj.integrated <- s8.integration.and.clustering(s.obj = s.obj, 
                                                    path.to.output = path.to.07.output, 
                                                    save.RDS.s8 = FALSE,
                                                    PROJECT = sprintf("%s.subset_%s.rds", dataset.name, paste(selected.clusters, collapse = "_")), 
                                                    num.dim.integration = num.dim.integration,
                                                    num.PCA = num.PCA,
                                                    num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                    num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                    cluster.resolution = cluster.resolution,
                                                    my_random_seed = 42,
                                                    umap.method = "uwot",
                                                    genes.to.not.run.PCA = NULL,
                                                    inte_pca_reduction_name = "INTE_PCA", 
                                                    inte_umap_reduction_name = "INTE_UMAP",
                                                    with.TSNE = FALSE,
                                                    k.filter = 200)
  
  saveRDS(s.obj, file.path(path.to.07.output, sprintf("%s.subset_%s.rds", dataset.name, paste(selected.clusters, collapse = "_"))))
}

