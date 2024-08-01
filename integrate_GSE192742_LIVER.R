gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline"
source(file.path(path.to.pipeline.src, "/scRNA_GEX_pipeline/processes_src/s8_integration_and_clustering.R"))

outdir <- "/media/hieunguyen/HNSD01/outdir"

PROJECT <- "FHager_datasets"
output.version <- "20240723"
config.version <- "v0.1"

path.to.integrate.output <- file.path(outdir, PROJECT, output.version, sprintf("integrate_GSE192742_LIVER_%s", config.version))
dir.create(path.to.integrate.output, showWarnings = FALSE, recursive = TRUE)

if (file.exists(file.path(path.to.integrate.output, "finished_saving_GSE192742.csv")) == FALSE){
  
  data.list <- list()
  for (sample.id in c("GSM5764245", "GSM5764259", "GSM5764288")){
    tmp <- readRDS(file.path(outdir, PROJECT, output.version, sprintf("%s_%s", sample.id, config.version), "s8a_output", sprintf("%s_%s.output.s8a.rds", sample.id, config.version)))
    tmp <- RenameCells(object = tmp, add.cell.id = sample.id)
    data.list[[sample.id]] <- tmp
  }
  
  s.obj <- merge(data.list[[1]], data.list[2:length(data.list)])
  
  num.dim.integration <- 30
  num.PCA <- 30
  num.PC.used.in.UMAP <- 30
  num.PC.used.in.Clustering <- 30
  cluster.resolution <- 0.5
  
  s.obj <- NormalizeData(s.obj) # ---> use Log Normalized
  s.obj <- FindVariableFeatures(s.obj, selection.method = "vst")
  s.obj <- ScaleData(s.obj, features = rownames(s.obj))
  
  s.obj.integrated <- s8.integration.and.clustering(s.obj = s.obj, 
                                                    path.to.output = path.to.integrate.output, 
                                                    save.RDS.s8 = TRUE,
                                                    PROJECT = sprintf("%s_%s", "integrate_GSE192742_LIVER", config.version), 
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
  
  write.csv(data.frame(status = c("finished integrating GSE192742")), file.path(path.to.integrate.output, "finished_saving_GSE192742.csv"))
}
