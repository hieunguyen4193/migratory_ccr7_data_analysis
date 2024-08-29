gc()
rm(list = ls())
my_random_seed <- 42
# __________VDJ DATA ANYLYSIS PIPELINE__________
path.to.main.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"
source(file.path(path.to.main.src, "config.R"))

if ("hdf5r" %in% installed.packages() == FALSE){
  install.packages("hdf5r")
}
config.version <- "v0.1"
output.version <- "20240828"
path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

set.seed(my_random_seed)

source(file.path(path2src, "import_libraries.R"))

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

save.RDS <- list(s1 = TRUE,
                 s2 = TRUE,
                 s3 = TRUE,
                 s4 = TRUE,
                 s5 = TRUE,
                 s6 = TRUE,
                 s7 = FALSE,
                 s8 = TRUE,
                 s8a = TRUE,
                 s9 = TRUE)

sw <- list(s1 = "on",
           s2 = "on",
           s3 = "on",
           s4 = "on",
           s5 = "on",
           s6 = "on",
           s7 = "off",
           s8 = "off",
           s8a = "on",
           s9 = "on")

rerun <- list(s1 = FALSE, 
              s2 = FALSE,
              s3 = FALSE,
              s4 = FALSE,
              s5 = FALSE,
              s6 = FALSE,
              s7 = FALSE,
              s8 = FALSE,
              s8a = FALSE,
              s9 = FALSE)

filter.thresholds <- config.params[[config.version]]

remove_doublet <- FALSE
num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25

# PROJECT <- "GSM5764245"
for (PROJECT in c("GSM5764259", 
                  "GSM5764245", 
                  "GSM5764288")){
  path.to.storage <- "/media/hieunguyen/HNSD01/storage"
  path.to.10X.doublet.estimation <- file.path(path.to.storage, "DoubletEstimation10X.csv")
  outdir <- "/media/hieunguyen/HNSD01/outdir"
  path.to.main.input <- file.path(path.to.storage, "FHager_datasets", "GSE192742_LIVER", PROJECT)
  path.to.main.output <- file.path(outdir, "FHager_datasets", output.version, sprintf("%s_%s", PROJECT, config.version))
  dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
  
  path2input <- path.to.main.input
  stage_lst <- c(PROJECT)
  names(stage_lst) <- c(PROJECT)
  
  MINCELLS  <- 0
  MINGENES  <- 0
  
  if (PROJECT == "GSM5764245"){
    input.method <- "readH5"  
  } else {
    input.method <- "from_txt_new"
  }
  
  print(PROJECT)
  print(stage_lst)
  
  s.obj <- run_pipeline_GEX(path2src=path2src,
                            path2input=path2input,
                            path.to.logfile.dir=file.path(path.to.main.output, "logs"),
                            stage_lst=stage_lst,
                            path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                            MINCELLS=MINCELLS,
                            MINGENES=MINGENES,
                            PROJECT=sprintf("%s_%s", PROJECT, config.version),
                            remove_doublet=remove_doublet,
                            save.RDS=save.RDS,
                            path.to.output=file.path(path.to.main.output),
                            rerun=rerun,
                            DE.test="wilcox",
                            num.PCA=num.PCA,
                            num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                            num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                            use.sctransform=FALSE,
                            filtered.barcodes=filtered.barcodes,
                            filter.thresholds=filter.thresholds,
                            path.to.anno.contigs=path.to.anno.contigs,
                            path.to.count.clonaltype=path.to.count.clonaltype,
                            input.method = input.method,
                            my_random_seed = my_random_seed,
                            sw = sw,
                            with.VDJ = FALSE)
}

##### After finishing the downstream analysis pipeline for each sample, we run the
##### integration. 

gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline"
source(file.path(path.to.pipeline.src, "/scRNA_GEX_pipeline/processes_src/s8_integration_and_clustering.R"))

outdir <- "/media/hieunguyen/HNSD01/outdir"

PROJECT <- "FHager_datasets"
output.version <- "20240828"
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

