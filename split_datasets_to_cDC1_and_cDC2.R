gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))
source(file.path(path.to.project.src, "list_cluster_cDC1_cDC2.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####


outdir <- "/media/hieunguyen/HNSD01/outdir"
PROJECT <- "FHager_datasets"
org.dataset <- "220907_FH"
config.version <- "v0.1"
output.version <- "20240806"
  
input.config <- config.params[[config.version]]
dataset.name <- sprintf("%s_%s", org.dataset, config.version)

if (file.exists(file.path(outdir,
                          PROJECT,
                          output.version, "finished_splitting_cDC1_cDC2.csv")) == FALSE){
  ##### read the saved object in 01_output
  s.obj <- readRDS(file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis", "01_output", sprintf("%s.rds", dataset.name)))
  
  cdc1.clusters <- selected.clusters[[org.dataset]][[config.version]]$cDC1
  cdc2.clusters <- selected.clusters[[org.dataset]][[config.version]]$cDC2
  
  s.obj.cdc1 <- subset(s.obj, seurat_clusters %in% cdc1.clusters)
  s.obj.cdc2 <- subset(s.obj, seurat_clusters %in% cdc2.clusters)
  
  re_process_subset_dataset <- function(input.s.obj){
    num.PCA <- 25
    num.PC.used.in.UMAP <- 25
    num.PC.used.in.Clustering <- 25
    chosen.assay <- "RNA"
    my_random_seed <- 42
    cluster.resolution <- 1
    
    input.s.obj <- DietSeurat(input.s.obj)
    input.s.obj <- NormalizeData(input.s.obj) # ---> use Log Normalized
    input.s.obj <- FindVariableFeatures(input.s.obj, selection.method = "vst")
    input.s.obj <- ScaleData(input.s.obj, features = rownames(input.s.obj))
    
    DefaultAssay(input.s.obj) <- chosen.assay
    
    input.s.obj <- RunPCA(input.s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
    
    input.s.obj <- RunUMAP(input.s.obj, reduction = sprintf("%s_PCA", chosen.assay), 
                           dims = 1:num.PC.used.in.UMAP, reduction.name=sprintf("%s_UMAP", chosen.assay),
                           seed.use = my_random_seed, umap.method = "uwot")
    
    # clustering 
    input.s.obj <- FindNeighbors(input.s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
    input.s.obj <- FindClusters(input.s.obj, resolution = 0.5, random.seed = 0)
    return(input.s.obj)
  }
  s.obj.cdc1 <- re_process_subset_dataset(s.obj.cdc1)
  s.obj.cdc2 <- re_process_subset_dataset(s.obj.cdc2)
  
  #### save the cDC1 and cDC2 objects to the same folder structure as the main datasets
  dir.create(file.path(outdir,
             PROJECT,
             output.version, 
             sprintf("%s_%s_%s", org.dataset, "cDC1", config.version), 
             "s8a_output"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(outdir,
             PROJECT,
             output.version, 
             sprintf("%s_%s_%s", org.dataset, "cDC2", config.version), 
             "s8a_output"), showWarnings = FALSE, recursive = TRUE)
  
  saveRDS(s.obj.cdc1, file.path(outdir,
                                PROJECT,
                                output.version, 
                                sprintf("%s_%s_%s", org.dataset, "cDC1", config.version), 
                                "s8a_output",
                                sprintf("%s_%s_%s.output.s8a.rds", org.dataset, "cDC1", config.version)))
  saveRDS(s.obj.cdc2, file.path(outdir,
                                PROJECT,
                                output.version, 
                                sprintf("%s_%s_%s", org.dataset, "cDC2", config.version), 
                                "s8a_output",
                                sprintf("%s_%s_%s.output.s8a.rds", org.dataset, "cDC2", config.version)))
  
  write.csv(data.frame(status = c("finished")), file.path(outdir,
                                                          PROJECT,
                                                          output.version, "finished_splitting_cDC1_cDC2.csv"))
} else {
  print(file.path(outdir,
                  PROJECT,
                  output.version, "finished_splitting_cDC1_cDC2.csv"))
}