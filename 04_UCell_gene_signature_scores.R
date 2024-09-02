gc()
rm(list = ls())

# install UCell
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("UCell", update = FALSE)


##### clean up #####
# gc()
# rm(list = ls())

#> This script must be run with Seurat 4 and monocle2 installed. Any monocle2 
#> version would work.
#> 
path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD01/outdir"

all.datasets <- c("220907_FH",
                  "GSM5764259", 
                  "integrate_GSE192742_LIVER",
                  "230228_FH",
                  "GSM5764288",            
                  "GSM5764245",
                  "gutcellatlas_myeloid",
                  "220907_FH_cDC1",
                  "220907_FH_cDC2")

all.s.obj <- list()
config.version <- "v0.1"
output.version <- "20240828"
PROJECT <- "FHager_datasets"
sigdf <- hash()
for (orig.dataset in all.datasets){
  print(sprintf("reading dataset %s in ",orig.dataset))
  input.config <- config.params[[config.version]]
  dataset.name <- sprintf("%s_%s", orig.dataset, config.version)
    if (orig.dataset == "integrate_GSE192742_LIVER"){
      umap.reduction.name <- "INTE_UMAP"
    } else {
      umap.reduction.name <- "RNA_UMAP"
    }
    path.to.main.output <- file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis")
    path.to.s.obj <- file.path(path.to.main.output, "01_output", sprintf("%s.rds", dataset.name))
    path.to.04.output <- file.path(path.to.main.output, "04_output")
    dir.create(file.path(path.to.04.output, "04_output"), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(path.to.04.output, sprintf("%s_added_UCell_module_score.rds", dataset.name))) == FALSE){
      
    path.to.monocle2.input <- file.path(path.to.main.output, "monocle2_inputs")
    dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)
    
    all.s.obj[[dataset.name]] <- readRDS(path.to.s.obj)
    
    cluster9.signature.genes <- readxl::read_excel(file.path(path.to.project.src, "Cluster 9 DEGs.xlsx")) %>% head(20) %>% pull("gene")
    if (orig.dataset == "gutcellatlas_myeloid"){
      signature.genes <- list(cluster9 = toupper(cluster9.signature.genes))
    } else {
      signature.genes <- list(cluster9 = cluster9.signature.genes)    
    }
    all.s.obj[[dataset.name]] <- UCell::AddModuleScore_UCell(all.s.obj[[dataset.name]], features = signature.genes)
    
    # FeaturePlot(object = all.s.obj[[dataset.name]], reduction = umap.reduction.name, features = "cluster9_UCell")
    feature.plot <- FeaturePlot(object = all.s.obj[[dataset.name]], features = "cluster9_UCell", label = TRUE, reduction = umap.reduction.name)
    violin.plot <- VlnPlot(object = all.s.obj[[dataset.name]], features = "cluster9_UCell")
    ggsave(plot = feature.plot, filename = sprintf("signatureScore_UMAP_%s.svg", dataset.name), 
           path = file.path(path.to.04.output), device = "svg", width = 14, height = 10, dpi = 300)
    ggsave(plot = violin.plot, filename = sprintf("signatureScore_violinplot_%s.svg", dataset.name), 
           path = file.path(path.to.04.output), device = "svg", width = 14, height = 10, dpi = 300)
    saveRDS(all.s.obj[[dataset.name]]@meta.data %>% rownames_to_column("barcode") %>% as.data.frame(), 
            file.path(path.to.04.output, sprintf("%s_added_UCell_module_score.rds", dataset.name)))
  } else {
    print(sprintf("Gene signature score for dataset %s finished", dataset.name))
    sigdf[[dataset.name]] <- readRDS( file.path(path.to.04.output, sprintf("%s_added_UCell_module_score.rds", dataset.name)))
  }
}


##### check genes in human dataset and compare to mouse dataset
s.obj.human <- readRDS(file.path(outdir, 
                                 PROJECT, 
                                 output.version, 
                                 "gutcellatlas_myeloid_v0.1", 
                                 "data_analysis", 
                                 "01_output", 
                                 sprintf("%s.rds", "gutcellatlas_myeloid_v0.1")))
s.obj.mouse <- readRDS(file.path(outdir, 
                                 PROJECT, 
                                 output.version, 
                                 "220907_FH_v0.1", 
                                 "data_analysis", 
                                 "01_output", 
                                 sprintf("%s.rds", "220907_FH_v0.1")))
cluster9.signature.genes <- readxl::read_excel(file.path(path.to.project.src, "Cluster 9 DEGs.xlsx")) %>% head(20) %>% pull("gene")
upper.cluster9.signature.genes <- toupper(cluster9.signature.genes)

intersect(upper.cluster9.signature.genes, row.names(s.obj.human))

setdiff(upper.cluster9.signature.genes, row.names(s.obj.human))

##### plot signature scores in different dataset
selected.clusters <- list(
  `220907_FH_v0.1` = 9,
  `integrate_GSE192742_LIVER_v0.1` = 19,
  `gutcellatlas_myeloid_v0.1` = 17
)

selected.sigdf <- data.frame()
for (i in names(selected.clusters)){
  tmpdf <- sigdf[[i]] %>%
    subset(select = c(cluster9_UCell, seurat_clusters)) %>%
    subset(seurat_clusters == selected.clusters[[i]])
  tmpdf$dataset.name <- i
  selected.sigdf <- rbind(selected.sigdf,  tmpdf)
}
path.to.main.output <- file.path(outdir, PROJECT, output.version)
p <- ggplot(data = selected.sigdf, aes(x = dataset.name, y = cluster9_UCell, fill = dataset.name)) + geom_violin()
ggsave(plot = p, filename = "selected_clusters_UCell_SignatureScores.svg", path = path.to.main.output, device = "svg", dpi = 300, width = 14, height = 10)
