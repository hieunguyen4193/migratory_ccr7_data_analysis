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

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
# outdir <- "/media/hieunguyen/HD0/outdir/CRC1382"

outdir <- "/media/hieunguyen/HNSD01/outdir"

all.datasets <- c("220907_FH",
                  "GSM5764259", 
                  "integrate_GSE192742_LIVER",
                  "230228_FH",
                  "GSM5764288",            
                  "GSM5764245",
                  "gutcellatlas_myeloid")

all.s.obj <- list()

for (orig.dataset in all.datasets){
  print(sprintf("reading dataset %s in ",orig.dataset))
  config.version <- "v0.1"
  output.version <- "20240723"
  
  PROJECT <- "FHager_datasets"
  
  input.config <- config.params[[config.version]]
  dataset.name <- sprintf("%s_%s", orig.dataset, config.version)
  
  if (orig.dataset == "integrate_GSE192742_LIVER"){
    path.to.main.input <- file.path(outdir,
                                    PROJECT,
                                    output.version, 
                                    dataset.name, 
                                    "s8_output",
                                    sprintf("%s.output.s8.rds", dataset.name))
    umap.reduction.name <- "INTE_UMAP"
  } else {
    path.to.main.input <- file.path(outdir,
                                    PROJECT,
                                    output.version, 
                                    dataset.name, 
                                    "s8a_output",
                                    sprintf("%s.output.s8a.rds", dataset.name))
    umap.reduction.name <- "RNA_UMAP"
  }
  
  
  path.to.main.output <- file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis")
  path.to.04.output <- file.path(path.to.main.output, "04_output")
  dir.create(file.path(path.to.04.output, "04_output"), showWarnings = FALSE, recursive = TRUE)
  
  path.to.monocle2.input <- file.path(path.to.main.output, "monocle2_inputs")
  dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)
  
  all.s.obj[[dataset.name]] <- readRDS(path.to.main.input)
  
  cluster9.signature.genes <- readxl::read_excel(file.path(path.to.project.src, "Cluster 9 DEGs.xlsx")) %>% head(20) %>% pull("gene")
  signature.genes <- list(cluster9 = cluster9.signature.genes)
  all.s.obj[[dataset.name]] <- UCell::AddModuleScore_UCell(all.s.obj[[dataset.name]], features = signature.genes)
  
  FeaturePlot(object = all.s.obj[[dataset.name]], reduction = umap.reduction.name, features = "cluster9_UCell")
  saveRDS(all.s.obj[[dataset.name]]@meta.data %>% rownames_to_column("barcode") %>% as.data.frame(), 
          file.path(path.to.04.output, sprintf("%s_added_UCell_module_score.csv", dataset.name)))
}


