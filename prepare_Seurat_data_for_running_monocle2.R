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

config.version <- "v0.1"
output.version <- "20240828"
PROJECT <- "FHager_datasets"

for (orig.dataset in all.datasets){
  input.config <- config.params[[config.version]]
  dataset.name <- sprintf("%s_%s", orig.dataset, config.version)
  path.to.main.input <- file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis", "01_output", sprintf("%s.rds", dataset.name))
  path.to.main.output <- file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis")
  path.to.monocle2.input <- file.path(path.to.main.output, "monocle2_inputs")
  dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.monocle2.input, sprintf("%s.rds", dataset.name))) == FALSE){
    print(sprintf("working on %s", dataset.name))
    s.obj <- readRDS(path.to.main.input)
    
    data <- GetAssayData(s.obj, slot = "count", assay = "RNA")
    
    pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
    
    fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fd)
    
    library(monocle)
    monocle.obj <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())
    saveRDS(monocle.obj, file.path(path.to.monocle2.input, sprintf("%s.rds", dataset.name)))
  } else {
    print(sprintf("File %s existed", file.path(path.to.monocle2.input, sprintf("%s.rds", dataset.name))))
  }
}

