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

# orig.dataset <- "220907_FH"
for (orig.dataset in all.datasets){
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
  } else {
    path.to.main.input <- file.path(outdir,
                                    PROJECT,
                                    output.version, 
                                    dataset.name, 
                                    "s8a_output",
                                    sprintf("%s.output.s8a.rds", dataset.name))
  }
  
  
  path.to.main.output <- file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis")
  path.to.monocle2.input <- file.path(path.to.main.output, "monocle2_inputs")
  dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)
  
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
}

