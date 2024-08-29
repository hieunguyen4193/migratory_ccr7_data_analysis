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
  path.to.save.colors <- file.path(outdir, PROJECT, output.version, "colors", dataset.name)
  if (file.exists(file.path(path.to.save.colors, "colordf.csv")) == FALSE){
    dir.create(path.to.save.colors, showWarnings = FALSE, recursive = TRUE)
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
    s.obj <- readRDS(path.to.main.input)
    num.clusters <- length(unique(s.obj$seurat_clusters))
    library(scales)
    colors <- hue_pal()(num.clusters)
    write.csv(data.frame(color = colors), file.path(path.to.save.colors, "colordf.csv"))
  } else {
    print(sprintf("File %s existed", file.path(path.to.save.colors, "colordf.csv")))
  }
}