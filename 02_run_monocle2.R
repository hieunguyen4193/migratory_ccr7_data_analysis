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

outdir <- "/media/hieunguyen/HNSD01/outdir"

all.datasets <- c("220907_FH",
                  "GSM5764259",
                  "230228_FH",
                  "GSM5764288",            
                  "GSM5764245",
                  "integrate_GSE192742_LIVER",
                  "gutcellatlas_myeloid",
                  "220907_FH_cDC1",
                  "220907_FH_cDC2"
)

# to solve unable to access index for repository https://mran.microsoft.com/snapshot/2020-07-16/src/contrib
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)})

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

config.version <- "v0.1"
output.version <- "20240806"
PROJECT <- "FHager_datasets"

all.monocle.obj <- list()
for (orig.dataset in all.datasets){
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
  dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.01.output <- file.path(path.to.main.output, "01_output")
  path.to.02.output <- file.path(path.to.main.output, "02_output_monocle_2.12.0_v1")
  # path.to.02.output <- file.path(path.to.main.output, "02_output_monocle_2.12.0_v2")
  dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)
  path.to.monocle2.input <- file.path(path.to.main.output, "monocle2_inputs")
  
  if ("svglite" %in% installed.packages() == FALSE){
    install.packages("svglite")
  }
  
  library(monocle)
  
  if (file.exists(file.path(path.to.02.output, "monocle_obj.rds")) == FALSE){
    print(sprintf("running monocle2 on dataset %s", dataset.name))
    monocle.obj <- readRDS(file.path(path.to.monocle2.input, sprintf("%s.rds", dataset.name)))
    monocle.obj <- run_monocle2_from_presave_obj(monocle.obj, path.to.02.output)
  } else {
    print("monocle result exists, reading in...")
    monocle.obj <- readRDS(file.path(path.to.02.output, "monocle_obj.rds"))
  }
  all.monocle.obj[[dataset.name]] <- monocle.obj
  
  ##### plot cell trajectory, color by seurat clusters
  p <- plot_cell_trajectory(monocle.obj, color_by = "seurat_clusters")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.seurat_clsuters.svg", dataset.name), path = path.to.02.output, device = "svg", dpi = 300, width = 14, height = 10)
  
  ##### plot cell trajectory, color by monocle2 states
  p <- plot_cell_trajectory(monocle.obj, color_by = "State")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.State.svg", dataset.name), path = path.to.02.output, device = "svg", dpi = 300, width = 14, height = 10)
  
  ##### plot cell trajectory, color by pseudotime
  p <- plot_cell_trajectory(monocle.obj, color_by = "Pseudotime")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.pseudotime.svg", dataset.name), path = path.to.02.output, device = "svg", dpi = 300, width = 14, height = 10)
  
  ##### plot cell trajectory, color by pseudotime
  monocle.obj.reverse <- orderCells(monocle.obj, reverse = TRUE)
  p <- plot_cell_trajectory(monocle.obj.reverse, color_by = "Pseudotime")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.rev_Pseudotime.svg", dataset.name), path = path.to.02.output, device = "svg", dpi = 300, width = 14, height = 10)
  
  ##### save monocle data to csv file
  monocledf <- data.frame(
    barcode = colnames(monocle.obj),
    state = monocle.obj$State,
    pseudotime = monocle.obj$Pseudotime
  )
  monocle.reversedf <- data.frame(
    barcode = colnames(monocle.obj.reverse),
    state = monocle.obj.reverse$State,
    pseudotime = monocle.obj.reverse$Pseudotime
  )
  write.csv(monocledf, file.path(path.to.02.output, "monocledf.csv"))
  write.csv(monocle.reversedf, file.path(path.to.02.output, "monocledf.rev.csv"))
}


