##### clean up #####
gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD01/outdir"

config.version <- "v0.1"
output.version <- "20240723"
PROJECT <- "FHager_datasets"

input.config <- config.params[[config.version]]

# path.to.save.html <- file.path(outdir, PROJECT, output.version, "html_outputs", "02_output")
# path.to.save.html <- file.path(outdir, PROJECT, output.version, "html_outputs", "02_output_monocle_2.32.0")
path.to.save.html <- file.path(outdir, PROJECT, output.version, "html_outputs", "02_output_monocle_2.32.0_v2")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### 01 html
#####----------------------------------------------------------------------#####
all.datasets <- c("220907_FH",
                  "GSM5764259", 
                  "integrate_GSE192742_LIVER",
                  "230228_FH",
                  "GSM5764288",            
                  "GSM5764245",
                  "gutcellatlas_myeloid"
)
for (orig.dataset in all.datasets){
  print(sprintf("working on the dataset %s", orig.dataset))
  path.to.rmd <- file.path(path.to.project.src, "02_trajectory_analysis_with_monocle2.Rmd")
  save.html.name <- sprintf("%s_%s.html", orig.dataset, config.version)
  if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
    rmarkdown::render(path.to.rmd, params = list(
      config.version = config.version,
      orig.dataset = orig.dataset,
      output.version = output.version
    ),
    output_file = save.html.name, 
    output_dir = path.to.save.html)
  } else {
    print(sprintf("The file %s already exists at %s", save.html.name, path.to.save.html))
  }
}
