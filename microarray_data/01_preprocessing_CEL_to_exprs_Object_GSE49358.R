#####----------------------------------------------------------------------#####
# Raw .CEL data are downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71171
# Download SDRF and IDF data files from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-71170#
# GSE71171 cannot not be queried from by the GEOquery function

# Similarly, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49358
# and https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-49358?query=GSE49358
# GSE49358 can be queried by the GEOquery function
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis/microarray_data"

if ("mogene10sttranscriptcluster.db" %in% installed.packages() == FALSE){
  BiocManager::install("mogene10sttranscriptcluster.db")
}

library(mogene10sttranscriptcluster.db)
source(file.path(path.to.project.src, "import_libraries.R"))
PROJECT <- "FHager_datasets"
path.to.storage <- "/media/hieunguyen/HNSD01/storage"
path.to.main.input <- file.path(path.to.storage, PROJECT)

dataset_gse <- "GSE49358"

path.to.main.output <- file.path(outdir, PROJECT, "microarray_output")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.output <- file.path(path.to.main.output, dataset_gse)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE)

gset <- getGEO(dataset_gse, GSEMatrix =TRUE, getGPL=FALSE)

gset.exprs <- Biobase::exprs(gset)
