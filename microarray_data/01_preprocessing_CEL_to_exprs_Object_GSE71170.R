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
  BiocManager::install("mogene10sttranscriptcluster.db", update = FALSE)
  install.packages("jpeg")
  BiocManager::install("arrayQualityMetrics", update = FALSE)
}

library(mogene10sttranscriptcluster.db)
source(file.path(path.to.project.src, "import_libraries.R"))
PROJECT <- "FHager_datasets"
output.version <- "20240806"

path.to.storage <- "/media/hieunguyen/HNSD01/storage"
path.to.main.input <- file.path(path.to.storage, PROJECT)
dataset_gse <- "GSE71170"

outdir <- "/media/hieunguyen/HNSD01/outdir"

path.to.main.output <- file.path(outdir, PROJECT, output.version, "microarray_output")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.output <- file.path(path.to.main.output, dataset_gse)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE)

path.to.dataset <- file.path(path.to.main.input, dataset_gse)
cel.datadir <- file.path(path.to.dataset, "CEL")

all.CEL.files <- Sys.glob(file.path(cel.datadir, "*.CEL"))

GPL1261.samples <- to_vec(for (item in all.CEL.files) if (grepl("MoGene", basename(item)) == FALSE & grepl("MOGENE", basename(item)) == FALSE) item)
GPL6246.samples <- to_vec(for (item in all.CEL.files) if (grepl("MoGene", basename(item)) == TRUE | grepl("MOGENE", basename(item)) == TRUE) item)

# oligo::read.celfiles: MUST INPUT SAME PLATFORM CEL FILES
# raw_data.GPL1261 <- oligo::read.celfiles(filenames = GPL1261.samples, verbose = FALSE)
# raw_data.GPL6246 <- oligo::read.celfiles(filenames = GPL6246.samples, verbose = FALSE)

thymus.immature.samples <- c("GSM1828776", "GSM1828777", "GSM1828778")
thymus.mature.samples <- c("GSM1828779", "GSM1828780", "GSM1828781")

lung.immature.samples <- c("GSM1828813", "GSM1828814", "GSM1828815")
lung.mature.samples <- c("GSM1828816", "GSM1828817")

thymus.samples <- c()
thymus.immature.celfiles <- c()
thymus.mature.celfiles <- c()
for (sample.id in c(thymus.immature.samples, thymus.mature.samples)){
  thymus.samples <- c(thymus.samples, to_vec(for (item in GPL1261.samples) if(grepl(sample.id, basename(item)) == TRUE) item)) 
  if (sample.id %in% thymus.immature.samples){
    thymus.immature.celfiles <- c(thymus.immature.celfiles, to_vec(for (item in GPL1261.samples) if(grepl(sample.id, basename(item)) == TRUE) item))  
  } else if (sample.id %in% thymus.mature.samples){
    thymus.mature.celfiles <- c(thymus.mature.celfiles, to_vec(for (item in GPL1261.samples) if(grepl(sample.id, basename(item)) == TRUE) item)) 
  }
}

lung.samples <- c()
lung.immature.celfiles <- c()
lung.mature.celfiles <- c()

for (sample.id in c(lung.immature.samples, lung.mature.samples)){
  lung.samples <- c(lung.samples, to_vec(for (item in GPL6246.samples) if(grepl(sample.id, basename(item)) == TRUE) item))  
  if (sample.id %in% lung.immature.samples == TRUE){
    lung.immature.celfiles <- c(lung.immature.celfiles, to_vec(for (item in GPL6246.samples) if(grepl(sample.id, basename(item)) == TRUE) item))  
  } 
  if (sample.id %in% lung.mature.samples == TRUE){
    lung.mature.celfiles <- c(lung.mature.celfiles, to_vec(for (item in GPL6246.samples) if(grepl(sample.id, basename(item)) == TRUE) item)) 
  }
}

sdrf_file1 <- read.delim(file.path(path.to.storage, PROJECT, dataset_gse, "E-GEOD-71170/E-GEOD-71170.sdrf.txt"))
sdrf_file1 <- subset(sdrf_file1, sdrf_file1$Array.Data.File %in% basename(lung.samples)) 
sdrf_file1 <- sdrf_file1 %>% rowwise %>% mutate(Status = ifelse(Array.Data.File %in% basename(lung.mature.celfiles), "Lung_Mature", "Lung_Immature")) %>% as.data.frame()
sdrf_file1$Status <- factor(sdrf_file1$Status, levels = c("Lung_Immature", "Lung_Mature"))
row.names(sdrf_file1) <- sdrf_file1$Array.Data.File
sdrf_file1 <- AnnotatedDataFrame(sdrf_file1)

sdrf_file2 <- read.delim(file.path(path.to.storage, PROJECT, dataset_gse, "E-GEOD-71166/E-GEOD-71166.sdrf.txt"))
sdrf_file2 <- subset(sdrf_file2, sdrf_file2$Array.Data.File %in% basename(thymus.samples)) 

sdrf_file2 <- sdrf_file2 %>% rowwise %>% mutate(Status = ifelse(Array.Data.File %in% basename(thymus.mature.celfiles), "Thymus_Mature", "Thymus_Immature")) %>% as.data.frame()
sdrf_file2$Status <- factor(sdrf_file2$Status, levels = c("Thymus_Immature", "Thymus_Mature"))
row.names(sdrf_file2) <- sdrf_file2$Array.Data.File
sdrf_file2 <- AnnotatedDataFrame(sdrf_file2)


lung.raw_data <- oligo::read.celfiles(filenames = file.path(path.to.dataset, "CEL", row.names(sdrf_file1)), 
                                      verbose = FALSE, 
                                      phenoData = sdrf_file1)

thymus.raw_data <- oligo::read.celfiles(filenames = file.path(path.to.dataset, "CEL", row.names(sdrf_file2)), 
                                        verbose = FALSE, 
                                        phenoData = sdrf_file2)

thymus.rma.data <- oligo::rma(thymus.raw_data, normalize = TRUE)
thymus.rma.exprs <- Biobase::exprs(thymus.rma.data)

lung.rma.data <- oligo::rma(lung.raw_data, normalize = TRUE)
lung.rma.exprs <- Biobase::exprs(lung.rma.data)

path.to.platform.info.dir <- file.path(path.to.storage, PROJECT, "GPL")

##### CONVERT PROBE IDS OF THE PLATFORM TO GENE.SYMBOL
gpl1261 <- read.delim(file.path(path.to.platform.info.dir, "GPL1261-56135.filtered.txt"))
gpl1261 <- subset(gpl1261, select = c(ID, Gene.Symbol))

annotation.lung.rma.data <- AnnotationDbi::select(mogene10sttranscriptcluster.db,
                                                  keys = (featureNames(lung.rma.data)),
                                                  columns = c("SYMBOL", "GENENAME"),
                                                  keytype = "PROBEID")

annotation.lung.rma.data <- subset(annotation.lung.rma.data, !is.na(SYMBOL))

#####----------------------------------------------------------------------#####
# Save objects to 01_output
saveRDS(thymus.rma.data, file.path(path.to.01.output, "thymus_rma_calibrated_data_object.rds"))
saveRDS(thymus.rma.exprs, file.path(path.to.01.output, "thymus_rma_calibrated_data_exprs_table.rds"))

saveRDS(lung.rma.data, file.path(path.to.01.output, "lung_rma_calibrated_data_object.rds"))
saveRDS(lung.rma.exprs, file.path(path.to.01.output, "lung_rma_calibrated_data_exprs_table.rds"))

# convert probe ID to Gene symbol
saveRDS(annotation.lung.rma.data, file.path(path.to.01.output, "GPL6246_preprocessed.rds"))
saveRDS(gpl1261, file.path(path.to.01.output, "GPL1261_preprocessed.rds"))


