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

dataset.name <- "gutcellatlas_myeloid"

outdir <- "/media/hieunguyen/HNSD01/outdir"
config.version <- "v0.1"
output.version <- "20240828"

dataset.name <- sprintf("%s_%s", dataset.name, config.version)
PROJECT <- "FHager_datasets"
path.to.main.output <- file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis")
path.to.06.output <- file.path(path.to.main.output, "06_output")
dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)

path.to.s.obj <- file.path(path.to.main.output, "01_output", sprintf("%s.rds", dataset.name))
s.obj <- readRDS(path.to.s.obj)
DimPlot(object = s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE)
s.obj <- subset(s.obj, seurat_clusters %in% c(15, 17, 22))

num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25

pca_reduction_name <- "RNA_PCA"
umap_reduction_name <- "RNA_UMAP"

s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=pca_reduction_name)
s.obj <- RunUMAP(s.obj, reduction = pca_reduction_name, 
                 dims = 1:num.PC.used.in.UMAP, reduction.name=umap_reduction_name,
                 seed.use = my_random_seed, umap.method = "uwot")
# clustering 
s.obj <- FindNeighbors(s.obj, reduction = pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = 0)