gc()
rm(list = ls())
my_random_seed <- 42
# __________VDJ DATA ANYLYSIS PIPELINE__________
path.to.main.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"
source(file.path(path.to.main.src, "config.R"))

all.PROJECTS <- list( `220907_Kopplin_Pabst` = "220907_FH",
                      `230228_Kopplin_Pabst` = "230228_FH")

path.to.storage <- "/media/hieunguyen/HNSD01/storage"
output.version <- "20240828"
outdir <- "/media/hieunguyen/HNSD01/outdir"

if ("hdf5r" %in% installed.packages() == FALSE){
  install.packages("hdf5r")
}

for (config.version in c("v0.1")){
  for (dataset.name in names(all.PROJECTS)){
    PROJECT <- sprintf("%s_%s", all.PROJECTS[[dataset.name]], config.version)

    path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline"
    path2src <- file.path(path.to.pipeline.src, "processes_src")

    set.seed(my_random_seed)

    source(file.path(path2src, "import_libraries.R"))

    source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

    # _____ROUND_____
    analysis.round <- "1st" # <<<<< CHOSE THE ROUND OF DATA ANALYSIS HERE
    path2input <- file.path(path.to.storage, "FHager_datasets", sprintf("%s_with_EGFP_SnapGene/%s", dataset.name, dataset.name))

    stage_lst <- c("with_SnapGene_EGFP")
    names(stage_lst) <- c(dataset.name)

    MINCELLS  <- 0
    MINGENES  <- 0

    save.RDS <- list(s1 = TRUE,
                     s2 = TRUE,
                     s3 = TRUE,
                     s4 = TRUE,
                     s5 = TRUE,
                     s6 = TRUE,
                     s7 = FALSE,
                     s8 = TRUE,
                     s8a = TRUE,
                     s9 = TRUE)

    sw <- list(s1 = "on",
               s2 = "on",
               s3 = "on",
               s4 = "on",
               s5 = "on",
               s6 = "on",
               s7 = "off",
               s8 = "off",
               s8a = "on",
               s9 = "on")

    rerun <- list(s1 = FALSE,
                  s2 = FALSE,
                  s3 = FALSE,
                  s4 = FALSE,
                  s5 = FALSE,
                  s6 = FALSE,
                  s7 = FALSE,
                  s8 = FALSE,
                  s8a = FALSE,
                  s9 = FALSE)

    filter.thresholds <- config.params[[config.version]]

    remove_doublet <- FALSE
    path.to.10X.doublet.estimation <- file.path(path.to.storage, "DoubletEstimation10X.csv")

    num.PCA <- 25
    num.PC.used.in.UMAP <- 25
    num.PC.used.in.Clustering <- 25

    for (sample.id in names(stage_lst)){
      path.to.output <- file.path(outdir, "FHager_datasets", output.version, PROJECT)
      dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)

      s.obj <- run_pipeline_GEX(path2src=path2src,
                                path2input=path2input,
                                path.to.logfile.dir=file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round), "logs"),
                                stage_lst=stage_lst,
                                path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                                MINCELLS=MINCELLS,
                                MINGENES=MINGENES,
                                PROJECT=PROJECT,
                                remove_doublet=remove_doublet,
                                save.RDS=save.RDS,
                                path.to.output=file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round)),
                                rerun=rerun,
                                DE.test="wilcox",
                                num.PCA=num.PCA,
                                num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                                num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                                use.sctransform=FALSE,
                                filtered.barcodes=filtered.barcodes,
                                filter.thresholds=filter.thresholds,
                                path.to.anno.contigs=path.to.anno.contigs,
                                path.to.count.clonaltype=path.to.count.clonaltype,
                                input.method = "normal",
                                my_random_seed = my_random_seed,
                                sw = sw,
                                with.VDJ = FALSE)
    }
    writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))
  }
}

# EOF