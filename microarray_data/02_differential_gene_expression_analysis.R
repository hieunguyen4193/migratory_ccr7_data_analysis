gc()
rm(list = ls())
library(dplyr)

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis/microarray_data"
source(file.path(path.to.project.src, "import_libraries.R"))

dataset_gse <- "GSE71170"

outdir <- "/media/hieunguyen/HNSD01/outdir"
PROJECT <- "FHager_datasets"

path.to.main.output <- file.path(outdir, PROJECT, output.version, "microarray_output", dataset_gse)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")

path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

thymus.rma.data <- readRDS(file.path(path.to.01.output, "thymus_rma_calibrated_data_object.rds"))
thymus.rma.exprs <- readRDS(file.path(path.to.01.output, "thymus_rma_calibrated_data_exprs_table.rds"))

lung.rma.data <- readRDS(file.path(path.to.01.output, "lung_rma_calibrated_data_object.rds"))
lung.rma.exprs <- readRDS(file.path(path.to.01.output, "lung_rma_calibrated_data_exprs_table.rds"))

library(tidyverse)
library(dplyr)

GPL6246 <- readRDS(file.path(path.to.01.output, "GPL6246_preprocessed.rds"))
GPL1261 <- readRDS(file.path(path.to.01.output, "GPL1261_preprocessed.rds"))

#####----------------------------------------------------------------------#####
# LUNG SAMPLES
#####----------------------------------------------------------------------#####
lung.design.matrix <- model.matrix(~ 0 + pData(lung.rma.data)$Status)

colnames(lung.design.matrix) <- c("Immature", "Mature")

# logFC = Immature / Mature
lung.contrast_matrix <- makeContrasts(Immature-Mature, levels = lung.design.matrix)

lung_fit <- eBayes(contrasts.fit(lmFit(lung.rma.data,
                                       design = lung.design.matrix), lung.contrast_matrix))

lung_DE_table <- topTable(lung_fit, number = Inf)

lung_DE_table <- lung_DE_table %>% 
  rownames_to_column("Gene.ID") %>% 
  rowwise %>% mutate(abs_log2FC = abs(logFC)) %>% 
  arrange(desc(abs_log2FC)) %>% mutate(Gene = paste(subset(GPL6246, GPL6246$PROBEID == Gene.ID)$SYMBOL, collapse = ", "))

#####----------------------------------------------------------------------#####
# THYMUS SAMPLES
#####----------------------------------------------------------------------#####
thymus.design.matrix <- model.matrix(~ 0 + pData(thymus.rma.data)$Status)

colnames(thymus.design.matrix) <- c("Immature", "Mature")

# logFC = Immature / Mature
thymus.contrast_matrix <- makeContrasts(Immature-Mature, levels = thymus.design.matrix) 

thymus_fit <- eBayes(contrasts.fit(lmFit(thymus.rma.data,
                                         design = thymus.design.matrix), thymus.contrast_matrix))

thymus_DE_table <- topTable(thymus_fit, number = Inf)

thymus_DE_table <- thymus_DE_table %>%
  rownames_to_column("Gene.ID") %>%
  rowwise %>% mutate(abs_log2FC = abs(logFC)) %>%
  arrange(desc(abs_log2FC)) %>% mutate(Gene = paste(subset(GPL1261, GPL1261$ID == Gene.ID)$Gene.Symbol, collapse = ", "))

