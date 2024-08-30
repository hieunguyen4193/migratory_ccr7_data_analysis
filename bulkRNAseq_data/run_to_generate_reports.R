gc()
rm(list=ls())

outdir <- "/media/hieunguyen/HNSD01/outdir"
PROJECT <- "FHager_datasets"
data.version <- "20240828"

path.to.main.output <- file.path(outdir, PROJECT, data.version, "bulkRNAseq_data")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.save.html.output <- file.path(path.to.main.output, "html_output")
dir.create(path.to.save.html.output, showWarnings = FALSE, recursive = TRUE)

rmd.dir <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis/bulkRNAseq_data"

render_rmd_report <- function( 
    cell.type, 
    sample.type1, 
    sample.type2, 
    cutoff.logFC, 
    cutoff.adjp, 
    template.name = "01_analysis_report.Rmd",
    html.dir = path.to.save.html.output){
  dir.create(html.dir, showWarnings = F)
  rmarkdown::render(file.path(rmd.dir, template.name), 
                    params = list(
                      cell.type = cell.type,
                      sample.type1 = sample.type1,
                      sample.type2 = sample.type2,
                      cutoff.adjp = cutoff.adjp,
                      cutoff.logFC = cutoff.logFC
                    ),
                    output_file = sprintf("01_main_analysis_report_GSE160156_DESeq2_%s_%s_vs_%s_cutoff_%s", cell.type, sample.type1, sample.type2, cutoff.adjp),
                    output_dir = html.dir,
  )
}

#####-------------------------------------------------------------------------#####
##### Cut-off 0.01 on adjusted p value
#####-------------------------------------------------------------------------#####
render_rmd_report("CD11b", "LP", "mL", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)

render_rmd_report("CD11b", "LP", "mLN", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("CD11b", "mL", "mLN", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("CD103", "LP", "mL", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("CD103", "LP", "mLN", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("CD103", "mL", "mLN", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("DP", "LP", "mL", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("DP", "LP", "mLN", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("DP", "mL", "mLN", 1, 0.01, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)

#####-------------------------------------------------------------------------#####
##### Cut-off 0.05 on adjusted p value
#####-------------------------------------------------------------------------#####
render_rmd_report("CD11b", "LP", "mL", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)

render_rmd_report("CD11b", "LP", "mLN", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("CD11b", "mL", "mLN", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("CD103", "LP", "mL", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("CD103", "LP", "mLN", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("CD103", "mL", "mLN", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("DP", "LP", "mL", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("DP", "LP", "mLN", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)
render_rmd_report("DP", "mL", "mLN", 1, 0.05, 
                  template.name = "01_analysis_report.Rmd",
                  html.dir = path.to.save.html.output
)

