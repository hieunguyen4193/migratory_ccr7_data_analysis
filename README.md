# Data analysis 

## Pipeline

- Run the function `run_pipeline.R`: Run the usual downstream analysis pipeline using `Seurat` to process the dataset `220907_FH` and `230228_FH`. 

- Run the function `run_pipeline_GSE192742.R`: same as above, downstream analysis pipeline for the 3 samples `GSM5764259`, `GSM5764245`, `GSM5764288` from the dataset GSE192742. 

- Similarly, run the function `run_pipeline_GutCellAtlas_Myeloid.R` to analyze the dataset download from Gut Cell Atlas, myeloid cells only. 

- For `GSM5764259`, `GSM5764245`, `GSM5764288`, download SRA raw data from SRA portal using the `sra-toolkit`. Convert SRA file to FASTQ file using `fastq-dump`. Run `CellRanger` pipeline. Once we have the output from 
`CellRanger`, generate spliced and unspliced data fro RNA velocity analysis using `velocyto`, see `run_velocty.sh`. Install `velocyto` in a `conda` environment. 

## Data analysis pipeline
- Run `run_pipeline.R`: Run the downstream analysis pipeline (`Seurat`) for the dataset `220907_FH` and `230228_FH`.

- Run `run_pipeline_GSE192742_LIVER.R`: Run the downstream analysis pipeline for the dataset `GSE192742` including the three samples `GSM5764259`, `GSM5764245`, `GSM5764288`.

- Run `run_to_generate_html_01.R`

- Run `prepare_Seurat_data_for_running_monocle2.R`

- Run `generate_color_lists.R`



 