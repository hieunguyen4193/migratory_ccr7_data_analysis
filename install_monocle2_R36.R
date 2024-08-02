# to solve unable to access index for repository https://mran.microsoft.com/snapshot/2020-07-16/src/contrib
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)})

install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.9", ask = FALSE)

# https://cran.r-project.org/src/contrib/Archive/VGAM/VGAM_1.0-6.tar.gz
# https://cran.r-project.org/src/contrib/Archive/fastICA/fastICA_1.2-1.tar.gz
# https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.2-6.tar.gz

path.to.main.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

install.packages(file.path(path.to.main.src, "VGAM_1.0-6.tar.gz"), type = "source", repos = NULL)
install.packages(file.path(path.to.main.src, "Matrix_1.3-0.tar.gz"), type = "source", repos = NULL)
install.packages(file.path(path.to.main.src, "fastICA_1.2-1.tar.gz"), type = "source", repos = NULL)
install.packages(file.path(path.to.main.src, "FNN_1.0.tar.gz"), type = "source", repos = NULL)
install.packages(file.path(path.to.main.src, "densityClust_0.3.3.tar.gz"), type = "source", repos = NULL)
install.packages(file.path(path.to.main.src, "XML_3.98-1.9.tar.gz"), type = "source", repos = NULL)

BiocManager::install("biocViews", version = "3.9")
BiocManager::install("monocle", version = "3.9")

