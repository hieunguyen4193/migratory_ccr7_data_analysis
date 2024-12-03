new.pkgs <- c(
  "Homo.sapiens",
  "TxDb.Hsapiens.UCSC.hg38.refGene",
  "Mus.musculus",
  "Orthology.eg.db"
)

for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    BiocManager::install(pkg, update = FALSE)
  }
}

library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.refGene
library(Mus.musculus)
library(Orthology.eg.db)
mapFun <- function(genes) {
  humeg <- mapIds(org.Hs.eg.db, genes, "ENTREZID","SYMBOL")
  mouseg <- mapIds(Orthology.eg.db, humeg, "Mus.musculus","Homo.sapiens")
  humdat <- select(Homo.sapiens, humeg, c("SYMBOL","CDSCHROM","CDSSTART","CDSEND", "CDSSTRAND"), "GENEID")
  mousedat <- select(Mus.musculus, mouseg,  c("SYMBOL","CDSCHROM","CDSSTART","CDSEND", "CDSSTRAND"), "GENEID")
  return(list(humdat = humdat, mousedat = mousedat))
}

mapFun("Apol7c")
