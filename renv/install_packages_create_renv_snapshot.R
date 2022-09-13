install.packages("renv")

renv::init(
  project = "~/projects/2020_amel_bter_expression_nachrs",
  bare = TRUE, restart = TRUE
)

load_cran_pkgs <- function(pkg) {
  new_pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new_pkg)) {
    install.packages(new_pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

cran_pkgs <- c(
  "devtools", "XML", "BiocManager", "tidyverse", "GGally", "ggtext", "ggh4x",
  "ggforce", "egg", "patchwork", "RColorBrewer",  "styler", "here"
)

load_cran_pkgs(cran_pkgs)

load_bioconductor_pkgs <- function(pkg) {
  new_pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new_pkg)) {
    BiocManager::install(new_pkg, version = "3.10")
  }
  sapply(pkg, require, character.only = TRUE)
}

bioconductor_pkgs <- c(
  "biomaRt", "rhdf5", "tximport", "DESeq2", "apeglm", "RUVSeq", "limma", 
  "edgeR", "csaw", "vsn", "hexbin"
)

load_bioconductor_pkgs(bioconductor_pkgs)

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
devtools::install_github("cttobin/ggthemr")

# renv::settings$r.version("3.6.1")

renv::snapshot()
