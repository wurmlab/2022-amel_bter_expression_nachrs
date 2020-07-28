Expression of nAChR genes in *Apis mellifera* and *Bombus terrestris*
================
Federico Lopez
2020-07-28

  - [Load libraries](#load-libraries)
  - [Set a custom ggplot theme](#set-a-custom-ggplot-theme)
  - [Read `kallisto` results and table of samples and
    covariates.](#read-kallisto-results-and-table-of-samples-and-covariates.)
      - [*Apis mellifera* `kallisto`
        quantifications](#apis-mellifera-kallisto-quantifications)
      - [*Bombus terrestris* `kallisto`
        quantifications](#bombus-terrestris-kallisto-quantifications)
  - [Associate transcripts to genes](#associate-transcripts-to-genes)
  - [Combine input files and data
    tables](#combine-input-files-and-data-tables)
  - [Import `kallisto` quantifications using
    `tximport`](#import-kallisto-quantifications-using-tximport)
  - [Create `DESeqDataSet` objects and normalise
    counts](#create-deseqdataset-objects-and-normalise-counts)

# Load libraries

``` r
load_cran_pkgs <- function(pkg) {
  sapply(pkg, require, character.only = TRUE)
}
cran_pkgs <- c(
  "BiocManager", "tidyverse", "ggsignif", "gghighlight", "GGally", "ggfortify",
  "RColorBrewer", "pheatmap", "styler", "here"
)
load_cran_pkgs(cran_pkgs)
```

    ##  BiocManager    tidyverse     ggsignif  gghighlight       GGally    ggfortify 
    ##         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE 
    ## RColorBrewer     pheatmap       styler         here 
    ##         TRUE         TRUE         TRUE         TRUE

``` r
load_bioconductor_pkgs <- function(pkg) {
  sapply(pkg, require, character.only = TRUE)
}
bioconductor_pkgs <- c(
  "biomaRt", "rhdf5", "tximport", "DESeq2", "apeglm", "vsn", "hexbin"
)
load_bioconductor_pkgs(bioconductor_pkgs)
```

    ##  biomaRt    rhdf5 tximport   DESeq2   apeglm      vsn   hexbin 
    ##     TRUE     TRUE     TRUE     TRUE     TRUE     TRUE     TRUE

# Set a custom ggplot theme

``` r
# Set a custom ggplot theme
# devtools::install_github("cttobin/ggthemr")
library(ggthemr)
# Diverging
custom_palette <- c(
  "#9E9E9E", "#59BDEF", "#EBCC2A", "#E1AF00", "#F27C8D", "#7B73D1", "#7B9FE0",
  "#F5CC7F", "#66C2A5", "#28A7B6", "#F2CB57", "#F2A057", "#F2845C"
)
# ggthemr::colour_plot(custom_palette)

custom_theme <- ggthemr::define_palette(
  swatch = custom_palette,
  gradient = c(lower = "#59BDEF", upper = "#FF6B5A")
)
ggthemr::ggthemr(custom_theme)
```

``` r
# here::here()
# Create directory for results
dir.create(here::here("results", "amel_bter"), recursive = TRUE)
```

# Read `kallisto` results and table of samples and covariates.

## *Apis mellifera* `kallisto` quantifications

``` r
# Jasper et al. 2015
# Set sample identifiers using the names of the kallisto output directories
jasper_ids <- dir(here::here(
  "2015-jasper_worker_tissues", "2020-06-16-kallisto"
))

# Read metadata table
jasper_table <- read.csv(
  here::here(
    "2015-jasper_worker_tissues",
    "2015-jasper_worker_tissues_metadata.csv"
  ),
  header = TRUE,
  stringsAsFactors = TRUE
)

jasper_files <- here::here(
  "2015-jasper_worker_tissues", "2020-06-16-kallisto",
  jasper_table$sample,
  "abundance.h5"
)

names(jasper_files) <- jasper_table$sample

# Check that all files exist
all(file.exists(jasper_files))
```

    ## [1] TRUE

``` r
# Select worker "brain" samples
# jasper_table <- jasper_table %>%
#   dplyr::filter(tissue == "brain")
#
# jasper_files <- jasper_files[grep("brain", jasper_files)]

# Manfredini et al. 2015
# Set sample identifiers using the names of the kallisto output directories
manfredini_ids <- dir(here::here(
  "2015-manfredini_queen_brains", "2020-06-18-kallisto"
))

# Read metadata table
manfredini_table <- read.csv(
  here::here(
    "2015-manfredini_queen_brains",
    "2015-manfredini_queen_brains_metadata.csv"
  ),
  header = TRUE,
  stringsAsFactors = TRUE
)

manfredini_files <- here::here(
  "2015-manfredini_queen_brains", "2020-06-18-kallisto",
  manfredini_table$sample,
  "abundance.h5"
)

names(manfredini_files) <- manfredini_table$sample

# Check that all files exist
all(file.exists(manfredini_files))
```

    ## [1] TRUE

``` r
# Liberti et al. 2015
# Set sample identifiers using the names of the kallisto output directories
liberti_ids <- dir(here::here(
  "2015-liberti_queen_brains", "2020-06-18-kallisto"
))

# Read metadata table
liberti_table <- read.csv(
  here::here(
    "2015-liberti_queen_brains",
    "2015-liberti_queen_brains_metadata.csv"
  ),
  header = TRUE,
  stringsAsFactors = TRUE
)

liberti_files <- here::here(
  "2015-liberti_queen_brains", "2020-06-18-kallisto",
  liberti_table$sample,
  "abundance.h5"
)

names(liberti_files) <- liberti_table$sample

# Check that all files exist
all(file.exists(liberti_files))
```

    ## [1] TRUE

``` r
# Christen et al. 2015
# Set sample identifiers using the names of the kallisto output directories
christen_ids <- dir(here::here(
  "2018-christen_worker_brains", "2020-06-18-kallisto"
))

# Read metadata table
christen_table <- read.csv(
  here::here(
    "2018-christen_worker_brains",
    "2018-christen_worker_brains_metadata.csv"
  ),
  header = TRUE,
  stringsAsFactors = TRUE
)

christen_files <- here::here(
  "2018-christen_worker_brains", "2020-06-18-kallisto",
  christen_table$sample,
  "abundance.h5"
)

names(christen_files) <- christen_table$sample

# Check that all files exist
all(file.exists(christen_files))
```

    ## [1] TRUE

## *Bombus terrestris* `kallisto` quantifications

``` r
# Harrison et al. 2015
# Set sample identifiers using the names of the kallisto output directories
harrison_ids <- dir(here::here(
  "2015-harrison_queen_worker_whole_body",
  "2020-07-14-kallisto"
))

# Read metadata table
harrison_table <- read.csv(
  here::here(
    "2015-harrison_queen_worker_whole_body",
    "2015-harrison_queen_worker_whole_body_metadata.csv"
  ),
  header = TRUE,
  stringsAsFactors = TRUE
)

harrison_files <- here::here(
  "2015-harrison_queen_worker_whole_body",
  "2020-07-14-kallisto",
  harrison_table$sample,
  "abundance.h5"
)

names(harrison_files) <- harrison_table$sample

# Check that all files exist
all(file.exists(harrison_files))
```

    ## [1] TRUE

# Associate transcripts to genes

``` r
# Remove all biomaRt cached files
# biomaRt::biomartCacheInfo()
# biomaRt::biomartCacheClear()

# Retrieve Ensembl gene identifiers/names and create transcript-to-gene table
biomaRt::listMarts(host = "metazoa.ensembl.org")
```

    ##              biomart                       version
    ## 1       metazoa_mart      Ensembl Metazoa Genes 47
    ## 2 metazoa_variations Ensembl Metazoa Variations 47

``` r
# biomaRt::listMarts(host = "eg40-metazoa.ensembl.org") # version 40
# biomaRt::listEnsemblArchives()

amel_mart <- biomaRt::useMart(
  biomart = "metazoa_mart",
  dataset = "amellifera_eg_gene",
  host = "metazoa.ensembl.org" # version 47
)

amel_tx2gene <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
  mart = amel_mart
)

# Bombus terrestris
bter_mart <- biomaRt::useMart(
  biomart = "metazoa_mart",
  dataset = "bterrestris_eg_gene",
  host = "metazoa.ensembl.org" # version 47
)

bter_tx2gene <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
  mart = bter_mart
)

# biomaRt::listAttributes(metazoa_mart)
```

# Combine input files and data tables

``` r
# Apis mellifera
amel_files <- c(jasper_files, manfredini_files, liberti_files, christen_files)
amel_table <- dplyr::bind_rows(
  jasper_table, manfredini_table, liberti_table, christen_table
) %>%
  dplyr::select(sample, caste, tissue, study, species)

# Bombus terrestris
```

# Import `kallisto` quantifications using `tximport`

``` r
# Apis mellifera
amel_txi <- tximport::tximport(amel_files,
  type = "kallisto",
  tx2gene = amel_tx2gene,
  txOut = FALSE
)

# Bombus terrestris
bter_txi <- tximport::tximport(harrison_files,
  type = "kallisto",
  tx2gene = bter_tx2gene,
  txOut = FALSE
)
```

# Create `DESeqDataSet` objects and normalise counts

``` r
# Apis mellifera
amel_dds <- DESeq2::DESeqDataSetFromTximport(amel_txi,
  colData = amel_table,
  design = ~tissue
)

# Estimate the size factors
amel_dds <- DESeq2::estimateSizeFactors(amel_dds)
# Filtering
keep <- rowSums(DESeq2::counts(amel_dds) >= 10) >= 3
dds <- amel_dds[keep, ]
amel_counts <- DESeq2::counts(amel_dds, normalized = TRUE)

amel_nachrs <- c(
  "GB42850", "GB42644", "GB47845", "GB43275", "GB50159", "GB43416", "GB53053",
  "GB40923", "GB53427", "GB53055", "GB53428"
)

amel_nachrs_counts <- amel_counts %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ens_gene") %>%
  tidyr::gather(key = "sample", value = "norm_count", -ens_gene) %>%
  dplyr::mutate_if(is.numeric, round, digits = 4) %>%
  dplyr::filter(ens_gene %in% amel_nachrs) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(mean_count = mean(norm_count)) %>%
  dplyr::mutate(subunit = case_when(
    ens_gene == "GB42850" ~ "alpha1",
    ens_gene == "GB42644" ~ "alpha2",
    ens_gene == "GB47845" ~ "alpha3",
    ens_gene == "GB43275" ~ "alpha4",
    ens_gene == "GB50159" ~ "alpha5",
    ens_gene == "GB43416" ~ "alpha6",
    ens_gene == "GB53053" ~ "alpha7",
    ens_gene == "GB40923" ~ "alpha8",
    ens_gene == "GB53427" ~ "alpha9",
    ens_gene == "GB53055" ~ "beta1",
    ens_gene == "GB53428" ~ "beta2"
  )) %>%
  dplyr::mutate(study = case_when(
    stringr::str_detect(
      sample, paste(jasper_table$sample, collapse = "|")
    ) ~ "jasper",
    stringr::str_detect(
      sample, paste(manfredini_table$sample, collapse = "|")
    ) ~ "manfredini",
    stringr::str_detect(
      sample, paste(liberti_table$sample, collapse = "|")
    ) ~ "liberti",
    stringr::str_detect(
      sample, paste(christen_table$sample, collapse = "|")
    ) ~ "christen"
  )) %>%
  dplyr::mutate(tissue = case_when(
    grepl("brain", sample) ~ "brain",
    grepl("antenna", sample) ~ "antenna",
    grepl("midgut", sample) ~ "midgut",
    grepl("hypopharyngeal_gland", sample) ~ "hypopharyngeal_gland",
    grepl("malpighian_tubule", sample) ~ "malpighian_tubule",
    grepl("mandibular_gland", sample) ~ "mandibular_gland",
    grepl("muscle", sample) ~ "muscle",
    grepl("nasonov_gland", sample) ~ "nasonov_gland",
    grepl("thoracic_ganglion", sample) ~ "thoracic_ganglion",
    stringr::str_detect(
      sample, paste(manfredini_table$sample, collapse = "|")
    ) ~ "brain",
    stringr::str_detect(
      sample, paste(liberti_table$sample, collapse = "|")
    ) ~ "brain",
    stringr::str_detect(
      sample, paste(christen_table$sample, collapse = "|")
    ) ~ "brain"
  )) %>%
  dplyr::mutate(caste = case_when(
    stringr::str_detect(
      sample, paste(jasper_table$sample, collapse = "|")
    ) ~ "worker",
    stringr::str_detect(
      sample, paste(manfredini_table$sample, collapse = "|")
    ) ~ "queen",
    stringr::str_detect(
      sample, paste(liberti_table$sample, collapse = "|")
    ) ~ "queen",
    stringr::str_detect(
      sample, paste(christen_table$sample, collapse = "|")
    ) ~ "worker"
  ))

# For each study, calculate the mean of sum of counts per sample
amel_nachrs_mean <- amel_nachrs_counts %>%
  group_by(study, sample, tissue, caste) %>%
  summarise(sum_count = sum(norm_count)) %>%
  group_by(study, tissue, caste) %>%
  summarise(
    mean_sum = mean(sum_count),
    sd = sd(sum_count)
  )

head(amel_nachrs_mean)
```

    ## # A tibble: 6 x 5
    ## # Groups:   study, tissue [6]
    ##   study    tissue               caste  mean_sum     sd
    ##   <chr>    <chr>                <chr>     <dbl>  <dbl>
    ## 1 christen brain                worker    3859.  862. 
    ## 2 jasper   antenna              worker    1244.   86.8
    ## 3 jasper   brain                worker   14565. 2521. 
    ## 4 jasper   hypopharyngeal_gland worker    1284.  167. 
    ## 5 jasper   malpighian_tubule    worker    1335.  372. 
    ## 6 jasper   mandibular_gland     worker    2002.  627.

``` r
# Error bar plot using the mean of sum of normalised counts per sample
ggplot(amel_nachrs_mean, aes(
  x = tissue, y = mean_sum, color = tissue, fill = tissue
)) +
  geom_bar(
    stat = "identity", position = position_dodge(),
    width = 0.4, alpha = 0.9
  ) +
  geom_errorbar(aes(ymin = mean_sum - sd, ymax = mean_sum + sd),
    color = "#9E9E9E", width = 0.2, size = 0.6
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain"),
    strip.text = element_text(size = 10, face = "plain", hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    panel.spacing.x = unit(1, "mm")
  ) +
  facet_grid(. ~ study + caste, scales = "free", space = "free") +
  labs(
    x = "Tissue per caste per study",
    y = "Mean of sum of normalised nAChR counts"
  )
```

<img src="amel_bter_mean_expression_files/figure-gfm/unnamed-chunk-9-1.png" width=" extwidth" height="  extheight" />

``` r
ggsave(
  filename = here::here(
    "results", "amel_bter",
    "amel_bter_mean_norm_counts_bar_chart.pdf"
  ),
  width = 8,
  height = 5,
  units = "in"
)
```