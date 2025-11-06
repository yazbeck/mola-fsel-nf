#!/usr/bin/env Rscript
# Download & initial preparation (RNA-seq counts, MAF, CNV) and cache at "project level"

## Define global libraries and the paths
lib_candidates <- c(
  "/scratch/yazbecka/mola-fsel-nf/R/x86_64-pc-linux-gnu-library/4.0/",                 # R library (mostly used on clusters)
  file.path(Sys.getenv("CONDA_PREFIX", unset=""), "lib/R/library"),                    # active conda env (in case using conda env)
  Sys.getenv("R_LIBS_USER"),                                                           # user library
  Sys.getenv("R_LIBS_SITE")                                                            # site library
)
lib_candidates <- lib_candidates[nzchar(lib_candidates) & dir.exists(lib_candidates)]
.libPaths(unique(c(lib_candidates, .libPaths())))

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(stringr); library(fs)
})

option_list <- list(
  make_option("--project", type="character", help="TCGA project (e.g., TCGA-BRCA)"),
  make_option("--gene",    type="character", help="Gene symbol (e.g., TP53)"),
  make_option("--outdir",  type="character", default="data/real", help="Base output dir for per-gene folders [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
stopifnot(!is.na(opt$project), !is.na(opt$gene))


# Project-level cache dir 

cache_base <- file.path(dirname(opt$outdir), "raw")
cache_dir  <- file.path(cache_base, opt$project)
fs::dir_create(cache_dir, recurse=TRUE)

# Per-gene output dir 
out_dir <- file.path(opt$outdir, opt$project, opt$gene)
fs::dir_create(out_dir, recurse = TRUE)

short_tcga <- function(x) substr(as.character(x), 1, 15)

# Install libraries
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
biopkgs <- c("TCGAbiolinks","SummarizedExperiment","maftools")
for (p in biopkgs) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
suppressPackageStartupMessages({ library(TCGAbiolinks); library(SummarizedExperiment); library(maftools) })

# Helper
need <- function(fname) !file.exists(file.path(cache_dir, fname))

# 1) RNA-seq counts (will be cached)
if (need("expression_se.rds")) {
  message("[cache] fetching expression_se.rds for ", opt$project)
  q_expr <- GDCquery(project = opt$project,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     experimental.strategy = "RNA-Seq",
                     workflow.type = "STAR - Counts",
                     sample.type = "Primary Tumor")
  GDCdownload(q_expr)
  expr_se <- GDCprepare(q_expr)
  saveRDS(expr_se, file = file.path(cache_dir, "expression_se.rds"))
} else {
  message("[cache] hit expression_se.rds")
}

# 2) Masked MAF (will be cached)
if (need("mutations_df.rds") || need("mutations_maf.rds")) {
  message("[cache] fetching mutations for ", opt$project)
  q_mut <- GDCquery(project = opt$project,
                    data.category = "Simple Nucleotide Variation",
                    data.type = "Masked Somatic Mutation",
                    access = "open")
  GDCdownload(q_mut,files.per.chunk = 100)
  mut_df <- GDCprepare(q_mut)
  saveRDS(mut_df,  file = file.path(cache_dir, "mutations_df.rds"))
  maf_obj <- maftools::read.maf(mut_df)
  saveRDS(maf_obj, file = file.path(cache_dir, "mutations_maf.rds"))
} else {
  message("[cache] hit mutations_*")
}

# 3) CNV (GISTIC2 thresholds; will be cache)
proj_code <- sub("^TCGA-", "", opt$project)
cnv_gz <- file.path(cache_dir, paste0("TCGA.", proj_code, "_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"))
if (need("cnv_thresholded_by_genes.rds")) {
  url <- paste0("https://tcga.xenahubs.net/download/TCGA.", proj_code,
                ".sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz")
  message("[cache] downloading CNV from: ", url)
  utils::download.file(url, destfile = cnv_gz, mode = "wb", quiet = TRUE)
  if (file.exists(cnv_gz) && file.size(cnv_gz) > 0) {
    cnv_mat <- read.delim(gzfile(cnv_gz), stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
    colnames(cnv_mat) <- short_tcga(colnames(cnv_mat))
    saveRDS(cnv_mat, file = file.path(cache_dir, "cnv_thresholded_by_genes.rds"))
    readr::write_tsv(cbind(gene=rownames(cnv_mat), as.data.frame(cnv_mat, check.names=FALSE)),
                     file.path(cache_dir, "cnv_thresholded_by_genes.tsv"))
  } else {
    warning("[cache] CNV download failed: ", cnv_gz)
  }
} else {
  message("[cache] hit CNV")
}


# Copy/symlink cached files into per-gene folder (keeps downstream unchanged)

safe_link_or_copy <- function(src, dst){
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(src)) return(invisible(FALSE))
  if (file.exists(dst))  return(invisible(TRUE))
  ok <- tryCatch(file.link(src, dst), error=function(e) FALSE)  # symlink/hardlink when possible
  if (!isTRUE(ok)) ok <- file.copy(src, dst, overwrite = TRUE)  # fallback to copy
  invisible(ok)
}
for (f in c("expression_se.rds","mutations_df.rds","mutations_maf.rds","cnv_thresholded_by_genes.rds")){
  safe_link_or_copy(file.path(cache_dir, f), file.path(out_dir, f))
}

message("DONE DownloadPrepareMatrics -> cache: ", cache_dir, " | per-gene: ", out_dir)
