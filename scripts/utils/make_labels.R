#!/usr/bin/env Rscript

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
  if (!requireNamespace("optparse", quietly=TRUE)) install.packages("optparse", repos="https://cloud.r-project.org")
  if (!requireNamespace("readr", quietly=TRUE))    install.packages("readr",    repos="https://cloud.r-project.org")
  if (!requireNamespace("fs", quietly=TRUE))       install.packages("fs",       repos="https://cloud.r-project.org")
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  if (!requireNamespace("maftools", quietly=TRUE)) BiocManager::install("maftools", ask=FALSE, update=FALSE)
  library(optparse); library(readr); library(fs); library(maftools)
})

short_tcga <- function(x) substr(as.character(x), 1, 15)

opt_list <- list(
  optparse::make_option("--project", type="character"),
  optparse::make_option("--gene",    type="character"),
  optparse::make_option("--in_data", type="character", default="data/raw"),
  optparse::make_option("--factors", type="character", help="Path to factors.tsv to align samples"),
  optparse::make_option("--outdir",  type="character", default=NA)
)
opt <- optparse::parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.na(opt$project), !is.na(opt$gene))

in_dir <- file.path(opt$in_data, opt$project)
maf_obj <- readRDS(file.path(in_dir, "mutations_maf.rds"))
factors <- readr::read_tsv(opt$factors, show_col_types = FALSE)

samples <- short_tcga(as.character(maf_obj@clinical.data$Tumor_Sample_Barcode))
samples <- intersect(samples, factors$sample)

d <- maf_obj@data
d$sample <- short_tcga(d$Tumor_Sample_Barcode)
d$gene   <- toupper(as.character(d$Hugo_Symbol))
d$cls    <- as.character(d$Variant_Classification)
trunc_cls <- c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site",
               "Nonstop_Mutation","Translation_Start_Site")
mis_cls <- c("Missense_Mutation")

target <- toupper(opt$gene)
keep <- d$gene == target & d$sample %in% samples & (d$cls %in% c(trunc_cls, mis_cls))
mutated <- unique(d$sample[keep])

df <- data.frame(
  sample = samples,
  class  = ifelse(samples %in% mutated, paste0(target, "_mut"), paste0(target, "_wt")),
  stringsAsFactors = FALSE
)

out_dir <- if (is.na(opt$outdir)) file.path("results", opt$project, "mofa", target) else opt$outdir
fs::dir_create(out_dir, recurse = TRUE)
readr::write_tsv(df, file.path(out_dir, "labels.tsv"))
message("Wrote labels: ", file.path(out_dir, "labels.tsv"), " (", sum(samples %in% mutated), " mut / ",
        length(samples) - sum(samples %in% mutated), " wt)")
