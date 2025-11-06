#!/usr/bin/env Rscript
# Extract top weights (+predefined genes) from MOFA weights

# Define global libraries and the paths
lib_candidates <- c(
  "/scratch/yazbecka/mola-fsel-nf/R/x86_64-pc-linux-gnu-library/4.0/",                 # R library (mostly used on clusters)
  file.path(Sys.getenv("CONDA_PREFIX", unset=""), "lib/R/library"),                    # active conda env (in case using conda env)
  Sys.getenv("R_LIBS_USER"),                                                           # user library
  Sys.getenv("R_LIBS_SITE")                                                            # site library
)
lib_candidates <- lib_candidates[nzchar(lib_candidates) & dir.exists(lib_candidates)]
.libPaths(unique(c(lib_candidates, .libPaths())))

suppressPackageStartupMessages(library(reticulate))
use_condaenv("/scratch/yazbecka/mola-fsel-nf/CondaMofa/conda_envs/mofa", required = TRUE)

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr); library(purrr) })
source("scripts/00_globals.R")

opt_list <- list(
  make_option("--project", type="character"),
  make_option("--gene",    type="character"),
  make_option("--mofa_dir",type="character", default=NA),
  make_option("--outdir",  type="character", default=NA),
  make_option("--n_top",   type="integer",   default=20),
  make_option("--views",    type="character", default="Expression,Mutation,CNV"),
  make_option("--predef",      type="character", default=""),
  make_option("--predef_file", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.na(opt$project), !is.na(opt$gene))

mofa_dir <- ifelse(is.na(opt$mofa_dir), file.path("results", opt$project, "mofa", opt$gene), opt$mofa_dir)
out_dir  <- ifelse(is.na(opt$outdir),   file.path("results", opt$project, "weights", opt$gene), opt$outdir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


w_list <- readRDS(file.path(opt$mofa_dir, "weights_list.rds"))
# make sure matrices have names we can subset by 
fix_matrix <- function(W){
  W <- as.matrix(W)
  if (is.null(rownames(W))) rownames(W) <- paste0("feat_", seq_len(nrow(W)))
  if (is.null(colnames(W))) colnames(W) <- paste0("Factor", seq_len(ncol(W)))
  W
}

# predefined set from CLI gene(s); normalized per view (e.g., add _MIS/_TRUNC for Mutation)
predef_base <- unique(toupper(trimws(c(opt$gene))))
norm_predef_for_view <- function(view, rn){
  if (identical(view, "Mutation")) {
    intersect(unique(c(predef_base,
                       paste0(predef_base, "_MIS"),
                       paste0(predef_base, "_TRUNC"))), rn)
  } else {
    intersect(predef_base, rn)
  }
}

# n_top is how many top abs(weights) per factor to keep (plus predefined)
top_plus_predef <- function(W, factors, n_top, view, outdir){
  W  <- fix_matrix(W)
  rn <- rownames(W)
  factors <- intersect(factors, colnames(W))
  if (!length(factors)) return(invisible(NULL))

  dir.create(file.path(outdir, view), recursive = TRUE, showWarnings = FALSE)

  for (col in factors){
    w <- W[, col, drop = TRUE]
    names(w) <- rn

    ord  <- order(abs(w), decreasing = TRUE, na.last = NA)
    topN <- rn[ head(ord, n_top) ]

    pdef <- norm_predef_for_view(view, rn)
    keep <- intersect(unique(c(topN, pdef)), rn)

    df <- data.frame(
      view    = view,
      factor  = col,
      feature = keep,
      weight  = unname(w[keep]),
      stringsAsFactors = FALSE
    )

    readr::write_tsv(
      df,
      file.path(outdir, view, sprintf("top%d_plus_predef_%s_%s.tsv", n_top, view, col))
    )
  }
}

# run for each view
purrr::iwalk(w_list, function(W, view){
  top_plus_predef(
    W            = W,
    factors      = colnames(as.matrix(W)),
    n_top        = opt$n_top,
    view         = view,
    outdir       = out_dir  
  )
})
