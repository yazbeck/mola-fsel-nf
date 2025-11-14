#!/usr/bin/env Rscript

#  Define global libraries and the paths
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

# build MOFA inputs and run model then export factors and labels
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(stringr)
  library(SummarizedExperiment); library(MOFA2); library(maftools)
})
source("scripts/00_globals.R")

opt_list <- list(
  make_option("--project", type="character"),
  make_option("--gene",    type="character", default=NA),
  make_option("--in_data", type="character", default="data/real"),
  make_option("--out_mofa",type="character", default=NA),
  make_option("--labels",  type="character", default=NA, help="Optional path to labels TSV with 'sample' and 'class'"),
  make_option("--topn",  type="integer",  default = 50,
  help="Top mutated genes for mutation view (default 50; 0 = ALL genes)"),
  make_option("--predef", type="character", default = "",
  help="Comma-separated genes to ALWAYS include (e.g. 'PIK3CA,TP53,BRCA1,KEAP1')"),
# (optional) file input  a list file FYI: not tested by AY
  make_option("--predef_file", type="character", default = NA,
  help="Text file with one gene per line to always include")
)

opt <- parse_args(OptionParser(option_list=opt_list))
#stopifnot(!is.na(opt$project), !is.na(opt$gene))
stopifnot(!is.na(opt$project))
cfg <- get_cfg()

#in_dir  <- file.path(opt$in_data, opt$project, opt$gene)
in_dir  <- file.path(opt$in_data, opt$project)

#Check the existence of the files
req <- c("expression_se.rds","mutations_maf.rds","cnv_thresholded_by_genes.rds")
miss <- req[!file.exists(file.path(in_dir, req))]
if (length(miss)) {
  stop("Missing inputs in ", in_dir, ": ", paste(miss, collapse=", "),
       "\nDid you --all ?")
}

out_dir <- ifelse(is.na(opt$out_mofa), file.path("results", opt$project, "mofa", opt$gene), opt$out_mofa)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_se <- readRDS(file.path(in_dir, "expression_se.rds"))
maf_obj <- readRDS(file.path(in_dir, "mutations_maf.rds"))

# Build Expression matrix (genes x samples) using the counts  
#expr_mat <- as.matrix(SummarizedExperiment::assay(expr_se))
# Ensure samples in columns and have short IDs; Because samples are not the same length across the omics layers
short_tcga <- function(x) substr(as.character(x), 1, 15)
# Expression: log2(TPM + 1) 
expr_to_logtpm <- function(seObj) {
  # assumes TPM is stored in assay "tpm_unstrand"
  tpm <- SummarizedExperiment::assay(seObj, "tpm_unstrand")
  log2(tpm + 1)
}

# Convert ENSGIDs to GENE NAMES
set_rownames_to_genes <- function(expr_mat, seObj, column = "gene_name") {
  rd <- SummarizedExperiment::rowData(seObj)

  # named vector: names = ENSG IDs (rownames of SE), values = gene_name
  gene_map <- setNames(as.character(rd[[column]]), rownames(rd))

  # match by ENSG rownames
  new_names <- gene_map[rownames(expr_mat)]

  # fallback: if gene_name missing, keep original ENSG ID
  idx_na <- is.na(new_names) | new_names == ""
  new_names[idx_na] <- rownames(expr_mat)[idx_na]

  new_names <- make.unique(new_names)
  
  rownames(expr_mat) <- new_names
  expr_mat
}


expr_mat <- expr_to_logtpm(expr_se)
expr_mat <- set_rownames_to_genes(expr_mat, expr_se)
colnames(expr_mat) <- short_tcga(colnames(expr_mat))

# Extract the matrix of the mutations count of TRUN and MIS mutations
mut_mat_mis_trunc <- function(maf_obj, sample_ids, top_n = 50, predefined = character()) {
  d <- maf_obj@data
  d$sample <- substr(as.character(d$Tumor_Sample_Barcode), 1, 15)
  d$gene   <- as.character(d$Hugo_Symbol)
  d$cls    <- as.character(d$Variant_Classification)

  trunc_cls <- c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site",
                 "Nonstop_Mutation","Translation_Start_Site")
  mis_cls   <- c("Missense_Mutation")

  # choose genes: top_n by frequency + predefined (e.g., the target gene). Pass it as NULL then all genes are selected
  all_genes <- maftools::getGeneSummary(maf_obj)$Hugo_Symbol
  all_genes <- toupper(trimws(maftools::getGeneSummary(maf_obj)$Hugo_Symbol))
  predefined <- unique(toupper(trimws(predefined)))
  genes_top <- if (is.null(top_n) || is.na(top_n)) all_genes else head(all_genes, as.integer(top_n))
  genes <- unique(c(genes_top, predefined))

  
  
  # keep only requested samples/genes and classify MIS/TRUNC
  keep <- d$gene %in% genes & d$sample %in% sample_ids
  d <- d[keep, , drop = FALSE]
  d$kind <- ifelse(d$cls %in% mis_cls, "MIS",
                   ifelse(d$cls %in% trunc_cls, "TRUNC", NA))
  d <- d[!is.na(d$kind), , drop = FALSE]

  # counts per (gene, sample, kind)
  # build a matrix with rows "<GENE>_<KIND>" and columns "sample"
  d$row <- paste0(d$gene, "_", d$kind)
  M <- xtabs(~ row + sample, data = d)  # counts

  # All rows for requested genes & kinds, and all columns for sample_ids
  all_rows <- as.vector(outer(genes, c("MIS","TRUNC"), paste, sep = "_"))
  if (!length(M)) {
    M <- matrix(0, nrow = length(all_rows), ncol = length(sample_ids),
                dimnames = list(all_rows, sample_ids))
  } else {
    # add any missing rows/cols as zeros, and order deterministically
    miss_rows <- setdiff(all_rows, rownames(M))
    if (length(miss_rows)) M <- rbind(M, matrix(0, nrow = length(miss_rows), ncol = ncol(M),
                                                dimnames = list(miss_rows, colnames(M))))
    miss_cols <- setdiff(sample_ids, colnames(M))
    if (length(miss_cols)) M <- cbind(M, matrix(0, nrow = nrow(M), ncol = length(miss_cols),
                                                dimnames = list(rownames(M), miss_cols)))
    M <- M[all_rows, sample_ids, drop = FALSE]
  }
  return(M)
}

samples <- intersect(colnames(expr_mat), short_tcga(maftools::getSampleSummary(maf_obj)$Tumor_Sample_Barcode))
expr_mat <- expr_mat[, samples, drop=FALSE]
#predef  <- if (is.na(opt$gene)) character() else opt$gene

# Build predefined gene list: --predef_file + --predef + optional --gene
predef_vec <- character(0)

if (!is.na(opt$predef_file) && file.exists(opt$predef_file)) {
  predef_vec <- c(predef_vec, readr::read_lines(opt$predef_file))
}
if (!is.null(opt$predef) && nzchar(opt$predef)) {
  predef_vec <- c(predef_vec, unlist(strsplit(opt$predef, ",")))
}
if (!is.na(opt$gene)) {
  predef_vec <- c(predef_vec, opt$gene)
}

# normalize (trim + uppercase)
predef_vec <- unique(toupper(trimws(predef_vec)))


# topn handling: 0 → ALL genes (NULL), else integer
#mut_mat  <- mut_mat_mis_trunc(maf_obj, samples, top_n = 50, predefined = c(opt$gene))

tn <- if (is.null(opt$topn) || is.na(opt$topn) || opt$topn == 0) NULL else as.integer(opt$topn)
#predef_vec <- if (is.na(opt$gene)) character() else c(opt$gene)
mut_mat <- mut_mat_mis_trunc(maf_obj, samples, top_n = tn, predefined = predef_vec)

message("Mutation view: top_n = ", if (is.null(tn)) "ALL" else tn,
        " | genes = ", length(unique(sub("_.*$","", rownames(mut_mat)))),
        " | rows = ", nrow(mut_mat))

# CNV (Gistic) view
cnv_path <- file.path(in_dir, "cnv_thresholded_by_genes.rds")
has_cnv  <- file.exists(cnv_path)
if (has_cnv) {
  cnv_mat <- readRDS(cnv_path)
  if (!is.matrix(cnv_mat)) cnv_mat <- as.matrix(cnv_mat)
  colnames(cnv_mat) <- short_tcga(colnames(cnv_mat))
  # re-align all to CNV availability
  samples_common <- intersect(samples, colnames(cnv_mat))
  expr_mat <- expr_mat[, samples_common, drop=FALSE]
  mut_mat  <- mut_mat[,  samples_common, drop=FALSE]
  cnv_mat  <- cnv_mat[,  samples_common, drop=FALSE]
  samples  <- samples_common
}

# set seed
set.seed(42)
# Build views
views <- list(
  Expression = expr_mat,
  Mutation   = mut_mat,    # MIS/TRUNC counts
  CNV        = cnv_mat
)

likelihoods <- c(Expression = "gaussian", Mutation = "poisson", CNV = "gaussian") #poisson selected here for mutations as they are counts
saveRDS(views, file.path(out_dir, "mofa_inputs_views.rds"))  # Expression/Mutation/CNV matrices fed to MOFA

# Create MOFA object from matrices (features x samples) 
mofa <- create_mofa(views)

## For training a high K (25) is given, verbose, NO drop threshold (so it works with all TCGA cohort; hopefully)
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 25
model_opts$likelihoods <- likelihoods[names(views)]

train_opts <- get_default_training_options(mofa)
train_opts$seed    <- 42
train_opts$verbose <- TRUE
# train_opts$drop_factor_threshold <- 0.01  

message("# ---- Training initial model (25 factors) ----")
prep1      <- prepare_mofa(mofa, model_options = model_opts, training_options = train_opts)
model_raw  <- run_mofa(prep1, outfile = file.path(out_dir, "model_raw.hdf5"))
saveRDS(model_raw, file = file.path(out_dir, "model_raw.rds"))

# Save and plot VE (Variance Explained)
pdf(file.path(out_dir, "Model_raw_variance_explained.pdf"))
plot_variance_explained(model_raw, plot_total = TRUE)
dev.off()

## Choose k_best from VE (factor contributions)

ve_raw <- get_variance_explained(model_raw)

# r2_per_factor can be a list by group; take the single group
r2_pf <- ve_raw$r2_per_factor
if (is.list(r2_pf)) r2_pf <- r2_pf[[1]]          # e.g. $group1

# Ensure it's factors (rows) x views (cols)
if (nrow(r2_pf) > 0 && ncol(r2_pf) > 0) {
  if (!grepl("^Factor", rownames(r2_pf)[1]) && grepl("^Factor", colnames(r2_pf)[1])) {
    r2_pf <- t(r2_pf)
  }
}

# Thresholds are in PERCENT (%)
# keep if (≥1% in ≥2 views) OR (≥4% total across views)
total_per_factor <- rowSums(r2_pf, na.rm = TRUE)
keep <- which(rowSums(r2_pf > 1, na.rm = TRUE) >= 2 | total_per_factor > 4)

# Robust fallback in case nothing passes (rare)
if (length(keep) == 0) {
  keep <- order(total_per_factor, decreasing = TRUE)[1:min(2, nrow(r2_pf))]
}

keep_names <- rownames(r2_pf)[keep]
k_best <- length(keep)
message(sprintf("# ---- Keeping %d factors: %s",
                k_best, paste(keep_names, collapse = ", ")))

## Build with the best K selected above (from the training run)

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- k_best

# model_opts$likelihoods <- likelihoods[names(views)]

train_opts <- get_default_training_options(mofa)
train_opts$seed                  <- 42
train_opts$verbose               <- TRUE
train_opts$drop_factor_threshold <- 0.01   

prep2 <- prepare_mofa(
  mofa,
  model_options    = model_opts,
  training_options = train_opts
)
trained <- run_mofa(prep2, outfile = file.path(out_dir, "model_final.hdf5"))
saveRDS(trained, file.path(out_dir, "model_final.rds"))

pdf(file.path(out_dir, "Model_final_variance_explained.pdf"))
plot_variance_explained(trained, plot_total = TRUE)
dev.off()

ve_final <- get_variance_explained(trained)
saveRDS(ve_final, file.path(out_dir, "Model_final_variance_explained.rds"))

# Export factors (samples x factors)
factors <- as.data.frame(get_factors(trained, factors = "all", as.data.frame = FALSE)[[1]])
factors <- cbind(sample = rownames(factors), factors)
readr::write_tsv(factors, file.path(out_dir, "factors.tsv"))

#
if (!is.na(opt$labels) && file.exists(opt$labels)) {
  lab <- readr::read_tsv(opt$labels, show_col_types = FALSE)

} else if (!is.na(opt$gene)) {
  gene_rows <- grepl(paste0("^", opt$gene, "_"), rownames(mut_mat))
  status <- as.integer(colSums(mut_mat[gene_rows, samples, drop = FALSE]) > 0)
  lab <- data.frame(
    sample = samples,
    class  = ifelse(status == 1, paste0(opt$gene, "_mut"), paste0(opt$gene, "_wt"))
  )

} else {
  lab <- NULL
  message("Labels skipped (no --labels and no --gene).")
}

if (!is.null(lab)) {
  readr::write_tsv(lab, file.path(out_dir, "labels.tsv"))
}


# Save weights for later
w_list <- get_weights(trained, views = "all", factors = "all", as.data.frame = FALSE)
saveRDS(w_list, file.path(out_dir, "weights_list.rds"))
message("MOFA exported to: ", out_dir)
