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
  if (!requireNamespace("glmnet", quietly=TRUE))   install.packages("glmnet",   repos="https://cloud.r-project.org")
  if (!requireNamespace("pROC", quietly=TRUE))     install.packages("pROC",     repos="https://cloud.r-project.org")
  if (!requireNamespace("MOFA2", quietly=TRUE))    BiocManager::install("MOFA2", ask=FALSE, update=FALSE)
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  library(optparse)

suppressPackageStartupMessages({
  if (!"AnnotationDbi" %in% loadedNamespaces())   library(AnnotationDbi)
  if (!"org.Hs.eg.db" %in% loadedNamespaces())    library(org.Hs.eg.db)
})

# Fucntion to convert ENSG to symbols
map_ensg_to_symbol <- function(ids) {
  ids0 <- sub("\\..*$","", ids)            # remove the .version
  is_ensg <- grepl("^ENSG", ids0)
  if (!any(is_ensg)) return(ids)           # if nothing to map

  keys <- unique(ids0[is_ensg])
  res <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = keys, keytype = "ENSEMBL", columns = "SYMBOL"),
    error = function(e) NULL
  )

  out <- ids
  if (!is.null(res) && nrow(res) > 0) {
    sym <- stats::setNames(res$SYMBOL, res$ENSEMBL)
    out[is_ensg] <- sym[ids0[is_ensg]]
  }
  out[is.na(out) | out == ""] <- ids[is.na(out) | out == ""]
  out
}

; library(readr); library(fs); library(glmnet); library(pROC); library(MOFA2)
})

opt_list <- list(
  make_option("--project", type="character", help="TCGA project, e.g. TCGA-BRCA"),
  make_option("--genes",   type="character", default=NA, help="Comma list of genes (e.g. TP53,KRAS). If NA, use config genes."),
  make_option("--views",   type="character", default="Expression,Mutation,CNV",
              help="Views to include (comma-separated). Default: Expression,Mutation,CNV"),
  make_option("--topn_expr", type="integer", default=0, help="Top-N reconstructed Expression features by variance (0=ALL)"),
  make_option("--topn_mut",  type="integer", default=0, help="Top-N reconstructed Mutation features by variance (0=ALL)"),
  make_option("--topn_cnv",  type="integer", default=0, help="Top-N reconstructed CNV features by variance (0=ALL)"),
  make_option("--nfolds",    type="integer", default=5, help="CV folds for glmnet (default 5)"),
  make_option("--seed",      type="integer", default=42, help="Random seed"),
  make_option("--results",   type="character", default="results", help="Base results dir (default results)"),
  make_option("--in_data",   type="character", default="data/raw", help="Base input dir for maf (for labels helper)"),
  make_option("--rscript",   type="character", default="Rscript", help="Path to Rscript (for labels helper)")
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.na(opt$project))

# Helper functions
short_tcga <- function(x) substr(as.character(x), 1, 15)

reconstruct_view <- function(F_mat, W_view) {
  fac <- intersect(colnames(F_mat), colnames(W_view))
  F2  <- as.matrix(F_mat[, fac, drop = FALSE])    # samples x K
  W2  <- as.matrix(W_view[, fac, drop = FALSE])   # features x K
  Xhat <- as.matrix(F2 %*% t(W2))                 # samples x features
  rownames(Xhat) <- rownames(F2)
  colnames(Xhat) <- rownames(W2)
  Xhat
}

select_topN_by_var <- function(X, n = 0) {
  if (is.null(n) || is.na(n) || n <= 0 || n >= ncol(X)) return(X)
  v <- apply(X, 2, var)
  X[, order(v, decreasing = TRUE)[seq_len(n)], drop = FALSE]
}

ensure_labels <- function(project, gene, results_dir="results", in_data="data/raw", rscript="Rscript") {
  outdir <- file.path(results_dir, project, "mofa", gene)
  fs::dir_create(outdir, recurse = TRUE)
  lpath <- file.path(outdir, "labels.tsv")
  if (file.exists(lpath)) return(lpath)
  f_factors <- file.path(results_dir, project, "mofa", "factors.tsv")
  cmd <- sprintf('%s scripts/utils/make_labels.R --project %s --gene %s --in_data %s --factors %s --outdir %s',
                 rscript, shQuote(project), shQuote(gene), shQuote(in_data), shQuote(f_factors), shQuote(outdir))
  message("Running: ", cmd)
  ok <- system(cmd)
  if (ok != 0 || !file.exists(lpath)) stop("Failed to create labels at: ", lpath)
  lpath
}

# Load MOFA models and extract factors and weights
mofa_dir <- file.path(opt$results, opt$project, "mofa")
model_path <- file.path(mofa_dir, "model_final.rds")
if (!file.exists(model_path)) stop("model_final.rds not found at: ", model_path)
model <- readRDS(model_path)

F_list <- MOFA2::get_factors(model, as.data.frame = FALSE)
F_mat  <- F_list[[1]]  # samples x factors
W_list <- MOFA2::get_weights(model, views = "all", as.data.frame = FALSE)  # list of matrices (features x factors)

views_req <- trimws(strsplit(opt$views, ",")[[1]])
views_req <- intersect(names(W_list), views_req)
if (!length(views_req)) stop("No valid views requested. Available: ", paste(names(W_list), collapse=", "))

# Rebuild views and store in DFs
X_view <- list()
if ("Expression" %in% views_req) {
  X_expr_hat <- reconstruct_view(F_mat, W_list[["Expression"]])
  X_view[["Expression"]] <- X_expr_hat
}
if ("CNV" %in% views_req) {
  X_cnv_hat <- reconstruct_view(F_mat, W_list[["CNV"]])
  X_view[["CNV"]] <- X_cnv_hat
}
if ("Mutation" %in% views_req) {
  X_mut_hat <- reconstruct_view(F_mat, W_list[["Mutation"]])
  X_view[["Mutation"]] <- X_mut_hat
}

# Selected to genes to be studied
genes <- if (!is.na(opt$genes)) trimws(strsplit(opt$genes, ",")[[1]]) else {
  # fallback to config if available
  cfg_path <- "config/config.yml"
  if (file.exists(cfg_path)) {
    y <- yaml::yaml.load_file(cfg_path)
    y$genes
  } else {
    stop("--genes not provided and config/config.yml missing")
  }
}

# loop over the genes
for (gene in genes) {
  geneU <- toupper(gene)
  message(">>> [LASSO] ", opt$project, " / ", geneU)

  # labels (lable samples as Gene_mut or Gene_wt)
  lpath <- ensure_labels(opt$project, geneU, results_dir=opt$results, in_data=opt$in_data, rscript=opt$rscript)
  lab <- readr::read_tsv(lpath, show_col_types = FALSE)
  y <- setNames(as.integer(lab$class == paste0(geneU, "_mut")), lab$sample)

  # Select top features per view 
  X_expr_hat <- X_view[["Expression"]]; X_cnv_hat <- X_view[["CNV"]]; X_mut_hat <- X_view[["Mutation"]]

  if (!is.null(X_expr_hat)) X_expr_hat <- select_topN_by_var(X_expr_hat, opt$topn_expr)
  if (!is.null(X_cnv_hat))  X_cnv_hat  <- select_topN_by_var(X_cnv_hat,  opt$topn_cnv)

  # leakage control: drop the outcome rows from Mutation view (e.g., TP53_MIS/TRUNC)
  if (!is.null(X_mut_hat)) {
    drop_rows <- paste0(geneU, c("_MIS","_TRUNC"))
    keep <- setdiff(colnames(X_mut_hat), drop_rows)
    X_mut_no_leak <- X_mut_hat[, keep, drop = FALSE]
    X_mut_no_leak <- select_topN_by_var(X_mut_no_leak, opt$topn_mut)
  } else {
    X_mut_no_leak <- NULL
  }

  # align samples across views 
  sample_sets <- list(rownames(F_mat), names(y))
  if (!is.null(X_expr_hat)) sample_sets <- c(sample_sets, list(rownames(X_expr_hat)))
  if (!is.null(X_cnv_hat))  sample_sets <- c(sample_sets, list(rownames(X_cnv_hat)))
  if (!is.null(X_mut_no_leak)) sample_sets <- c(sample_sets, list(rownames(X_mut_no_leak)))
  samples_common <- Reduce(intersect, sample_sets)
  if (length(samples_common) < 20) warning("Very few samples in common: ", length(samples_common))

  sel <- function(M) if (is.null(M)) NULL else M[samples_common, , drop = FALSE]

  X_expr_hat <- sel(X_expr_hat); X_cnv_hat <- sel(X_cnv_hat); X_mut_no_leak <- sel(X_mut_no_leak)
  y_vec <- y[samples_common]

  # tag predictors with view name (so far here: _CNV, _MIS. _TRUNC or _Expr)
  tag <- function(M, suf) { if (is.null(M)) return(NULL); colnames(M) <- paste0(colnames(M), "_", suf); M }
  X_expr_hat <- tag(X_expr_hat, "Expression")
  X_cnv_hat  <- tag(X_cnv_hat,  "CNV")
  X_mut_no_leak <- tag(X_mut_no_leak, "Mutation")

  # build X the top Weights DF
  mats <- Filter(Negate(is.null), list(X_expr_hat, X_cnv_hat, X_mut_no_leak))
  if (!length(mats)) stop("No predictors available; check views/topN")
  X <- do.call(cbind, mats)
  # drop NA and zero-variance
  X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
  sd_ok <- apply(X, 2, function(z) sd(z) > 0)
  X <- X[, sd_ok, drop = FALSE]
  if (ncol(X) == 0) stop("All predictors zero-variance after filtering.")

  # LASSO (binomial) with CV (AUC)
  set.seed(opt$seed)
  w_pos <- length(y_vec) / (2 * sum(y_vec))
  w_neg <- length(y_vec) / (2 * sum(1 - y_vec))
  cvfit <- glmnet::cv.glmnet(X, y_vec, family="binomial", alpha=1, nfolds=opt$nfolds,
                             type.measure="auc", standardize=TRUE) # weights=ifelse(y_vec==1, w_pos, w_neg)
  prob <- as.numeric(glmnet:::predict.cv.glmnet(cvfit, newx=X, s=cvfit$lambda.min, type="response"))
  AUC <- as.numeric(pROC::auc(pROC::roc(y_vec, prob)))

  # export
  out_dir <- file.path(opt$results, opt$project, "lasso", geneU)
  fs::dir_create(out_dir, recurse=TRUE)

  write_tsv(data.frame(metric="roc_auc", value=AUC), file.path(out_dir, "metrics.tsv"))
  # nonzero coefs
  nz <- coef(cvfit, s=cvfit$lambda.min)
  nz <- nz[as.numeric(nz)!=0, , drop=FALSE]
  if (length(nz) > 0) {
    base_feature <- sub("_(Expression|CNV|Mutation)$", "", rownames(nz))
    view <- sub("^.*_(Expression|CNV|Mutation)$", "\\1", rownames(nz))
    res_df <- data.frame(
      feature = rownames(nz),
      beta = as.numeric(nz),
      view = view,
      base_feature = base_feature,
      stringsAsFactors = FALSE
    )
    write_tsv(res_df, file.path(out_dir, "coefficients_nonzero.tsv"))
  }
  # Add symbols (needs Expression_feature_map.tsv produced in PrepareMOFA)
  map_path <- file.path(opt$results, opt$project, "mofa", "Expression_feature_map.tsv")
  expr_map <- if (file.exists(map_path)) readr::read_tsv(map_path, show_col_types = FALSE) else NULL

  res_df$symbol <- NA_character_

  # Convert ENSG to symbol
  idxE <- res_df$view == "Expression"
  if (!is.null(expr_map) && any(idxE)) {
      res_df$symbol[idxE] <- expr_map$symbol[match(res_df$base_feature[idxE], expr_map$feature)]
  }

  # Mutation features: TP53_MIS / TP53_TRUNC 
  idxM <- res_df$view == "Mutation"
  if (any(idxM)) {
      res_df$symbol[idxM] <- sub("_.*$", "", res_df$base_feature[idxM])
  }

# CNV features: if they’re ENSG, try the same map; otherwise leave NA 
  idxC <- res_df$view == "CNV" & is.na(res_df$symbol)
  if (!is.null(expr_map) && any(idxC)) {
      res_df$symbol[idxC] <- expr_map$symbol[match(res_df$base_feature[idxC], expr_map$feature)]
  }

  # save design actually used 
  write_tsv(data.frame(sample=samples_common, class=ifelse(y_vec==1, paste0(geneU,"_mut"), paste0(geneU,"_wt"))),
            file.path(out_dir, "labels_used.tsv"))
  write_tsv(as.data.frame(X), file.path(out_dir, "X_used.tsv"))
}

message("DONE.")


## normalize names: I only convert expression ENSG ids; I leave CNV and Mutation as they are
if (exists("X_expr_hat") && !is.null(X_expr_hat) &&
    any(grepl("^ENSG", colnames(X_expr_hat)))) {
  colnames(X_expr_hat) <- map_ensg_to_symbol(colnames(X_expr_hat))
}


# Select by top-N per view (0 = keep all)
if (exists("select_topN_by_var")) {
  if (exists("X_expr_hat") && !is.null(X_expr_hat) && opt$topn_expr > 0)
    X_expr_hat <- select_topN_by_var(X_expr_hat, opt$topn_expr)
  if (exists("X_mut_hat")  && !is.null(X_mut_hat)  && opt$topn_mut  > 0)
    X_mut_hat  <- select_topN_by_var(X_mut_hat,  opt$topn_mut)
  if (exists("X_cnv_hat")  && !is.null(X_cnv_hat)  && opt$topn_cnv  > 0)
    X_cnv_hat  <- select_topN_by_var(X_cnv_hat,  opt$topn_cnv)
}

