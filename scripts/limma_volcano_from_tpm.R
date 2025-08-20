#!/usr/bin/env Rscript

# Volcano plot from TPM using limma-trend (fallback when counts are not available)
# Inputs:
#   - Gene-level TPM CSV with columns: gene_id, gene_name, then sample columns
#   - Metadata CSV mapping sample -> condition (control, case)
# Outputs:
#   - outputs/limma_results_tpm.csv
#   - outputs/volcano_tpm.pdf and outputs/volcano_tpm.png

suppressPackageStartupMessages({
  ensure_pkg <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) BiocManager::install(pkg, ask = FALSE) else install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  ensure_pkg("limma")
  ensure_pkg("EnhancedVolcano", bioc = TRUE)
})

# ---------- Configuration ----------
tpm_path    <- "/Users/alexcrescent/Documents/Coding/bioinformatics/[alex]salmon.merged.gene_tpm - salmon.merged.gene_tpm.csv"
samples_path <- "/Users/alexcrescent/Documents/Coding/bioinformatics/samples.csv"
output_dir   <- "/Users/alexcrescent/Documents/Coding/bioinformatics/outputs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load input data ----------
stopifnot(file.exists(tpm_path))
stopifnot(file.exists(samples_path))

raw <- utils::read.csv(tpm_path, check.names = FALSE)
req_cols <- c("gene_id", "gene_name")
if (!all(req_cols %in% colnames(raw))) stop("TPM CSV must contain 'gene_id' and 'gene_name' columns.")

annot <- raw[, req_cols, drop = FALSE]
expr  <- raw[, setdiff(colnames(raw), req_cols), drop = FALSE]
rownames(expr) <- annot$gene_id

samples <- utils::read.csv(samples_path, stringsAsFactors = FALSE)
if (!all(c("sample", "condition") %in% colnames(samples))) stop("samples.csv must have columns: sample, condition")
if (!all(samples$sample %in% colnames(expr))) {
  missing <- setdiff(samples$sample, colnames(expr))
  stop(sprintf("These samples in samples.csv are not present in TPM table: %s", paste(missing, collapse = ", ")))
}
expr <- expr[, samples$sample, drop = FALSE]

samples$condition <- factor(samples$condition)
if (!all(c("control", "case") %in% levels(samples$condition))) stop("'condition' must include both 'control' and 'case'.")
samples$condition <- stats::relevel(samples$condition, ref = "control")

# ---------- limma-trend analysis ----------
logTPM <- log2(as.matrix(expr) + 0.5)
design <- stats::model.matrix(~ samples$condition)
fit <- limma::lmFit(logTPM, design)
fit <- limma::eBayes(fit, trend = TRUE)

tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
tt$adj.P.Val <- stats::p.adjust(tt$P.Value, method = "BH")

# Attach gene_name labels
gene_label_map <- stats::setNames(annot$gene_name, annot$gene_id)
tt$gene_name <- gene_label_map[rownames(tt)]

utils::write.csv(tt, file = file.path(output_dir, "limma_results_tpm.csv"), row.names = TRUE)

# ---------- Volcano plot ----------
lab_vec <- ifelse(is.na(tt$gene_name) | tt$gene_name == "", rownames(tt), tt$gene_name)

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  tt,
  lab = lab_vec,
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 3.0,
  title = "Case vs Control (TPM, limma-trend)",
  subtitle = "log2FC vs FDR"
)

grDevices::pdf(file.path(output_dir, "volcano_tpm.pdf"), width = 8, height = 6)
print(volcano_plot)
grDevices::dev.off()

grDevices::png(file.path(output_dir, "volcano_tpm.png"), width = 1600, height = 1200, res = 200)
print(volcano_plot)
grDevices::dev.off()

message(sprintf("Done. Results: %s | Plots: %s",
                file.path(output_dir, "limma_results_tpm.csv"),
                file.path(output_dir, "volcano_tpm.{pdf,png}")))


