#!/usr/bin/env Rscript

# Differential expression (DESeq2) and volcano plot pipeline
# Inputs:
#   - Gene-level counts matrix TSV with columns: gene_id, optional gene_name, then sample columns
#   - Metadata CSV mapping sample -> condition (levels: control, case)
# Outputs:

#   - outputs/deseq2_results.csv
#   - outputs/volcano.pdf and outputs/volcano.png

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  ensure_pkg <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) BiocManager::install(pkg, ask = FALSE) else install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  ensure_pkg("DESeq2", bioc = TRUE)
  ensure_pkg("EnhancedVolcano", bioc = TRUE)
  # Optional, for ID→gene symbol mapping (human)
  ensure_pkg("AnnotationDbi", bioc = TRUE)
  ensure_pkg("org.Hs.eg.db", bioc = TRUE)
})

# ---------- Configuration ----------
counts_path <- "/Users/alexcrescent/Documents/Coding/bioinformatics/salmon.merged.gene_counts.tsv"  # <- update if needed
samples_path <- "/Users/alexcrescent/Documents/Coding/bioinformatics/samples.csv"
output_dir   <- "/Users/alexcrescent/Documents/Coding/bioinformatics/outputs"

# Plot tuning (adjust to taste)
use_shrunk_lfc_for_x <- FALSE   # Use unshrunk LFC on x for better spread while keeping shrunk stats
y_metric <- "pvalue"            # switch to raw p-values for better volcano shape
fc_cutoff <- 1.0                # 2x fold-change (industry standard)
xlim_range <- c(-5, 5)          # symmetric x-axis limits (wider view)
base_mean_min <- 20             # filter out very low-expression genes by baseMean

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load input data ----------
stopifnot(file.exists(counts_path))
stopifnot(file.exists(samples_path))

raw <- utils::read.delim(counts_path, check.names = FALSE)
if (!"gene_id" %in% colnames(raw)) stop("Input counts file must contain a 'gene_id' column.")
has_gene_name <- "gene_name" %in% colnames(raw)

annotation_cols <- c("gene_id", if (has_gene_name) "gene_name")
annot <- raw[, annotation_cols, drop = FALSE]
count_matrix <- raw[, setdiff(colnames(raw), annotation_cols), drop = FALSE]
rownames(count_matrix) <- raw$gene_id

# Coerce to numeric matrix and round (DESeq2 expects counts)
count_matrix <- as.matrix(data.frame(lapply(count_matrix, function(x) as.numeric(as.character(x))), check.names = FALSE))
# Restore gene_id rownames lost during coercion
rownames(count_matrix) <- raw$gene_id
if (anyNA(count_matrix)) stop("Counts matrix contains non-numeric values after coercion.")
count_matrix <- round(count_matrix)

coldata <- utils::read.csv(samples_path, stringsAsFactors = FALSE)
req_cols <- c("sample", "condition")
if (!all(req_cols %in% colnames(coldata))) stop("samples.csv must have columns: sample, condition")

# Check that all listed samples are present in counts
if (!all(coldata$sample %in% colnames(count_matrix))) {
  missing <- setdiff(coldata$sample, colnames(count_matrix))
  stop(sprintf("These samples in samples.csv are not present in counts: %s", paste(missing, collapse = ", ")))
}

# Align columns to coldata order
count_matrix <- count_matrix[, coldata$sample, drop = FALSE]

# Validate conditions
if (any(is.na(coldata$condition)) || any(coldata$condition == "")) {
  stop("Please fill 'condition' for all samples in samples.csv (use 'control' or 'case').")
}
coldata$condition <- factor(coldata$condition)
if (!all(c("control", "case") %in% levels(coldata$condition))) {
  stop("'condition' must include both 'control' and 'case'.")
}
coldata$condition <- relevel(coldata$condition, ref = "control")

# ---------- DESeq2 analysis ----------
dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)

# Filter low counts: keep genes with at least 10 counts in at least 2 samples
keep <- rowSums(DESeq2::counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# Preserve gene IDs after filtering; DESeq2 results will follow this order
gene_ids <- rownames(dds)

dds <- DESeq2::DESeq(dds)

# LFC shrinkage for better effect size estimates
coef_name <- "condition_case_vs_control"
res_shrunk <- NULL
if (requireNamespace("apeglm", quietly = TRUE)) {
  res_shrunk <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
} else if (requireNamespace("ashr", quietly = TRUE)) {
  res_shrunk <- DESeq2::lfcShrink(dds, coef = coef_name, type = "ashr")
} else {
  message("apeglm/ashr not installed; using normal shrinkage.")
  res_shrunk <- DESeq2::lfcShrink(dds, coef = coef_name, type = "normal")
}

# Also keep unshrunk results (for optional x-axis spread and baseMean)
res_unshrunk <- as.data.frame(DESeq2::results(dds, name = coef_name))

res_df <- as.data.frame(res_shrunk)

# Preserve gene_id explicitly using result rownames
res_df$gene_id <- rownames(res_shrunk)

# Optionally use unshrunk LFC on x for visualization while keeping FDR from shrunk fit
if (!use_shrunk_lfc_for_x) {
  res_df$log2FoldChange <- res_unshrunk$log2FoldChange[match(res_df$gene_id, rownames(res_unshrunk))]
}

# Attach baseMean for filtering if present
if ("baseMean" %in% colnames(res_unshrunk)) {
  res_df$baseMean <- res_unshrunk$baseMean[match(res_df$gene_id, rownames(res_unshrunk))]
}

# Attach gene labels from the original annotation by gene_id
if (has_gene_name) {
  res_df$gene_name <- annot$gene_name[match(res_df$gene_id, annot$gene_id)]
} else {
  res_df$gene_name <- NA_character_
}

# Order by adjusted p-value
res_df <- res_df[order(res_df$padj, res_df$pvalue), ]

# Filter low expression genes for plotting (keeps DESeq2 stats intact)
if ("baseMean" %in% colnames(res_df)) {
  res_df <- res_df[is.na(res_df$baseMean) | res_df$baseMean >= base_mean_min, ]
}

# Write results
utils::write.csv(res_df, file = file.path(output_dir, "deseq2_results.csv"), row.names = TRUE)

# ---------- Volcano plot ----------
lab_vec <- res_df$gene_name  # labels strictly use gene names

# Select a limited set of labels (most significant among those passing FC)
sig_idx <- which(!is.na(res_df[[y_metric]]) & res_df[[y_metric]] < 0.05 & abs(res_df$log2FoldChange) >= fc_cutoff)
order_sig <- sig_idx[order(res_df[[y_metric]][sig_idx], -abs(res_df$log2FoldChange)[sig_idx])]
select_labels <- lab_vec[order_sig]
select_labels <- unique(select_labels[select_labels != ""]) 
select_labels <- head(select_labels, 30)   # Show top 30 labels for readability

# Prepare classification and label tracking for export
fc_sig <- abs(res_df$log2FoldChange) >= fc_cutoff
p_sig  <- !is.na(res_df[[y_metric]]) & res_df[[y_metric]] < 0.05
category <- ifelse(fc_sig & p_sig, "p-value and log2 FC",
             ifelse(fc_sig & !p_sig, "log2 FC",
             ifelse(!fc_sig & p_sig, "p-value", "NS")))
labeled_flag <- lab_vec %in% select_labels

class_df <- data.frame(
  gene_id = rownames(res_df),
  gene_name = res_df$gene_name,
  log2FoldChange = res_df$log2FoldChange,
  pvalue = res_df$pvalue,
  padj = res_df$padj,
  baseMean = if ("baseMean" %in% colnames(res_df)) res_df$baseMean else NA_real_,
  category = category,
  labeled = labeled_flag,
  stringsAsFactors = FALSE
)
utils::write.csv(class_df, file = file.path(output_dir, "volcano_points_classification.csv"), row.names = FALSE)
red_df <- class_df[class_df$category == "p-value and log2 FC", ]
ord <- order(red_df[[y_metric]], -abs(red_df$log2FoldChange), na.last = TRUE)
red_df <- red_df[ord, ]
utils::write.csv(red_df, file = file.path(output_dir, "volcano_red_points.csv"), row.names = FALSE)

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  res_df,
  lab = lab_vec,
  x = "log2FoldChange",
  y = y_metric,
  pCutoff = 0.05,
  FCcutoff = fc_cutoff,
  pointSize = 2.0,
  labSize = 3.0,
  title = "Case vs Control (DESeq2)",
  subtitle = if (identical(y_metric, "padj")) "log2FC vs FDR" else "log2FC vs p-value",
  xlim = xlim_range,
  ylim = {
    neglog <- -log10(res_df[[y_metric]]);
    ymax <- suppressWarnings(quantile(neglog, probs = 0.995, na.rm = TRUE));
    if (is.na(ymax) || !is.finite(ymax)) ymax <- max(neglog, na.rm = TRUE);
    c(0, max(5, ymax + 0.5))
  },
  colAlpha = 0.7,
  borderWidth = 0.5,
  legendPosition = "right",
  col = c("#B0B0B0", "#7BC66F", "#4E79A7", "#E15759"),
  selectLab = select_labels,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = 100
)

# Save plots
grDevices::pdf(file.path(output_dir, "volcano.pdf"), width = 8, height = 6)
print(volcano_plot)
grDevices::dev.off()

grDevices::png(file.path(output_dir, "volcano.png"), width = 1600, height = 1200, res = 200)
print(volcano_plot)
grDevices::dev.off()

message(sprintf("Done. Results: %s | Plots: %s", 
                file.path(output_dir, "deseq2_results.csv"), 
                file.path(output_dir, "volcano.{pdf,png}")))


