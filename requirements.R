# R requirements for DESeq2 volcano plot pipeline
# Install with: source("requirements.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Core packages
BiocManager::install(c(
  "DESeq2",           # Differential expression analysis
  "EnhancedVolcano",  # Volcano plots
  "AnnotationDbi",    # Gene annotation interface
  "org.Hs.eg.db"     # Human gene annotations
), ask = FALSE)

# Optional shrinkage estimators (recommended)
BiocManager::install(c(
  "apeglm",           # Adaptive prior estimation for log fold changes
  "ashr"              # Adaptive shrinkage estimator
), ask = FALSE)

# Base R packages (should be included)
# - utils (read.csv, write.csv)
# - stats (setNames)
# - grDevices (pdf, png, dev.off)

cat("All packages installed successfully!\n")
cat("Session info:\n")
sessionInfo()
