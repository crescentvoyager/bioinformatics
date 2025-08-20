## Volcano plot pipeline (DESeq2)

### Current Analysis Parameters

**Data filtering:**
- Pre-analysis: Keep genes with ≥10 counts in ≥2 samples (standard DESeq2 filter)
- Plot display: Show genes with baseMean ≥20 (removes very low-expression noise)

**Statistical thresholds:**
- **p-value cutoff:** p < 0.05 (raw p-value, not FDR-adjusted)
- **Fold-change cutoff:** |log2FC| ≥ 0.58 (≥1.5× change up or down)
- **Why 1.5×:** More sensitive than industry standard 2×, captures moderate but meaningful changes

**Point classification:**
- **Red ("p-value and log2 FC"):** Both p < 0.05 AND ≥1.5× change (highest confidence)
- **Blue ("p-value"):** p < 0.05 but <1.5× change (significant but small effect)
- **Green ("log2 FC"):** ≥1.5× change but p ≥ 0.05 (large but not significant)
- **Grey ("NS"):** Neither threshold met

**Labeling:** Top 30 red genes by significance (p-value), then fold-change

### Files provided:
- `scripts/deseq2_volcano.R`: DESeq2 analysis with apeglm shrinkage
- `samples.csv`: 20 samples mapped to control/case groups
- `requirements.R`: Install all needed R packages

### Outputs:
- `outputs/volcano.pdf/.png`: Volcano plot visualization
- `outputs/deseq2_results.csv`: Full differential expression results
- `outputs/volcano_red_points.csv`: Red genes sorted by significance (follow-up candidates)
- `outputs/volcano_points_classification.csv`: All genes with categories and label flags

### Usage:
```bash
Rscript scripts/deseq2_volcano.R
```

### Setup:
1) Place your gene counts file: `salmon.merged.gene_counts.tsv`
2) Configure sample groups in `samples.csv` (control vs case)
3) Install packages: `source("requirements.R")` in R console
4) Run analysis: `Rscript scripts/deseq2_volcano.R`