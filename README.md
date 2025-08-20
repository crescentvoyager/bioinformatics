## Volcano plot pipeline (DESeq2)

Files provided:
- `scripts/deseq2_volcano.R`: runs DESeq2 and saves a volcano plot and results.
- `samples.csv`: template listing your 20 samples. Fill `condition` with `control` or `case`.

Expected input (you provide):
- Gene-level counts matrix TSV, e.g. `salmon.merged.gene_counts.tsv`, with columns: `gene_id`, optional `gene_name`, then sample columns matching `samples.csv`.

Steps:
1) Place your counts file at a known path. By default the script expects:
   - `/Users/alexcrescent/Documents/Coding/bioinformatics/salmon.merged.gene_counts.tsv`
   Update the `counts_path` at the top of `scripts/deseq2_volcano.R` if your filename differs.
2) Open `samples.csv` and set each sample’s `condition` to `control` or `case`.
3) Run the script in R:
```
Rscript scripts/deseq2_volcano.R
```
4) Outputs will appear in `outputs/`:
   - `deseq2_results.csv`
   - `volcano.pdf` and `volcano.png`

Notes:
- The script installs required packages if missing (DESeq2, EnhancedVolcano).
- It filters lowly expressed genes (>=10 counts in at least 2 samples).
- Log2 fold-changes are shrunk with `apeglm`/`ashr` if available for better volcano visuals.
- Keep all samples in one counts matrix; use `samples.csv` to mark groups.


