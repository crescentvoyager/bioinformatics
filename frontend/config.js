// Simple config. Update defaults here or via UI controls.
window.APP_CONFIG = {
  // Absolute or relative path to a results CSV (summary or full export is fine)
  // Defaults to the latest independent filtering full export if present, otherwise a summary file.
  sources: [
    "outputs/deseq2_full_data_no_filter_fc1.5x_p0.05_pvalue.csv",
    "outputs/deseq2_full_data_independent_filter_fc1.5x_p0.05_pvalue.csv",
    "outputs/deseq2_results_independent_filter_fc1.5x_p0.05_pvalue.csv",
    "outputs/deseq2_results_manual_filter_10c2s_bm20_fc1.4x_p0.05_pvalue.csv"
  ],
  // Path to samples metadata (sample,condition)
  samplesCsv: "samples.csv",
  // Default patient/sample ID to highlight
  defaultPatientId: "AMATICA1219191846",
  // Column names to look for in result files
  columns: {
    geneId: "gene_id",
    geneName: "gene_name",
    log2fc: "log2FoldChange",
    pvalue: "pvalue",
    padj: "padj",
  }
};


