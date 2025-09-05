# 🧬 Bioinformatics Analysis Pipeline

A comprehensive bioinformatics analysis pipeline featuring DESeq2 differential expression analysis, interactive gene visualization, and custom gene panel reporting with PDF export capabilities.

## 🚀 Quick Start

### 1. Start the Backend Server
```bash
python start_server.py
```

### 2. Start the Frontend (in a separate terminal)
```bash
python start_frontend.py
```

### 3. Access the Application
- Frontend: http://localhost:8080
- Backend API: http://localhost:5001
- API Documentation: http://localhost:5001/docs

## 📊 Features

### Gene Analysis Dashboard
- **Interactive Gene Plots**: Search and visualize individual gene expression across samples
- **Multiple Data Sources**: Switch between different DESeq2 analysis results
- **Patient-Specific Views**: Highlight specific patient values in context of all samples
- **Real-time Search**: Type-ahead gene search with suggestions

### Custom Gene Panels
- **Pre-defined Panels**: Ready-to-use gene panels for specific analysis (e.g., Infection Panel 1)
- **Patient-Specific Reports**: Generate detailed reports for any patient/panel combination
- **PDF Export**: Export professional reports with patient data and gene expression values
- **Interactive Tables**: Sortable tables with expression data, statistics, and comparisons

### Available Panels
- **Infection Panel 1 test**: 33 genes related to infection response analysis

## 🔬 Volcano Plot Pipeline (DESeq2)

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

### R Analysis Files:
- `scripts/deseq2_volcano_independent_filtering.R`: DESeq2 analysis with independent filtering
- `scripts/deseq2_volcano_manual_filtering.R`: DESeq2 analysis with manual filtering
- `scripts/deseq2_volcano_no_filtering.R`: DESeq2 analysis with no filtering
- `samples.csv`: 20 samples mapped to control/case groups
- `requirements.R`: Install all needed R packages

### Analysis Outputs:
- `outputs/volcano.pdf/.png`: Volcano plot visualization
- `outputs/deseq2_results.csv`: Full differential expression results
- `outputs/deseq2_full_data_no_filter_fc1.5x_p0.05_pvalue.csv`: Complete dataset for panel analysis
- `outputs/volcano_red_points.csv`: Red genes sorted by significance (follow-up candidates)
- `outputs/volcano_points_classification.csv`: All genes with categories and label flags

## 📁 Project Structure

```
bioinformatics/
├── 🐍 Backend (Python/FastAPI)
│   ├── panel_processor_fastapi.py  # Main FastAPI backend server
│   ├── panel_processor.py          # Legacy Flask version
│   ├── start_server.py             # Backend startup script
│   ├── requirements.txt            # Python dependencies
│   └── venv/                       # Virtual environment
├── 🌐 Frontend (HTML/CSS/JS)
│   ├── index.html             # Main dashboard interface
│   ├── app.js                 # Gene analysis functionality
│   ├── panels.js              # Panel management system
│   ├── config.js              # Configuration settings
│   ├── style.css              # Styling and themes
│   └── start_frontend.py      # Frontend server script
├── 📊 Data Files
│   ├── samples.csv            # Patient metadata
│   ├── salmon.merged.gene_counts.tsv  # Raw count data
│   └── outputs/               # Analysis results
└── 📝 R Scripts
    └── scripts/               # DESeq2 analysis scripts
```

## 🔧 Technical Details

### Backend (Python/FastAPI)
- **Framework**: FastAPI with automatic API documentation and type validation
- **Performance**: High-performance async framework with uvicorn server
- **Data Processing**: Pandas for data manipulation
- **Panel System**: Modular design for easy panel creation
- **PDF Generation**: HTML-to-PDF export functionality
- **API Documentation**: Interactive docs at `/docs` and `/redoc`
- **API Endpoints**:
  - `GET /api/panels` - List available panels
  - `GET /api/patients` - List available patients
  - `GET /api/panel/{panel_id}/patient/{patient_id}` - Get panel data
  - `GET /panel/{panel_id}/patient/{patient_id}` - HTML report view
  - `GET /docs` - Interactive API documentation (Swagger UI)
  - `GET /redoc` - Alternative API documentation (ReDoc)

### Frontend (Vanilla JavaScript)
- **No Framework Dependencies**: Pure HTML/CSS/JavaScript
- **Interactive Plotting**: Plotly.js for gene expression visualization
- **Real-time Search**: Dynamic gene filtering and suggestions
- **Responsive Design**: Mobile-friendly interface
- **Tab Navigation**: Separate views for gene analysis and panel management

### Data Processing
- **Tidy Functions**: Modular design for easy panel creation
- **Statistical Calculations**: Mean, median, standard deviation comparisons
- **Data Validation**: Comprehensive error handling and data quality checks
- **Performance Optimized**: Efficient data filtering and processing

## 🎯 Usage Examples

### Creating a New Gene Panel
1. Edit `panel_processor.py`
2. Add new panel to `PREDEFINED_PANELS` dictionary:
```python
"my_custom_panel": {
    "name": "My Custom Panel",
    "genes": ["ENSG00000123456", "ENSG00000789012", ...]
}
```
3. Restart the backend server

### Analyzing a Specific Patient
1. Open the frontend at http://localhost:8080
2. Navigate to "Gene Panels" tab
3. Select your panel and patient
4. Click "Generate Panel Report"
5. View results and export PDF as needed

### Running R Analysis
```bash
# Run different filtering approaches
Rscript scripts/deseq2_volcano_independent_filtering.R
Rscript scripts/deseq2_volcano_manual_filtering.R
Rscript scripts/deseq2_volcano_no_filtering.R
```

## 📋 Requirements

### Python Dependencies
- FastAPI >=0.104.0
- Uvicorn[standard] >=0.24.0  
- Pandas >=2.0.0
- Jinja2 >=3.1.0
- NumPy >=1.24.0
- Python-multipart >=0.0.6

### R Dependencies
- DESeq2
- ggplot2
- dplyr
- readr

## 🛠️ Development Setup

### Initial Setup:
1. Ensure you have Python 3.7+ installed
2. Place your gene counts file: `salmon.merged.gene_counts.tsv`
3. Configure sample groups in `samples.csv` (control vs case)
4. Run `python start_server.py` to install dependencies and start backend
5. Run `python start_frontend.py` to start the frontend

### Adding New Features:
- **New Panels**: Edit `PREDEFINED_PANELS` in `panel_processor.py`
- **Frontend Features**: Add to `panels.js` or `app.js`
- **Styling**: Update `style.css`
- **API Endpoints**: Add routes to `panel_processor.py`

## 📄 License

This project is for research and educational purposes.