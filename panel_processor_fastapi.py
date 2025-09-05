#!/usr/bin/env python3
"""
Gene Panel Processor for Patient-Specific Expression Analysis (FastAPI)
Handles filtering and organizing expression data for custom gene panels.
"""

import pandas as pd
import numpy as np
import os
from typing import List, Dict, Optional, Any
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from pydantic import BaseModel
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="Gene Panel Analysis API",
    description="API for processing gene panels and generating patient-specific reports",
    version="1.0.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Pydantic models for request/response validation
class PanelInfo(BaseModel):
    name: str
    genes: List[str]

class GeneData(BaseModel):
    gene_id: str
    gene_name: str
    patient_value: Any
    log2FoldChange: Any
    pvalue: Any
    padj: Any
    baseMean: Any
    control_range: str
    relative_position: str

class ProteinData(BaseModel):
    marker_name: str
    patient_value: Any
    control_range: str
    relative_position: str

class PanelDataResponse(BaseModel):
    success: bool
    patient_id: str
    patient_condition: str
    gene_count: int
    panel_name: str
    data: List[GeneData]

class ProteinDataResponse(BaseModel):
    success: bool
    patient_id: str
    patient_condition: str
    marker_count: int
    data: List[ProteinData]
    other_patients_data: List[ProteinData]

class PatientsResponse(BaseModel):
    patients: List[str]

class PanelProcessor:
    """Handles gene panel data processing and filtering."""
    
    def __init__(self, data_file: str, samples_file: str):
        self.data_file = data_file
        self.samples_file = samples_file
        self.expression_data = None
        self.samples_meta = None
        self.load_data()
    
    def load_data(self):
        """Load expression data and sample metadata."""
        try:
            logger.info(f"Loading expression data from {self.data_file}")
            self.expression_data = pd.read_csv(self.data_file)
            
            logger.info(f"Loading sample metadata from {self.samples_file}")
            self.samples_meta = pd.read_csv(self.samples_file)
            
            # Clean up samples metadata
            self.samples_meta = self.samples_meta.dropna()
            
            logger.info(f"Loaded {len(self.expression_data)} genes and {len(self.samples_meta)} samples")
            
        except Exception as e:
            logger.error(f"Error loading data: {e}")
            raise
    
    def get_available_patients(self) -> List[str]:
        """Get list of available patient IDs (excluding healthy controls)."""
        # Only return disease samples (condition != 'HC')
        disease_samples = self.samples_meta[self.samples_meta['condition'] != 'HC']
        return disease_samples['sample'].tolist()
    
    def get_patient_condition(self, patient_id: str) -> Optional[str]:
        """Get the condition for a specific patient."""
        patient_row = self.samples_meta[self.samples_meta['sample'] == patient_id]
        if not patient_row.empty:
            return patient_row['condition'].iloc[0]
        return None
    
    def filter_panel_data(self, gene_ids: List[str], patient_id: str) -> Dict[str, Any]:
        """
        Filter expression data for specific genes and patient.
        
        Args:
            gene_ids: List of gene IDs to include in panel
            patient_id: Target patient ID
            
        Returns:
            Dictionary containing filtered data and metadata
        """
        try:
            # Filter data to only include specified genes
            panel_data = self.expression_data[
                self.expression_data['gene_id'].isin(gene_ids)
            ].copy()
            
            if panel_data.empty:
                return {"error": "No genes found in dataset", "data": []}
            
            # Get patient condition
            patient_condition = self.get_patient_condition(patient_id)
            if not patient_condition:
                return {"error": f"Patient {patient_id} not found", "data": []}
            
            # Get sample columns (exclude metadata columns)
            metadata_cols = ['gene_id', 'gene_name', 'baseMean', 'log2FoldChange', 
                           'lfcSE', 'pvalue', 'padj']
            sample_cols = [col for col in panel_data.columns if col not in metadata_cols]
            
            # Get healthy control samples from metadata
            healthy_controls = self.samples_meta[
                self.samples_meta['condition'] == 'HC'
            ]['sample'].tolist()
            
            # Helper function to handle NaN values
            def clean_value(val):
                if pd.isna(val) or val is None:
                    return 'N/A'
                if isinstance(val, (int, float)):
                    if pd.isna(val) or not np.isfinite(val):
                        return 'N/A'
                    return float(val)
                return str(val)
            
            # Prepare results
            results = []
            for _, row in panel_data.iterrows():
                gene_result = {
                    'gene_id': str(row['gene_id']),
                    'gene_name': str(row.get('gene_name', 'Unknown')),
                    'patient_value': clean_value(row.get(patient_id, 'N/A')),
                    'log2FoldChange': clean_value(row.get('log2FoldChange', 'N/A')),
                    'pvalue': clean_value(row.get('pvalue', 'N/A')),
                    'padj': clean_value(row.get('padj', 'N/A')),
                    'baseMean': clean_value(row.get('baseMean', 'N/A'))
                }
                
                # Calculate control range (only healthy controls)
                control_samples = [col for col in sample_cols if col != patient_id and col in healthy_controls]
                control_values = []
                
                for sample in control_samples:
                    val = row.get(sample)
                    if pd.notna(val) and isinstance(val, (int, float)):
                        control_values.append(val)
                
                # Get patient value for relative position calculation
                patient_val = row.get(patient_id)
                patient_is_numeric = pd.notna(patient_val) and isinstance(patient_val, (int, float))
                
                if control_values:
                    min_control = min(control_values)
                    max_control = max(control_values)
                    gene_result['control_range'] = f"{clean_value(min_control)} - {clean_value(max_control)}"
                    
                    # Calculate relative position
                    if patient_is_numeric:
                        if patient_val < min_control:
                            gene_result['relative_position'] = 'Low'
                        elif patient_val > max_control:
                            gene_result['relative_position'] = 'High'
                        else:
                            gene_result['relative_position'] = 'Normal'
                    else:
                        gene_result['relative_position'] = 'N/A'
                else:
                    gene_result['control_range'] = 'N/A'
                    gene_result['relative_position'] = 'N/A'
                
                results.append(gene_result)
            
            return {
                "success": True,
                "patient_id": patient_id,
                "patient_condition": patient_condition,
                "gene_count": len(results),
                "data": results
            }
            
        except Exception as e:
            logger.error(f"Error filtering panel data: {e}")
            return {"error": str(e), "data": []}

class ProteinProcessor:
    """Handles protein marker data processing and filtering."""
    
    def __init__(self, data_file: str, samples_file: str):
        self.data_file = data_file
        self.samples_file = samples_file
        self.protein_data = None
        self.samples_meta = None
        self.load_data()
    
    def load_data(self):
        """Load protein data and sample metadata."""
        try:
            logger.info(f"Loading protein data from {self.data_file}")
            self.protein_data = pd.read_csv(self.data_file)
            
            logger.info(f"Loading protein sample metadata from {self.samples_file}")
            self.samples_meta = pd.read_csv(self.samples_file)
            
            # Clean up samples metadata and strip whitespace
            self.samples_meta.columns = self.samples_meta.columns.str.strip()
            self.protein_data.columns = self.protein_data.columns.str.strip()
            self.samples_meta = self.samples_meta.dropna()
            
            logger.info(f"Loaded {len(self.protein_data)} patients and {len(self.protein_data.columns)-1} protein markers")
            
        except Exception as e:
            logger.error(f"Error loading protein data: {e}")
            raise
    
    def get_available_patients(self) -> List[str]:
        """Get list of available patient IDs (excluding healthy controls)."""
        # Only return disease samples (condition != 'HC')
        disease_samples = self.samples_meta[self.samples_meta['condition'] != 'HC']
        return disease_samples['sample'].astype(str).tolist()
    
    def get_patient_condition(self, patient_id: str) -> Optional[str]:
        """Get the condition for a specific patient."""
        patient_row = self.samples_meta[self.samples_meta['sample'].astype(str) == str(patient_id)]
        if not patient_row.empty:
            return patient_row['condition'].iloc[0]
        return None
    
    def filter_protein_data(self, patient_id: str) -> Dict[str, Any]:
        """
        Filter protein data for specific patient.
        
        Args:
            patient_id: Target patient ID
            
        Returns:
            Dictionary containing filtered data and metadata
        """
        try:
            # Find patient row
            patient_row = self.protein_data[self.protein_data['Patient ID'].astype(str) == str(patient_id)]
            if patient_row.empty:
                return {"error": f"Patient {patient_id} not found", "data": []}
            
            patient_row = patient_row.iloc[0]
            
            # Get patient condition
            patient_condition = self.get_patient_condition(patient_id)
            if not patient_condition:
                return {"error": f"Patient {patient_id} not found in metadata", "data": []}
            
            # Get healthy control samples from metadata
            healthy_controls = self.samples_meta[
                self.samples_meta['condition'] == 'HC'
            ]['sample'].astype(str).tolist()
            
            # Get other patients (non-HC, excluding current patient)
            other_patients = self.samples_meta[
                (self.samples_meta['condition'] != 'HC') & 
                (self.samples_meta['sample'].astype(str) != str(patient_id))
            ]['sample'].astype(str).tolist()
            
            # Get protein marker columns (exclude Patient ID)
            protein_cols = [col for col in self.protein_data.columns if col != 'Patient ID']
            
            # Helper function to handle NaN values
            def clean_value(val):
                if pd.isna(val) or val is None:
                    return 'N/A'
                if isinstance(val, (int, float)):
                    if pd.isna(val) or not np.isfinite(val):
                        return 'N/A'
                    return float(val)
                return str(val)
            
            # Prepare results
            results = []
            other_patients_results = []
            
            for marker in protein_cols:
                marker_result = {
                    'marker_name': marker,
                    'patient_value': clean_value(patient_row[marker])
                }
                
                # Calculate control range (only healthy controls)
                control_values = []
                for control_id in healthy_controls:
                    control_row = self.protein_data[self.protein_data['Patient ID'].astype(str) == control_id]
                    if not control_row.empty:
                        val = control_row.iloc[0][marker]
                        if pd.notna(val) and isinstance(val, (int, float)):
                            control_values.append(val)
                
                # Calculate other patients values for this marker
                other_patient_values = []
                for other_id in other_patients:
                    other_row = self.protein_data[self.protein_data['Patient ID'].astype(str) == other_id]
                    if not other_row.empty:
                        val = other_row.iloc[0][marker]
                        if pd.notna(val) and isinstance(val, (int, float)):
                            other_patient_values.append(val)
                
                # Get patient value for relative position calculation
                patient_val = patient_row[marker]
                patient_is_numeric = pd.notna(patient_val) and isinstance(patient_val, (int, float))
                
                if control_values:
                    min_control = min(control_values)
                    max_control = max(control_values)
                    marker_result['control_range'] = f"{clean_value(min_control)} - {clean_value(max_control)}"
                    
                    # Calculate relative position for current patient
                    if patient_is_numeric:
                        if patient_val < min_control:
                            marker_result['relative_position'] = 'Low'
                        elif patient_val > max_control:
                            marker_result['relative_position'] = 'High'
                        else:
                            marker_result['relative_position'] = 'Normal'
                    else:
                        marker_result['relative_position'] = 'N/A'
                else:
                    marker_result['control_range'] = 'N/A'
                    marker_result['relative_position'] = 'N/A'
                
                # Calculate other patients average and relative position
                other_marker_result = {
                    'marker_name': marker,
                    'control_range': marker_result['control_range']
                }
                
                if other_patient_values:
                    other_avg = sum(other_patient_values) / len(other_patient_values)
                    other_marker_result['patient_value'] = clean_value(other_avg)
                    
                    # Calculate relative position for other patients average
                    if control_values:
                        if other_avg < min_control:
                            other_marker_result['relative_position'] = 'Low'
                        elif other_avg > max_control:
                            other_marker_result['relative_position'] = 'High'
                        else:
                            other_marker_result['relative_position'] = 'Normal'
                    else:
                        other_marker_result['relative_position'] = 'N/A'
                else:
                    other_marker_result['patient_value'] = 'N/A'
                    other_marker_result['relative_position'] = 'N/A'
                
                results.append(marker_result)
                other_patients_results.append(other_marker_result)
            
            return {
                "success": True,
                "patient_id": patient_id,
                "patient_condition": patient_condition,
                "marker_count": len(results),
                "data": results,
                "other_patients_data": other_patients_results
            }
            
        except Exception as e:
            logger.error(f"Error filtering protein data: {e}")
            return {"error": str(e), "data": []}

# Predefined panels
PREDEFINED_PANELS = {
    "infection_panel_1": {
        "name": "Infection Panel 1 test",
        "genes": [
            "ENSG00000165949", "ENSG00000137959", "ENSG00000185745", "ENSG00000134321",
            "ENSG00000119917", "ENSG00000115155", "ENSG00000119922", "ENSG00000133106",
            "ENSG00000149131", "ENSG00000089127", "ENSG00000126709", "ENSG00000229391",
            "ENSG00000130656", "ENSG00000055332", "ENSG00000204010", "ENSG00000179639",
            "ENSG00000156265", "ENSG00000250361", "ENSG00000206177", "ENSG00000198692",
            "ENSG00000223609", "ENSG00000100225", "ENSG00000156113", "ENSG00000153208",
            "ENSG00000105246", "ENSG00000182118", "ENSG00000100024", "ENSG00000174837",
            "ENSG00000204179", "ENSG00000183160", "ENSG00000124107", "ENSG00000163993",
            "ENSG00000124102"
        ]
    }
}

# Initialize processors with default data
DATA_FILE = "outputs/deseq2_full_data_no_filter_fc1.5x_p0.05_pvalue.csv"
SAMPLES_FILE = "samples.csv"
PROTEIN_DATA_FILE = "outputs/31-dataset.csv"
PROTEIN_SAMPLES_FILE = "31markersamples.csv"

# Global processor instances
processor = None
protein_processor = None

@app.on_event("startup")
async def startup_event():
    """Initialize the processors on startup."""
    global processor, protein_processor
    try:
        processor = PanelProcessor(DATA_FILE, SAMPLES_FILE)
        logger.info("Panel processor initialized successfully")
        
        protein_processor = ProteinProcessor(PROTEIN_DATA_FILE, PROTEIN_SAMPLES_FILE)
        logger.info("Protein processor initialized successfully")
    except Exception as e:
        logger.error(f"Failed to initialize processors: {e}")
        raise

@app.get("/", response_class=HTMLResponse)
async def root():
    """Root endpoint with API information."""
    return """
    <html>
        <head>
            <title>Gene Panel Analysis API</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; }
                .container { max-width: 800px; margin: 0 auto; }
                .endpoint { background: #f5f5f5; padding: 10px; margin: 10px 0; border-radius: 5px; }
                .method { font-weight: bold; color: #2563eb; }
            </style>
        </head>
        <body>
            <div class="container">
                <h1>🧬 Gene Panel Analysis API</h1>
                <p>FastAPI backend for gene panel processing and reporting.</p>
                
                <h2>Available Endpoints:</h2>
                
                <div class="endpoint">
                    <div class="method">GET /docs</div>
                    <div>Interactive API documentation (Swagger UI)</div>
                </div>
                
                <div class="endpoint">
                    <div class="method">GET /redoc</div>
                    <div>Alternative API documentation (ReDoc)</div>
                </div>
                
                <div class="endpoint">
                    <div class="method">GET /api/panels</div>
                    <div>Get list of available gene panels</div>
                </div>
                
                <div class="endpoint">
                    <div class="method">GET /api/patients</div>
                    <div>Get list of available patients</div>
                </div>
                
                <div class="endpoint">
                    <div class="method">GET /api/panel/{panel_id}/patient/{patient_id}</div>
                    <div>Get panel data for specific patient</div>
                </div>
                
                <div class="endpoint">
                    <div class="method">GET /panel/{panel_id}/patient/{patient_id}</div>
                    <div>Generate HTML report for panel data</div>
                </div>
                
                <p><strong>Frontend:</strong> <a href="http://localhost:8080">http://localhost:8080</a></p>
            </div>
        </body>
    </html>
    """

@app.get("/api/panels", response_model=Dict[str, PanelInfo])
async def get_panels():
    """Get list of available panels."""
    return PREDEFINED_PANELS

@app.get("/api/patients", response_model=PatientsResponse)
async def get_patients():
    """Get list of available patients."""
    if not processor:
        raise HTTPException(status_code=500, detail="Processor not initialized")
    
    try:
        patients = processor.get_available_patients()
        return PatientsResponse(patients=patients)
    except Exception as e:
        logger.error(f"Error getting patients: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/panel/{panel_id}/patient/{patient_id}")
async def get_panel_data(panel_id: str, patient_id: str):
    """Get panel data for specific patient."""
    if not processor:
        raise HTTPException(status_code=500, detail="Processor not initialized")
    
    if panel_id not in PREDEFINED_PANELS:
        raise HTTPException(status_code=404, detail="Panel not found")
    
    try:
        panel = PREDEFINED_PANELS[panel_id]
        result = processor.filter_panel_data(panel['genes'], patient_id)
        
        if "error" in result:
            raise HTTPException(status_code=400, detail=result["error"])
        
        result['panel_name'] = panel['name']
        return result
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting panel data: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/panel/{panel_id}/patient/{patient_id}", response_class=HTMLResponse)
async def panel_report(panel_id: str, patient_id: str):
    """Generate HTML report for panel data."""
    if not processor:
        raise HTTPException(status_code=500, detail="Processor not initialized")
    
    if panel_id not in PREDEFINED_PANELS:
        raise HTTPException(status_code=404, detail="Panel not found")
    
    try:
        panel = PREDEFINED_PANELS[panel_id]
        result = processor.filter_panel_data(panel['genes'], patient_id)
        
        if "error" in result:
            raise HTTPException(status_code=400, detail=result["error"])
        
        # Generate HTML report with webapp-style design
        html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ panel_name }} - Patient {{ patient_id }}</title>
    <style>
        :root{
          --bg:#0b0f14;--panel:#121821;--text:#e6edf3;--muted:#98a2b3;--accent:#7c9cff;--accent2:#f97316;--grid:#1f2937;
        }
        * { box-sizing: border-box; }
        html, body { height: 100%; }
        @import url('https://fonts.googleapis.com/css2?family=Open+Sans:wght@400;600&display=swap');
        body {
            margin: 0;
            background: var(--bg);
            color: var(--text);
            font-family: 'Open Sans', ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial;
            padding: 20px;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
        }
        .header {
            text-align: center;
            padding: 24px 20px;
            border-bottom: 1px solid #111826;
            margin-bottom: 20px;
        }
        .header h1 {
            margin: 0;
            font-size: 28px;
            letter-spacing: 0.3px;
            color: var(--accent);
        }
        .header h2 {
            margin: 6px 0 0;
            color: var(--muted);
            font-weight: normal;
        }
        .controls {
            display: flex;
            flex-wrap: wrap;
            gap: 14px;
            align-items: center;
            justify-content: center;
            padding: 12px 20px 8px;
            border-bottom: 1px solid #111826;
            background: var(--panel);
            margin-bottom: 20px;
            border-radius: 8px;
        }
        .control {
            display: flex;
            align-items: center;
            gap: 8px;
            color: var(--muted);
            font-size: 14px;
        }
        .control strong {
            color: var(--text);
        }
        .export-btn {
            position: fixed;
            top: 20px;
            right: 20px;
            background: linear-gradient(135deg, #28a745, #20c997);
            color: white;
            border: none;
            padding: 12px 24px;
            border-radius: 25px;
            cursor: pointer;
            font-size: 16px;
            font-weight: 600;
            box-shadow: 0 2px 10px rgba(0,0,0,0.3);
            transition: all 0.3s ease;
            z-index: 1000;
        }
        .export-btn:hover {
            background: linear-gradient(135deg, #218838, #1ea085);
            transform: translateY(-2px);
            box-shadow: 0 4px 15px rgba(0,0,0,0.4);
        }
        .panel-table {
            width: 100%;
            border-collapse: collapse;
            background: var(--panel);
            border-radius: 8px;
            overflow: hidden;
            border: 1px solid #1f2937;
        }
        .panel-table th {
            background: #1f2937;
            color: var(--text);
            padding: 12px;
            text-align: left;
            font-weight: 600;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .panel-table td {
            padding: 12px;
            border-bottom: 1px solid #111826;
            color: var(--text);
        }
        .panel-table tr:nth-child(even) {
            background: rgba(255, 255, 255, 0.02);
        }
        .panel-table tr:hover {
            background: rgba(124, 156, 255, 0.1);
        }
        .panel-table .numeric {
            text-align: right;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
        }
        .panel-table .gene-name {
            font-weight: 600;
            color: var(--accent);
        }
        .status {
            text-align: center;
            padding: 16px;
            color: var(--muted);
            margin-bottom: 16px;
        }
        @media print {
            :root {
                --bg: white;
                --panel: white;
                --text: #333;
                --muted: #666;
                --accent: #007acc;
            }
            .export-btn { display: none; }
            body { background: white; }
            .panel-table { border: 1px solid #ddd; }
            .panel-table th { background: #f8f9fa; color: #333; }
            .panel-table td { color: #333; border-color: #ddd; }
        }
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/html2pdf.js/0.10.1/html2pdf.bundle.min.js"></script>
</head>
<body>
    <button class="export-btn" onclick="exportToPDF()">📄 Export PDF</button>
    
    <div class="container" id="report-content">
        <div class="header">
            <h1>🧬 {{ panel_name }}</h1>
            <p class="subtitle">Gene expression analysis results</p>
        </div>
        
        <div class="controls">
            <div class="control">
                Panel: <strong>{{ panel_name }}</strong>
            </div>
            <div class="control">
                Patient: <strong>{{ patient_id }} ({{ patient_condition }})</strong>
            </div>
            <div class="control">
                Genes: <strong>{{ gene_count }}</strong>
            </div>
            <div class="control">
                Date: <strong id="report-date"></strong>
            </div>
        </div>
        
        <div class="status">Panel Results</div>
        
        <table class="panel-table">
            <thead>
                <tr>
                    <th>Gene ID</th>
                    <th>Gene Name</th>
                    <th>Patient Value</th>
                    <th>Log2 FC</th>
                    <th>P-value</th>
                    <th>Adj. P-value</th>
                    <th>Mean (Others)</th>
                    <th>Base Mean</th>
                </tr>
            </thead>
            <tbody>
                {% for gene in data %}
                <tr>
                    <td>{{ gene.gene_id }}</td>
                    <td class="gene-name">{{ gene.gene_name }}</td>
                    <td class="numeric">{{ format_number(gene.patient_value) }}</td>
                    <td class="numeric">{{ format_number(gene.log2FoldChange) }}</td>
                    <td class="numeric">{{ format_scientific(gene.pvalue) }}</td>
                    <td class="numeric">{{ format_scientific(gene.padj) }}</td>
                    <td class="numeric">{{ format_number(gene.mean_others) }}</td>
                    <td class="numeric">{{ format_number(gene.baseMean) }}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
    </div>
    
    <script>
        // Set current date
        document.getElementById('report-date').textContent = new Date().toLocaleDateString();
        
        function exportToPDF() {
            const element = document.getElementById('report-content');
            const opt = {
                margin: 0.75,
                filename: '{{ panel_name }}-{{ patient_id }}-report.pdf',
                image: { type: 'jpeg', quality: 0.98 },
                html2canvas: { 
                    scale: 2,
                    useCORS: true,
                    backgroundColor: '#0b0f14'
                },
                jsPDF: { 
                    unit: 'in', 
                    format: 'letter', 
                    orientation: 'landscape' 
                }
            };
            
            html2pdf().set(opt).from(element).save();
        }
    </script>
</body>
</html>
        """
        
        # Helper functions for formatting
        def format_number(val):
            if val == 'N/A' or pd.isna(val):
                return 'N/A'
            try:
                return f"{float(val):.3f}"
            except:
                return str(val)
        
        def format_scientific(val):
            if val == 'N/A' or pd.isna(val):
                return 'N/A'
            try:
                num = float(val)
                if num < 0.001:
                    return f"{num:.2e}"
                else:
                    return f"{num:.4f}"
            except:
                return str(val)
        
        # Render template
        from jinja2 import Template
        template = Template(html_template)
        
        return template.render(
            panel_name=panel['name'],
            patient_id=patient_id,
            patient_condition=result['patient_condition'],
            gene_count=result['gene_count'],
            data=result['data'],
            panel_id=panel_id,
            format_number=format_number,
            format_scientific=format_scientific
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error generating report: {e}")
        raise HTTPException(status_code=500, detail=f"Error generating report: {e}")

# Protein Analysis Endpoints
@app.get("/api/protein/patients", response_model=PatientsResponse)
async def get_protein_patients():
    """Get list of available protein analysis patients."""
    if not protein_processor:
        raise HTTPException(status_code=500, detail="Protein processor not initialized")
    
    try:
        patients = protein_processor.get_available_patients()
        return PatientsResponse(patients=patients)
    except Exception as e:
        logger.error(f"Error getting protein patients: {e}")
        raise HTTPException(status_code=500, detail=f"Error getting protein patients: {e}")

@app.get("/api/protein/patient/{patient_id}", response_model=ProteinDataResponse)
async def get_protein_data(patient_id: str):
    """Get protein analysis data for a specific patient."""
    if not protein_processor:
        raise HTTPException(status_code=500, detail="Protein processor not initialized")
    
    try:
        result = protein_processor.filter_protein_data(patient_id)
        
        if "error" in result:
            raise HTTPException(status_code=400, detail=result["error"])
        
        return ProteinDataResponse(**result)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting protein data: {e}")
        raise HTTPException(status_code=500, detail=f"Error getting protein data: {e}")

if __name__ == '__main__':
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=5001)
