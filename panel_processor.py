#!/usr/bin/env python3
"""
Gene Panel Processor for Patient-Specific Expression Analysis
Handles filtering and organizing expression data for custom gene panels.
"""

import pandas as pd
import json
import os
from typing import List, Dict, Optional, Any
from flask import Flask, request, jsonify, render_template_string
from flask_cors import CORS
import logging

app = Flask(__name__)
CORS(app)
logging.basicConfig(level=logging.INFO)

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
            logging.info(f"Loading expression data from {self.data_file}")
            self.expression_data = pd.read_csv(self.data_file)
            
            logging.info(f"Loading sample metadata from {self.samples_file}")
            self.samples_meta = pd.read_csv(self.samples_file)
            
            # Clean up samples metadata
            self.samples_meta = self.samples_meta.dropna()
            
            logging.info(f"Loaded {len(self.expression_data)} genes and {len(self.samples_meta)} samples")
            
        except Exception as e:
            logging.error(f"Error loading data: {e}")
            raise
    
    def get_available_patients(self) -> List[str]:
        """Get list of available patient IDs."""
        return self.samples_meta['sample'].tolist()
    
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
            
            # Prepare results
            results = []
            for _, row in panel_data.iterrows():
                gene_result = {
                    'gene_id': row['gene_id'],
                    'gene_name': row.get('gene_name', 'Unknown'),
                    'patient_value': row.get(patient_id, 'N/A'),
                    'log2FoldChange': row.get('log2FoldChange', 'N/A'),
                    'pvalue': row.get('pvalue', 'N/A'),
                    'padj': row.get('padj', 'N/A'),
                    'baseMean': row.get('baseMean', 'N/A')
                }
                
                # Calculate summary statistics for other samples
                other_samples = [col for col in sample_cols if col != patient_id]
                other_values = []
                
                for sample in other_samples:
                    val = row.get(sample)
                    if pd.notna(val) and isinstance(val, (int, float)):
                        other_values.append(val)
                
                if other_values:
                    gene_result['mean_others'] = sum(other_values) / len(other_values)
                    gene_result['median_others'] = pd.Series(other_values).median()
                    gene_result['std_others'] = pd.Series(other_values).std()
                else:
                    gene_result['mean_others'] = 'N/A'
                    gene_result['median_others'] = 'N/A'
                    gene_result['std_others'] = 'N/A'
                
                results.append(gene_result)
            
            return {
                "success": True,
                "patient_id": patient_id,
                "patient_condition": patient_condition,
                "gene_count": len(results),
                "data": results
            }
            
        except Exception as e:
            logging.error(f"Error filtering panel data: {e}")
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

# Initialize processor with default data
DATA_FILE = "outputs/deseq2_full_data_no_filter_fc1.5x_p0.05_pvalue.csv"
SAMPLES_FILE = "samples.csv"

processor = None

@app.before_first_request
def initialize():
    global processor
    try:
        processor = PanelProcessor(DATA_FILE, SAMPLES_FILE)
        logging.info("Panel processor initialized successfully")
    except Exception as e:
        logging.error(f"Failed to initialize processor: {e}")

@app.route('/api/panels', methods=['GET'])
def get_panels():
    """Get list of available panels."""
    return jsonify(PREDEFINED_PANELS)

@app.route('/api/patients', methods=['GET'])
def get_patients():
    """Get list of available patients."""
    if not processor:
        return jsonify({"error": "Processor not initialized"}), 500
    
    try:
        patients = processor.get_available_patients()
        return jsonify({"patients": patients})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/api/panel/<panel_id>/patient/<patient_id>', methods=['GET'])
def get_panel_data(panel_id, patient_id):
    """Get panel data for specific patient."""
    if not processor:
        return jsonify({"error": "Processor not initialized"}), 500
    
    if panel_id not in PREDEFINED_PANELS:
        return jsonify({"error": "Panel not found"}), 404
    
    try:
        panel = PREDEFINED_PANELS[panel_id]
        result = processor.filter_panel_data(panel['genes'], patient_id)
        result['panel_name'] = panel['name']
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/panel/<panel_id>/patient/<patient_id>')
def panel_report(panel_id, patient_id):
    """Generate HTML report for panel data."""
    if not processor:
        return "Processor not initialized", 500
    
    if panel_id not in PREDEFINED_PANELS:
        return "Panel not found", 404
    
    try:
        panel = PREDEFINED_PANELS[panel_id]
        result = processor.filter_panel_data(panel['genes'], patient_id)
        
        if "error" in result:
            return f"Error: {result['error']}", 400
        
        # Generate HTML report
        html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ panel_name }} - Patient {{ patient_id }}</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
            color: #333;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .header {
            text-align: center;
            border-bottom: 3px solid #007acc;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }
        .header h1 {
            color: #007acc;
            margin: 0;
            font-size: 2.2em;
        }
        .header h2 {
            color: #666;
            margin: 10px 0 0 0;
            font-weight: normal;
        }
        .summary {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        .summary-item {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            border-left: 4px solid #007acc;
        }
        .summary-item h3 {
            margin: 0 0 5px 0;
            color: #007acc;
            font-size: 1.1em;
        }
        .summary-item p {
            margin: 0;
            font-size: 1.2em;
            font-weight: bold;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            background: white;
        }
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #007acc;
            color: white;
            font-weight: 600;
            text-transform: uppercase;
            font-size: 0.9em;
            letter-spacing: 0.5px;
        }
        tr:nth-child(even) {
            background-color: #f8f9fa;
        }
        tr:hover {
            background-color: #e8f4f8;
        }
        .numeric {
            text-align: right;
            font-family: 'Courier New', monospace;
        }
        .export-btn {
            position: fixed;
            top: 20px;
            right: 20px;
            background: #28a745;
            color: white;
            border: none;
            padding: 12px 24px;
            border-radius: 25px;
            cursor: pointer;
            font-size: 16px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.2);
            transition: all 0.3s ease;
        }
        .export-btn:hover {
            background: #218838;
            transform: translateY(-2px);
        }
        @media print {
            .export-btn { display: none; }
            body { background: white; }
            .container { box-shadow: none; }
        }
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/html2pdf.js/0.10.1/html2pdf.bundle.min.js"></script>
</head>
<body>
    <button class="export-btn" onclick="exportToPDF()">📄 Export PDF</button>
    
    <div class="container" id="report-content">
        <div class="header">
            <h1>{{ panel_name }}</h1>
            <h2>Patient: {{ patient_id }} ({{ patient_condition }})</h2>
        </div>
        
        <div class="summary">
            <div class="summary-item">
                <h3>Total Genes</h3>
                <p>{{ gene_count }}</p>
            </div>
            <div class="summary-item">
                <h3>Patient ID</h3>
                <p>{{ patient_id }}</p>
            </div>
            <div class="summary-item">
                <h3>Condition</h3>
                <p>{{ patient_condition }}</p>
            </div>
            <div class="summary-item">
                <h3>Report Date</h3>
                <p id="report-date"></p>
            </div>
        </div>
        
        <table>
            <thead>
                <tr>
                    <th>Gene ID</th>
                    <th>Gene Name</th>
                    <th>Patient Value</th>
                    <th>Log2 Fold Change</th>
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
                    <td><strong>{{ gene.gene_name }}</strong></td>
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
                margin: 1,
                filename: 'panel-report-{{ patient_id }}-{{ panel_id }}.pdf',
                image: { type: 'jpeg', quality: 0.98 },
                html2canvas: { scale: 2 },
                jsPDF: { unit: 'in', format: 'letter', orientation: 'portrait' }
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
        
    except Exception as e:
        logging.error(f"Error generating report: {e}")
        return f"Error generating report: {e}", 500

if __name__ == '__main__':
    app.run(debug=True, port=5001)
