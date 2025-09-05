/**
 * Protein Analysis Module
 * Handles 31 Marker Protein Analysis functionality
 */
(function() {
    'use strict';
    
    const API_BASE = 'http://localhost:5001/api';
    
    // DOM elements
    const elements = {
        proteinPatientInput: document.getElementById('protein-patient-input'),
        proteinExportPdfBtn: document.getElementById('protein-export-pdf-btn'),
        proteinStatus: document.getElementById('protein-status'),
        proteinPreview: document.getElementById('protein-preview'),
        proteinDataBody: document.getElementById('protein-data-body'),
        proteinDescription: document.getElementById('protein-description')
    };
    
    // State management
    let currentPatientId = null;
    let currentProteinData = null;
    let radarChart = null;

    // Protein marker categories based on Amatica Health functional groupings
    const PROTEIN_CATEGORIES = {
        'Vascular & Viral': {
            markers: ['ROCK1 (pg/ml)', 'ROCK2 (pg/ml)', 'ARG1 (pg/ml)'],
            color: '#E74C3C'
        },
        'RAAS System': {
            markers: ['Ang(1-7) (pg/ml)', 'ACE (pg/ml)', 'ACE2 (pg/ml)', 'ANG I (pg/ml)', 'ANG II (pg/ml)'],
            color: '#3498DB'
        },
        'TGF-β System': {
            markers: ['TGFB1 (pg/ml)', 'f.TGFB1 (pg/ml)', 'TGFB2 (pg/ml)', 'TGFB3 (pg/ml)', 'GDF-15 (pg/ml)', 'Activin B (pg/ml)', 'Follistatin (pg/ml)'],
            color: '#9B59B6'
        },
        'Mitochondrial': {
            markers: ['DRP1 (pg/ml)', 'PINK1 (pg/ml)', 'SIRT1 pg/ml', 'HIF1a (pg/ml)'],
            color: '#F39C12'
        },
        'Immune System': {
            markers: ['PGE2 (pg/ml)', 'IFNa (pg/ml)', 'IFNB (pg/ml)', 'IFNG (pg/ml)', 'IFNL1 (pg/ml)', 'aNAGA (pg/ml)'],
            color: '#27AE60'
        },
        'Neurological': {
            markers: ['NEFL (pg/ml)', 'S100B (pg/ml)', 'Serotonin (pg/ml)', 'ATG13 (pg/ml)', 'BH2 (pg/ml)', 'BH4 (pg/ml)'],
            color: '#E67E22'
        },
        'Exercise & Muscle': {
            markers: ['TWEAK (pg/ml)'],
            color: '#1ABC9C'
        }
    };
    
    /**
     * Initialize protein analysis module
     */
    function init() {
        console.log('Initializing protein analysis module...');
        
        // Load available patients
        loadProteinPatients();
        
        // Setup event listeners
        setupEventListeners();
        
        console.log('Protein analysis module initialized successfully');
    }
    
    /**
     * Setup event listeners
     */
    function setupEventListeners() {
        // Patient selection change
        elements.proteinPatientInput.addEventListener('change', (e) => {
            currentPatientId = e.target.value;
            loadProteinData();
        });
        
        // Export PDF button
        elements.proteinExportPdfBtn.addEventListener('click', () => {
            if (currentPatientId) {
                // Open the report page in a new tab
                const reportUrl = `protein-report.html?patient=${currentPatientId}&autoprint=true`;
                window.open(reportUrl, '_blank');
            }
        });
    }
    
    /**
     * Load available protein patients
     */
    async function loadProteinPatients() {
        try {
            setProteinStatus('Loading patients...');
            
            const response = await fetch(`${API_BASE}/protein/patients`);
            if (!response.ok) throw new Error('Failed to fetch patients');
            
            const data = await response.json();
            console.log('Protein patients loaded:', data);
            
            // Populate patient dropdown
            elements.proteinPatientInput.innerHTML = '<option value="">Select patient...</option>';
            data.patients.forEach(patient => {
                const option = document.createElement('option');
                option.value = patient;
                option.textContent = `Patient ${patient}`;
                elements.proteinPatientInput.appendChild(option);
            });
            
            setProteinStatus('');
            
        } catch (error) {
            console.error('Error loading protein patients:', error);
            setProteinStatus(`Error loading patients: ${error.message}`);
        }
    }
    
    /**
     * Load protein data for selected patient
     */
    async function loadProteinData() {
        if (!currentPatientId) {
            elements.proteinPreview.style.display = 'none';
            elements.proteinExportPdfBtn.style.display = 'none';
            elements.proteinDescription.style.display = 'none';
            return;
        }
        
        try {
            setProteinStatus('Loading protein analysis data...');
            
            const response = await fetch(`${API_BASE}/protein/patient/${currentPatientId}`);
            if (!response.ok) throw new Error('Failed to fetch protein data');
            
            const data = await response.json();
            console.log('Protein data loaded:', data);
            
            // Store data globally
            currentProteinData = data;
            
            // Display the data
            displayProteinData(data);
            
            // Create radar chart
            createRadarChart(data);
            
            // Show description, preview and export button
            elements.proteinDescription.style.display = 'block';
            elements.proteinPreview.style.display = 'block';
            elements.proteinExportPdfBtn.style.display = 'inline-block';
            
            setProteinStatus(`Loaded ${data.marker_count} protein markers for Patient ${currentPatientId}`);
            
        } catch (error) {
            console.error('Error loading protein data:', error);
            setProteinStatus(`Error: ${error.message}`);
        }
    }
    
    /**
     * Display protein data in table
     */
    function displayProteinData(data) {
        // Clear existing table data
        elements.proteinDataBody.innerHTML = '';
        
        // Populate table with protein data
        data.data.forEach(protein => {
            const row = document.createElement('tr');
            row.innerHTML = `
                <td class="protein-name">${protein.marker_name}</td>
                <td class="numeric">${formatNumber(protein.patient_value)}</td>
                <td class="numeric">${protein.control_range || 'N/A'}</td>
                <td class="relative-position ${getRelativePositionClass(protein.relative_position)}">${protein.relative_position || 'N/A'}</td>
            `;
            elements.proteinDataBody.appendChild(row);
        });
    }
    
    /**
     * Set status message for protein operations
     */
    function setProteinStatus(message) {
        if (elements.proteinStatus) {
            elements.proteinStatus.textContent = message;
        }
    }
    
    /**
     * Format numeric values for display
     */
    function formatNumber(value) {
        if (value === 'N/A' || value === null || value === undefined) {
            return 'N/A';
        }
        
        try {
            const num = parseFloat(value);
            if (isNaN(num)) return 'N/A';
            
            // Format large numbers with commas
            if (num >= 1000) {
                return num.toLocaleString(undefined, {
                    minimumFractionDigits: 0,
                    maximumFractionDigits: 2
                });
            }
            return num.toFixed(2);
        } catch {
            return 'N/A';
        }
    }
    
    /**
     * Get CSS class for relative position styling
     */
    function getRelativePositionClass(position) {
        switch(position) {
            case 'Low': return 'position-low';
            case 'High': return 'position-high';
            case 'Normal': return 'position-normal';
            default: return 'position-na';
        }
    }

    /**
     * Create radar chart for functional categories
     */
    function createRadarChart(data) {
        const ctx = document.getElementById('protein-radar-chart');
        if (!ctx) return;

        // Destroy existing chart
        if (radarChart) {
            radarChart.destroy();
        }

        // Calculate category aggregations
        const categoryData = calculateCategoryAggregations(data.data, data.other_patients_data);
        
        const chartData = {
            labels: Object.keys(PROTEIN_CATEGORIES),
            datasets: [
                {
                    label: 'You',
                    data: categoryData.patient,
                    backgroundColor: 'rgba(142, 68, 173, 0.2)',
                    borderColor: '#8E44AD',
                    borderWidth: 3,
                    pointBackgroundColor: '#8E44AD',
                    pointBorderColor: '#8E44AD',
                    pointRadius: 6,
                    pointHoverRadius: 8
                },
                {
                    label: 'Healthy Controls',
                    data: categoryData.controls,
                    backgroundColor: 'rgba(68, 114, 196, 0.2)',
                    borderColor: '#4472C4',
                    borderWidth: 3,
                    pointBackgroundColor: '#4472C4',
                    pointBorderColor: '#4472C4',
                    pointRadius: 6,
                    pointHoverRadius: 8
                },
                {
                    label: 'Other Patients',
                    data: categoryData.otherPatients,
                    backgroundColor: 'rgba(243, 156, 18, 0.2)',
                    borderColor: '#F39C12',
                    borderWidth: 3,
                    pointBackgroundColor: '#F39C12',
                    pointBorderColor: '#F39C12',
                    pointRadius: 6,
                    pointHoverRadius: 8
                }
            ]
        };

        const config = {
            type: 'radar',
            data: chartData,
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: {
                        display: false // We have custom legend
                    },
                    tooltip: {
                        backgroundColor: 'rgba(17, 17, 17, 0.9)',
                        titleColor: '#fff',
                        bodyColor: '#fff',
                        borderColor: '#333',
                        borderWidth: 1,
                        callbacks: {
                            label: function(context) {
                                const label = context.dataset.label;
                                const value = context.parsed.r;
                                if (label === 'You') {
                                    return `You: ${value.toFixed(1)}% vs control average`;
                                } else if (label === 'Healthy Controls') {
                                    return `Control Average: ${value.toFixed(1)}% (baseline)`;
                                } else if (label === 'Other Patients') {
                                    return `Other Patients: ${value.toFixed(1)}% vs control average`;
                                } else {
                                    return `${label}: ${value.toFixed(1)}% vs control average`;
                                }
                            }
                        }
                    }
                },
                scales: {
                                            r: {
                            beginAtZero: true,
                            max: 200, // 200% of normal range
                            min: 0,
                            ticks: {
                                stepSize: 50,
                                color: '#98a2b3',
                                font: {
                                    size: 12
                                },
                                display: false // Hide the radial labels (50%, 100%, 150%)
                            },
                        grid: {
                            color: '#1f2937',
                            lineWidth: 1
                        },
                        angleLines: {
                            color: '#1f2937',
                            lineWidth: 1
                        },
                        pointLabels: {
                            color: '#e6edf3',
                            font: {
                                size: 13,
                                weight: '600'
                            }
                        }
                    }
                },
                elements: {
                    line: {
                        tension: 0.1
                    }
                }
            }
        };

        radarChart = new Chart(ctx, config);
        
        // Setup clickable legend
        setupLegendToggles();
    }
    
    /**
     * Setup clickable legend items to toggle dataset visibility
     */
    function setupLegendToggles() {
        const legendItems = document.querySelectorAll('.legend-item');
        
        legendItems.forEach(item => {
            item.addEventListener('click', function() {
                const datasetIndex = parseInt(this.dataset.dataset);
                
                if (radarChart && radarChart.data.datasets[datasetIndex]) {
                    // Toggle visibility
                    const dataset = radarChart.data.datasets[datasetIndex];
                    dataset.hidden = !dataset.hidden;
                    
                    // Update legend item appearance
                    if (dataset.hidden) {
                        this.classList.add('disabled');
                    } else {
                        this.classList.remove('disabled');
                    }
                    
                    // Update chart
                    radarChart.update();
                }
            });
        });
    }

    /**
     * Calculate aggregated values for each functional category
     */
    function calculateCategoryAggregations(patientData, otherPatientsData) {
        const categoryNames = Object.keys(PROTEIN_CATEGORIES);
        const patientValues = [];
        const controlValues = [];
        const otherPatientsValues = [];

        categoryNames.forEach(categoryName => {
            const category = PROTEIN_CATEGORIES[categoryName];
            const categoryMarkers = category.markers;
            
            let patientSum = 0;
            let controlSum = 0;
            let otherPatientsSum = 0;
            let validMarkersCount = 0;

            categoryMarkers.forEach(markerName => {
                const patientMarker = patientData.find(p => p.marker_name === markerName);
                const otherPatientsMarker = otherPatientsData.find(p => p.marker_name === markerName);
                
                if (patientMarker && patientMarker.patient_value !== 'N/A' && patientMarker.control_range !== 'N/A') {
                    // Extract control range values
                    const controlRange = patientMarker.control_range.split(' - ');
                    if (controlRange.length === 2) {
                        const controlMin = parseFloat(controlRange[0]);
                        const controlMax = parseFloat(controlRange[1]);
                        const controlMean = (controlMin + controlMax) / 2;
                        
                        // Calculate as percentage of control mean
                        const patientValue = parseFloat(patientMarker.patient_value);
                        const otherPatientsValue = otherPatientsMarker && otherPatientsMarker.patient_value !== 'N/A' 
                            ? parseFloat(otherPatientsMarker.patient_value) 
                            : null;
                        
                        if (!isNaN(patientValue) && !isNaN(controlMean) && controlMean > 0) {
                            const patientPercent = (patientValue / controlMean) * 100;
                            const controlPercent = 100; // Controls are 100% by definition
                            
                            patientSum += patientPercent;
                            controlSum += controlPercent;
                            
                            // Add other patients percentage if available
                            if (otherPatientsValue !== null && !isNaN(otherPatientsValue)) {
                                const otherPatientsPercent = (otherPatientsValue / controlMean) * 100;
                                otherPatientsSum += otherPatientsPercent;
                            } else {
                                otherPatientsSum += 100; // Default to control level if no data
                            }
                            
                            validMarkersCount++;
                        }
                    }
                }
            });

            // Calculate average for this category
            if (validMarkersCount > 0) {
                patientValues.push(patientSum / validMarkersCount);
                controlValues.push(controlSum / validMarkersCount);
                otherPatientsValues.push(otherPatientsSum / validMarkersCount);
            } else {
                patientValues.push(0);
                controlValues.push(100);
                otherPatientsValues.push(100);
            }
        });

        return {
            patient: patientValues,
            controls: controlValues,
            otherPatients: otherPatientsValues
        };
    }
    
    // Expose global functions for external access
    window.ProteinAnalysis = {
        init: init,
        loadProteinData: loadProteinData,
        getCurrentData: () => currentProteinData
    };
    
    // Initialize when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }
    
})();
