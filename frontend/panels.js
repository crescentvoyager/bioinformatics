/**
 * Gene Panel Management System
 * Handles panel selection, data fetching, and PDF export functionality
 */

(function() {
    'use strict';
    
    // API Configuration
    const API_BASE = 'http://localhost:5001/api';
    
    // DOM Elements
    const elements = {
        // Tab navigation
        tabBtns: document.querySelectorAll('.tab-btn'),
        tabContents: document.querySelectorAll('.tab-content'),
        
        // Panel controls
        panelSelect: document.getElementById('panel-select'),
        panelPatientInput: document.getElementById('panel-patient-input'),
        exportPdfBtn: document.getElementById('export-pdf-btn'),
        
        // Panel description
        panelDescription: document.getElementById('panel-description'),
        descriptionTitle: document.getElementById('description-title'),
        descriptionText: document.getElementById('description-text'),
        
        // Panel preview
        panelPreview: document.getElementById('panel-preview'),
        panelStatus: document.getElementById('panel-status'),
        panelDataBody: document.getElementById('panel-data-body')
    };
    
    // State management
    let availablePanels = {};
    let availablePatients = [];
    let currentPanelId = null;
    let currentPatientId = null;
    let currentPanelData = null;
    
    // Panel descriptions
    const PANEL_DESCRIPTIONS = {
        "infection_panel_1": {
            title: "Infection Panel 1 Test",
            description: "This comprehensive infection panel contains 33 carefully selected genes that are key indicators of immune response and infection status. These genes are involved in innate and adaptive immunity, inflammatory responses, and pathogen recognition pathways. The panel is designed to provide insights into the patient's immune system activation and response to potential infections, making it valuable for clinical assessment and research applications."
        }
    };
    
    /**
     * Initialize the panel system
     */
    async function init() {
        try {
            setupTabNavigation();
            await loadAvailablePanels();
            await loadAvailablePatients();
            setupEventListeners();
            console.log('Panel system initialized successfully');
        } catch (error) {
            console.error('Failed to initialize panel system:', error);
            setPanelStatus('Failed to initialize panel system. Please check the backend connection.');
        }
    }
    
    /**
     * Setup tab navigation functionality
     */
    function setupTabNavigation() {
        elements.tabBtns.forEach(btn => {
            btn.addEventListener('click', () => {
                const targetTab = btn.getAttribute('onclick').match(/'([^']+)'/)[1];
                showTab(targetTab);
            });
        });
    }
    
    /**
     * Show specific tab and update navigation
     */
    window.showTab = function(tabId) {
        // Update tab buttons
        elements.tabBtns.forEach(btn => {
            btn.classList.remove('active');
            if (btn.getAttribute('onclick').includes(tabId)) {
                btn.classList.add('active');
            }
        });
        
        // Update tab content
        elements.tabContents.forEach(content => {
            content.classList.remove('active');
            if (content.id === tabId) {
                content.classList.add('active');
            }
        });
    };
    
    /**
     * Load available panels from the backend
     */
    async function loadAvailablePanels() {
        try {
            const response = await fetch(`${API_BASE}/panels`);
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            
            availablePanels = await response.json();
            populatePanelDropdown();
            
        } catch (error) {
            console.error('Error loading panels:', error);
            setPanelStatus('Failed to load available panels');
        }
    }
    
    /**
     * Load available patients from the backend
     */
    async function loadAvailablePatients() {
        try {
            const response = await fetch(`${API_BASE}/patients`);
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            
            const data = await response.json();
            availablePatients = data.patients || [];
            populatePatientDropdown();
            
        } catch (error) {
            console.error('Error loading patients:', error);
            setPanelStatus('Failed to load available patients');
        }
    }
    
    /**
     * Populate the panel selection dropdown
     */
    function populatePanelDropdown() {
        elements.panelSelect.innerHTML = '<option value="">Choose a panel...</option>';
        
        Object.entries(availablePanels).forEach(([panelId, panel]) => {
            const option = document.createElement('option');
            option.value = panelId;
            option.textContent = panel.name;
            elements.panelSelect.appendChild(option);
        });
    }
    
    /**
     * Populate the patient selection dropdown
     */
    function populatePatientDropdown() {
        elements.panelPatientInput.innerHTML = '<option value="">Select patient...</option>';
        
        availablePatients.forEach(patientId => {
            const option = document.createElement('option');
            option.value = patientId;
            option.textContent = patientId;
            elements.panelPatientInput.appendChild(option);
        });
    }
    
    /**
     * Setup event listeners for panel functionality
     */
    function setupEventListeners() {
        // Panel selection change
        elements.panelSelect.addEventListener('change', (e) => {
            const panelId = e.target.value;
            currentPanelId = panelId || null;
            
            // Show panel description
            if (panelId && PANEL_DESCRIPTIONS[panelId]) {
                const desc = PANEL_DESCRIPTIONS[panelId];
                elements.descriptionTitle.textContent = desc.title;
                elements.descriptionText.textContent = desc.description;
                elements.panelDescription.style.display = 'block';
            } else {
                elements.panelDescription.style.display = 'none';
            }
            
            checkAndLoadPanelData();
        });
        
        // Patient selection change
        elements.panelPatientInput.addEventListener('change', (e) => {
            currentPatientId = e.target.value;
            checkAndLoadPanelData();
        });
        
        // Export PDF button - open report page in new tab and auto-print
        elements.exportPdfBtn.addEventListener('click', () => {
            if (currentPanelId && currentPatientId) {
                setPanelStatus('Opening report for PDF export...');
                
                // Open the report page in a new tab (not popup)
                const reportUrl = `report.html?panel=${currentPanelId}&patient=${currentPatientId}&autoprint=true`;
                window.open(reportUrl, '_blank');
                
                // Set status back
                setTimeout(() => {
                    setPanelStatus('Report opened in new tab. Print dialog should appear automatically.');
                }, 1000);
            }
        });
    }
    
    /**
     * Check if both panel and patient are selected, then load data automatically
     */
    function checkAndLoadPanelData() {
        if (currentPanelId && currentPatientId) {
            generatePanelReport(currentPanelId, currentPatientId);
            elements.exportPdfBtn.style.display = 'block';
        } else {
            elements.panelPreview.style.display = 'none';
            elements.exportPdfBtn.style.display = 'none';
            if (!currentPanelId) {
                elements.panelDescription.style.display = 'none';
            }
        }
    }
    
    /**
     * Generate panel report and display preview
     */
    async function generatePanelReport(panelId, patientId) {
        setPanelStatus('Loading panel data...');
        
        try {
            const response = await fetch(`${API_BASE}/panel/${panelId}/patient/${patientId}`);
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            
            const data = await response.json();
            
            if (data.error) {
                throw new Error(data.error);
            }
            
            currentPanelData = data;
            displayPanelPreview(data);
            setPanelStatus(`Successfully loaded ${data.gene_count} genes`);
            
        } catch (error) {
            console.error('Error loading panel data:', error);
            setPanelStatus(`Error: ${error.message}`);
        }
    }
    
    /**
     * Display panel data preview table
     */
    function displayPanelPreview(data) {
        // Store panel data globally for PDF export
        currentPanelData = data;
        
        elements.panelPreview.style.display = 'block';
        
        // Clear existing table data
        elements.panelDataBody.innerHTML = '';
        
        // Populate table with gene data
        data.data.forEach(gene => {
            const row = document.createElement('tr');
            row.innerHTML = `
                <td>${gene.gene_id}</td>
                <td class="gene-name">
                    <a href="#" class="gene-link" data-gene="${gene.gene_name}" data-patient="${data.patient_id}">${gene.gene_name}</a>
                </td>
                <td class="numeric">${formatNumber(gene.patient_value)}</td>
                <td class="numeric">${gene.control_range || 'N/A'}</td>
                <td class="relative-position ${getRelativePositionClass(gene.relative_position)}">${gene.relative_position || 'N/A'}</td>
            `;
            elements.panelDataBody.appendChild(row);
        });
        
        // Add click handlers for gene links
        elements.panelDataBody.querySelectorAll('.gene-link').forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const geneName = e.target.dataset.gene;
                const patientId = e.target.dataset.patient;
                openGeneExplorer(geneName, patientId);
            });
        });
    }
    
    /**
     * Open gene explorer with specific gene and patient
     */
    function openGeneExplorer(geneName, patientId) {
        // Switch to Gene Analysis tab
        showTab('gene-analysis');
        
        // Set the gene input
        const geneInput = document.getElementById('gene-input');
        const patientInput = document.getElementById('patient-input');
        const plotBtn = document.getElementById('plot-btn');
        
        if (geneInput && patientInput && plotBtn) {
            geneInput.value = geneName;
            patientInput.value = patientId;
            
            // Trigger the plot
            plotBtn.click();
        }
    }
    
    /**
     * Set status message for panel operations
     */
    function setPanelStatus(message) {
        if (elements.panelStatus) {
            elements.panelStatus.textContent = message;
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
            return num.toFixed(3);
        } catch {
            return 'N/A';
        }
    }
    
    /**
     * Format scientific notation values for display
     */
    function formatScientific(value) {
        if (value === 'N/A' || value === null || value === undefined) {
            return 'N/A';
        }
        
        try {
            const num = parseFloat(value);
            if (isNaN(num)) return 'N/A';
            
            if (num < 0.001) {
                return num.toExponential(2);
            } else {
                return num.toFixed(4);
            }
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
    
    // Initialize when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }
    
})();
